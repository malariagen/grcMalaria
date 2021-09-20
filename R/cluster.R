
###############################################################################
# Cluster Analysis 
###############################################################################
#
# Main entry point- load stored clusters, or identify them from scratch if they are not stored
#
cluster.findbyIdentity <- function (ctx, datasetName, analysisName, thresholdValue, params) {
    dataset <- ctx[[datasetName]]    			#; print(thresholdValue)
    clustersData <- cluster.getClustersData (ctx, analysisName, thresholdValue)
    if (!is.null(clustersData)) {
        return (clustersData)
    }
    thresholdLabel <- cluster.getIdentityLevelLabel (thresholdValue)

    # Create graph from connectivity data
    distData <- dataset$distance			#; print(nrow(distData))
    edgeData <- graph.getWeightedEdgeData (distData, thresholdValue, params)	#; print(nrow(edgeData))
    nodeNames <- unique(c(as.character(edgeData$Sample1),as.character(edgeData$Sample2)))
    nodeCount <- length(nodeNames)			#; print (paste(length(nodeNames),length(unique(edgeData$Sample1)),length(unique(edgeData$Sample2))))
    nodeData <- data.frame(NodeName=nodeNames, Count=rep(1,nodeCount), NodeType=rep("sample",nodeCount))
    gr <- igraph::graph_from_data_frame(edgeData, directed=FALSE, vertices=nodeData)	#; print(paste("processClusters",thresholdValue))

    # Identify all clusters of sufficient size
    communityMethods <- c("louvain")
    method <- analysis.getParam ("cluster.method", params)
    if (method == "allNeighbours") {
        clustersList <- cluster.findAllNeighbourClusters (gr, params)
    } else if (method %in% communityMethods) {
        clustersList <- cluster.findGraphCommunities (gr, method, params)
    } else {
        stop (paste("Invalid method specified in parameter 'cluster.method':", method))
    }

    clCount <- length(clustersList)					#; print(clCount)
    if (clCount == 0) {
        print(paste("No clusters of desired minimum side were found at identity threshold", thresholdValue))
        return (NULL)
    }
    
    # Turn the list of clusters into a data frame    
    sampleCounts <- integer(clCount)
    sampleLists <- character(clCount)
    for (clIdx in 1 : clCount) {
        clSampleNames <- clustersList[[clIdx]]
        sampleCounts[clIdx] <- length(clSampleNames)
        sampleLists[clIdx] <- paste(clSampleNames, collapse=",")
    }
    clustersData <- data.frame(Count=sampleCounts, SampleList=sampleLists, stringsAsFactors=FALSE)	#; print(head(clustersData))
    
    # Sort clusters by descending size, and label them in order
    clustersData <- clustersData[rev(order(clustersData$Count)),]	#; print(head(clustersData))
    clustersData$ClusterId <- paste0("HG",formatC(seq(1,clCount), width=3, format="d", flag="0"))
    rownames(clustersData) <- clustersData$ClusterId			#; print(head(clustersData))

    # Write out sample/cluster association
    clusterMembers <- cluster.getMemberData (clustersData)
    clMembersFile <- cluster.getClustersDataFile(ctx, analysisName, "clusterMembers", thresholdValue)
    writeSampleData(clusterMembers, clMembersFile)

    # Write out cluster definitions
    clustersData <- cluster.estimateClusterStats (ctx, clustersData, clusterMembers, dataset$meta)
    clustersDataFile <- cluster.getClustersDataFile(ctx, analysisName, "clusters", thresholdValue)
    utils::write.table(clustersData, file=clustersDataFile, quote=TRUE, sep="\t", row.names=FALSE, col.names=TRUE)

    clustersData
}
#
# Get the connectivity thresholds, sorted from highest to lowest
#
cluster.getIdentityLevels <- function (params) {
    identityLevels <- analysis.getParam ("cluster.identity.thresholds", params)
    identityLevels <- identityLevels[order(identityLevels, decreasing=TRUE)]	#; print(identityLevels)
    identityLevels
}
#
cluster.getIdentityLevelLabel <- function (thresholdValue) {
    label <- paste("ge", format(thresholdValue, digits=2, nsmall=2), sep="")
    label
}
#
# Retrieve clusters from file storage
#
cluster.getClustersData <- function(ctx, analysisName, thresholdValue) {
    clustersDataFile <- cluster.getClustersDataFile(ctx, analysisName, "clusters", thresholdValue, createFolder=FALSE)
    if (!file.exists(clustersDataFile)) {
        return(NULL)
    }
    clustersData <- utils::read.table(clustersDataFile, as.is=TRUE, header=TRUE, sep="\t", quote="\"", check.names=FALSE)
    rownames(clustersData) <- clustersData$ClusterId
    clustersData
}
#
cluster.getClustersDataFile <- function(ctx, analysisName, filePrefix, thresholdValue, createFolder=TRUE) {
    thresholdLabel <- cluster.getIdentityLevelLabel (thresholdValue)
    dataFolder  <- getOutFolder(ctx, analysisName, c("cluster", "data", thresholdLabel), createFolder)
    clustersDataFile <- paste(dataFolder, "/", filePrefix, "-", analysisName, "-", thresholdLabel, ".tab", sep="")
    clustersDataFile
}
#
cluster.getMemberData <- function(clustersData) {
    cIds <- c()
    sIds <- c()
    for (clIdx in 1 : nrow(clustersData)) {
         # Get members of the cluster and label them
         clusterName <- clustersData$ClusterId[clIdx]
         clSampleNames <- unlist(strsplit(clustersData$SampleList[clIdx], split=","))
         sIds <- c(sIds, clSampleNames)
         cIds <- c(cIds, rep(clusterName, length(clSampleNames)))
    }
    #print(sIds)
    #print(cIds)
    clusterMembers <- data.frame(Sample=sIds, Cluster=cIds, stringsAsFactors=FALSE)
    rownames(clusterMembers) <- sIds
    clusterMembers
}
#
# #######################################################################################
# Descriptive data about the clusters (e.g. prevalence of mutations, etc.)
#
cluster.estimateClusterStats <- function(ctx, clustersData, clusterMembers, sampleMeta) {
    config <- ctx$config
    clNames <- rownames(clustersData)			#; print(head(clustersData))	#; print(clNames)
    
    # Create a table of stats data
    statsData <- NULL
    
    # Get the prevalence/counts columns to be reported- if nothing, don't bother doing calculations
    prevCols  <- c()
    if(!is.null(config$cluster.prevalenceColumns)) { 
        prevCols  <- config$cluster.prevalenceColumns
    }

    countCols  <- c()
    if(!is.null(config$cluster.countColumns)) {
        countCols <- config$cluster.countColumns
    }
    cNames <- c(prevCols, countCols)

    if (length(cNames) == 0) {
        return(clustersData)
    }

    for (clIdx in 1 : length(clNames)) {
        clName <- clNames[clIdx]
        clSampleNames <- clusterMembers$Sample[which(clusterMembers$Cluster==clName)]
        clSampleMeta <- sampleMeta[clSampleNames,]
        
        statValues <- c()
        # Get the prevalence of drug resistances for this cluster- result is a vector of numeric
        if (length(prevCols) > 0) {
            clPrevalence <- meta.getResistancePrevalence (ctx, clSampleMeta, prevCols)
            clPrevalence <- round(as.numeric(clPrevalence), digits=2)
            statValues <- c(statValues, clPrevalence)
	}

        # Get the allele counts for this cluster- result is a vector of integers
        if (length(countCols) > 0) {			#; print (countCols)
            clCounts <- meta.getValueCounts (clSampleMeta, countCols)
            statValues <- c(statValues, clCounts)
        }
        # Stick the joined row data to the Stats table
        statsData <- rbind(statsData, statValues)
    }
    colnames(statsData) <- cNames
    clustersData <- cbind(clustersData, statsData)
    clustersData
}
#
cluster.getClusterStatsText <- function(ctx, clusterName, clustersData) {
    config <- ctx$config
    # Get the prevalence/counts columns to be reported- if nothing, don't bother doing calculations
    prevCols  <- c()
    if(!is.null(config$cluster.prevalenceColumns)) {
        prevCols  <- config$cluster.prevalenceColumns
    }
    countCols  <- c()
    if(!is.null(config$cluster.countColumns)) { 
        countCols <- config$cluster.countColumns
    }
    inclCols <- c("Count", prevCols, countCols)
    inclVals <- clustersData[clusterName,inclCols]
    statTextLines <- paste(inclCols, inclVals, sep=": ")
    statText <- paste(unlist(statTextLines), collapse = "\n")
    statText
}
#
# #######################################################################################
#
# Clustering Method: identification by traversing All Neighbours.
#
# To find each cluster, we pick a node, and find all neighbouts recursively.
#
cluster.findAllNeighbourClusters <- function (gr, params) {
    # Make a working copy of the graph and count the nodes
    wg <- igraph::induced_subgraph(gr, igraph::V(gr))	#; plot(wg)
    nodeCnt <- length(igraph::V(wg))    		#; print(paste("nodeCnt",nodeCnt))

    # Identify all clusters with >= minCount samples
    minCount <- analysis.getParam ("cluster.identity.minCount", params)
    clList <- list()
    while (nodeCnt > 0) {				#; print(nodeCnt)
        #
        # Get the fist node in the graph, and build a cluster by traversing to all neighbours
        #
        firstNode <- igraph::V(wg)[1]
        clNodeIdxs <- cluster.findAllNeighbourConnectedNodes (wg, firstNode, as.integer(firstNode))
        						#; print(paste("new cluster - nodes:",paste(clNodeIdxs,collapse=",")))
        clSampleCount <- length(clNodeIdxs)		#; print(clSampleCount)
        if (clSampleCount >= minCount) {
            clSampleNames <- names(igraph::V(wg))[clNodeIdxs]
            clIdx <- length(clList)+1
            clList[[clIdx]] <- clSampleNames
        }
        # Remove the samples from the cluster we found
        keepSamples <- seq(1, nodeCnt)[-clNodeIdxs]
        wg <- igraph::induced_subgraph(wg, keepSamples, impl="copy_and_delete")
        nodeCnt <- length(igraph::V(wg))
    }
    clList
}

cluster.findAllNeighbourConnectedNodes <- function (gr, curr, nodes) {
    ns <- igraph::neighbors(gr, curr)
    if (length(ns) > 0) {
        foundNew <- FALSE
        for (nsIdx in 1 : length(ns)) {
            n <- ns[nsIdx]
            nl <- as.integer(n)
            if (!(nl %in% nodes)) {
                foundNew <- TRUE
                nodes <- c (nodes, nl)
                nodes <- cluster.findAllNeighbourConnectedNodes (gr, n, nodes)
            }
        }
    }
    nodes
}
#
# #######################################################################################
#
# Clustering Method: Community Analysis.
#
cluster.findGraphCommunities <- function (gr, method, params) {

    partition <- igraph::cluster_louvain(gr, weights=igraph::E(gr)$weight) 
    nodeComms <- partition$membership
    commIds <- sort(as.integer(unique(nodeComms)))
    sampleNames <- names(igraph::V(gr))
    #names(nodeComms) <- sampleNames
    
    minCount <- analysis.getParam ("cluster.identity.minCount", params)
    clList <- list()
    for (clIdx in 1:length(commIds)) {
        commId <- commIds[clIdx]
        clNodeIdxs <- which(nodeComms == commId)
        clSampleCount <- length(clNodeIdxs)		#; print(clSampleCount)
        if (clSampleCount < minCount) {
            next
        }
        clSampleNames <- sampleNames[clNodeIdxs]
        clIdx <- length(clList)+1
        clList[[clIdx]] <- clSampleNames
    }
    clList
}
