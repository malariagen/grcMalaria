###############################################################################
# 
###############################################################################
default.cluster.identity.thresholds <- c(1.0)
default.cluster.identity.minCount <- 5

###############################################################################
# Cluster Analysis 
###############################################################################
#
cluster.findbyIdentity <- function (analysisContext, analysisName, thresholdValue, params) {
    clustersData <- cluster.getClustersData (analysisName, thresholdValue)
    if (!is.null(clustersData)) {
        return (clustersData)
    }

    thresholdLabel <- cluster.getIdentityLevelLabel (thresholdValue)
    minCount <- analysis.getParam ("cluster.identity.minCount", params, default.cluster.identity.minCount)

    # Create graph from connectivity data
    distData <- analysisContext$distance		#; print(nrow(distData))
    edgeData <- cluster.getPairwiseIdentityData (distData, thresholdValue, params)	#; print(nrow(edgeData))
    nodeNames <- unique(c(as.character(edgeData$Sample1),as.character(edgeData$Sample2)))
    nodeCount <- length(nodeNames)			#; print (paste(length(nodeNames),length(unique(edgeData$Sample1)),length(unique(edgeData$Sample2))))
    nodeData <- data.frame(NodeName=nodeNames, Count=rep(1,nodeCount), NodeType=rep("sample",nodeCount))
    gr <- igraph::graph_from_data_frame(edgeData, directed=FALSE, vertices=nodeData)	#; print(paste("processClusters",thresholdValue))

    # Make a working copy of the graph and count the nodes
    wg <- igraph::induced_subgraph(gr, igraph::V(gr))	#; plot(wg)
    nodeCnt <- length(igraph::V(wg))    		#; print(paste("nodeCnt",nodeCnt))

    # Identify all clusters with >= minCount samples
    clustersData <- NULL
    while (nodeCnt > 0) {				#; print(nodeCnt)
        clNodeIdxs <- cluster.findSampleNodes(wg)	#; print(paste("new cluster - nodes:",paste(clNodeIdxs,collapse=",")))
        clSampleCount <- length(clNodeIdxs)		#; print(clSampleCount)
        if (clSampleCount >= minCount) {
            clSampleNames <- names(igraph::V(wg))[clNodeIdxs]
            clSampleNameList <- paste(clSampleNames, collapse=",")
            clustersData <- rbind(clustersData, c(clSampleCount, clSampleNameList))
        }
        # Remove the samples from the cluster we found
        keepSamples <- seq(1, nodeCnt)[-clNodeIdxs]
        wg <- igraph::induced_subgraph(wg, keepSamples, impl="copy_and_delete")
        nodeCnt <- length(igraph::V(wg))
    }									#; print(head(clustersData))
    if (is.null(clustersData)) {
        print(paste("No clusters of desired minimum side were found at identity threshold", thresholdValue))
        return (NULL)
    }
    clustersData <- data.frame(clustersData, stringsAsFactors=FALSE)
    colnames(clustersData) <- c("Count", "SampleList")			#; print(head(clustersData))
    clustersData$Count      <- as.integer(clustersData$Count)		#; print(head(clustersData))
    clustersData$SampleList <- as.character(clustersData$SampleList)	#; print(head(clustersData))
    
    # Sort clusters by descending size
    clustersData <- clustersData[rev(order(clustersData$Count)),]	#; print(head(clustersData))

    clCount <- nrow(clustersData)					#; print(clCount)
    clustersData$ClusterId <- paste("HG",formatC(seq(1,clCount), width=3, format="d", flag="0"),sep="")
    rownames(clustersData) <- clustersData$ClusterId			#; print(head(clustersData))

    # Write out sample/cluster association
    clusterMembers <- cluster.getMemberData (clustersData)
    clMembersFile <- cluster.getClustersDataFile("clusterMembers", analysisName, thresholdValue)
    writeSampleData(clusterMembers, clMembersFile)

    # Write out cluster definitions
    clustersData <- cluster.estimateClusterStats (clustersData, clusterMembers, analysisContext$meta) 
    clustersDataFile <- cluster.getClustersDataFile("clusters", analysisName, thresholdValue)
    utils::write.table(clustersData, file=clustersDataFile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

    clustersData
}

cluster.estimateClusterStats <- function(clustersData, clusterMembers, sampleMeta) {
    clNames <- rownames(clustersData)			#; print(head(clustersData)); print(clNames)
    
    # Create a table of stats data
    statsData <- NULL
    
    # Get the prevalence/counts columns to be reported- if nothing, don't bother doing calculations
    prevCols  <- c()
    if(!is.null(cluster.prevalenceColumns)) { prevCols  <- cluster.prevalenceColumns}
    countCols  <- c()
    if(!is.null(cluster.countColumns)) { countCols <- cluster.countColumns}
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
            clPrevalence <- meta.getResistancePrevalence (clSampleMeta, prevCols)
            clPrevalence <- round(as.numeric(clPrevalence), digits=2)
            statValues <- c(statValues, clPrevalence)
	}
        # Get the allele counts for this cluster- result is a vector of integers
        if (length(countCols) > 0) {
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

cluster.getClusterStatsText <- function(clusterName, clustersData) {
    # Get the prevalence/counts columns to be reported- if nothing, don't bother doing calculations
    prevCols  <- c()
    if(!is.null(cluster.prevalenceColumns)) { prevCols  <- cluster.prevalenceColumns}
    countCols  <- c()
    if(!is.null(cluster.countColumns)) { countCols <- cluster.countColumns}

    inclCols <- c("Count", prevCols, countCols)
    inclVals <- clustersData[clusterName,inclCols]
    statTextLines <- paste(inclCols, inclVals, sep=": ")
    statText <- paste(unlist(statTextLines), collapse = "\n")
    statText
}

cluster.getClustersData <- function(analysisName, thresholdValue) {
    clustersDataFile <- cluster.getClustersDataFile("clusters", analysisName, thresholdValue, createFolder=FALSE)
    if (!file.exists(clustersDataFile)) {
        return(NULL)
    }
    clustersData <- utils::read.table(clustersDataFile, as.is=TRUE, header=TRUE, quote="", sep="\t", check.names=FALSE)
    rownames(clustersData) <- clustersData$ClusterId
    clustersData
}

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

cluster.getClustersDataFile <- function(filePrefix, analysisName, thresholdValue, createFolder=TRUE) {
    thresholdLabel <- cluster.getIdentityLevelLabel (thresholdValue)
    dataFolder  <- getOutFolder(analysisName, c("cluster", "data", thresholdLabel), createFolder)
    clustersDataFile <- paste(dataFolder, "/", filePrefix, "-", analysisName, "-", thresholdLabel, ".tab", sep="")
    clustersDataFile
}

cluster.findSampleNodes <- function (gr) {
    first <- igraph::V(gr)[1]
    nodes <- c(as.integer(first))
    #print(nodes)
    nodes <- cluster.findNodes (gr, first, nodes)
    nodes
}

cluster.findNodes <- function (gr, curr, nodes) {
    ns <- igraph::neighbors(gr, curr)
    if (length(ns) > 0) {
        foundNew <- FALSE
        for (nsIdx in 1 : length(ns)) {
            n <- ns[nsIdx]
            nl <- as.integer(n)
            if (!(nl %in% nodes)) {
                foundNew <- TRUE
                nodes <- c (nodes, nl)
                nodes <- cluster.findNodes (gr, n, nodes)
            }
        }
    }
    nodes
}

#
# Get the connectivity thresholds, sorted from highest to lowest
#
cluster.getIdentityLevels <- function (params) {
    identityLevels <- analysis.getParam ("cluster.identity.thresholds", params, default.cluster.identity.thresholds)
    identityLevels <- identityLevels[order(identityLevels, decreasing=TRUE)]	#; print(identityLevels)
    identityLevels
}

cluster.getIdentityLevelLabel <- function (thresholdValue) {
    label <- paste("ge", format(thresholdValue, digits=2, nsmall=2), sep="")
    label
}

#
# Get Unweighted Graph Edge Data
# Turn the distance matrix into a table of pairwise identity with three columns: "Sample1", "Sample2", "Distance", "Identity"
#
cluster.getPairwiseIdentityData <- function (distData, minIdentity, params) {
    mat <- as.matrix(distData)
    mat[lower.tri(mat,diag=TRUE)] = NA
    pairData <- as.data.frame(as.table(mat))
    pairData <- stats::na.omit(pairData) # remove NA
    colnames(pairData) <- c("Sample1", "Sample2", "Distance")
    
    # Convert genetic distance to barcode identity
    pairData$Identity <- 1.0 - pairData$Distance
    
    # Keep only pairs above  threshold
    pairData <- pairData[which(pairData$Identity >= minIdentity),]
    pairData
}
