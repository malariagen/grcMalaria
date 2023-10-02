###############################################################################
# Cluster Analysis 
###############################################################################
#
# Main entry point- load stored clusters, or identify them from scratch if they are not stored
#
cluster.findClusters <- function (userCtx, sampleSetName, params) {
    if (!sampleSetName %in% names(userCtx$sampleSets)) {
        stop(paste("Sample set not initialized:", sampleSetName))
    }
    sampleSet <- userCtx$sampleSets[[sampleSetName]]
    ctx       <- sampleSet$ctx
    config    <- ctx$config

    clusterSetName    <- param.getParam ("cluster.clusterSet.name", params)	#; print(clusterSetName)
    minIdentityLevels <- param.getParam ("cluster.identity.min", params) 	#; print(minIdentityLevels)
    imputeBarcodes    <- param.getParam ("cluster.impute", params)		#; print(imputeBarcodes)
    method            <- param.getParam ("cluster.method", params)		#; print(method)
    print(paste("Clustering Method:",method))
    
    clusterSetInfos <- list()
    for (idIdx in 1:length(minIdentityLevels)) {
        minIdentity <- minIdentityLevels[idIdx]
        minIdentityLabel <- getMinIdentityLabel (minIdentity)

        # Determine what clustering approach must be used.
        # TODO - At the moment, only graph-based clustering methods are implemented, but this may be extended later.
        clustersData <- NULL
        if (method %in% cluster.graphMethods) {
            clustersData <- cluster.findClustersFromGraph (ctx, clusterSetName, method, minIdentity, imputeBarcodes, params)
        } else {
            stop (paste("Invalid method specified in parameter 'cluster.method':", method))
        }
    
        # Sort clusters by descending size, and label them in order, appending the sequence number to the cluster set name
        clustersData <- clustersData[rev(order(clustersData$Count)),]		#; print(head(clustersData))
        clusterSetName <- param.getParam ("cluster.clusterSet.name", params)	#; print(clusterSetName)
        clustersData$ClusterId <- paste(clusterSetName, 
                                        formatC(seq(1,nrow(clustersData)), width=3, format="d", flag="0"), sep="-")
        rownames(clustersData) <- clustersData$ClusterId			#; print(head(clustersData))
        clustersData <- as.data.frame(clustersData[,c("ClusterId","Count","SampleList")])

        # Write out result files
        dataFileSuffix <- paste("", sampleSetName, clusterSetName, minIdentityLabel, sep="-")
        dataFolder  <- getOutFolder(ctx$config, sampleSetName, 
                                    c("cluster", "data", clusterSetName, minIdentityLabel))
    
        # Write out cluster definitions
        clustersDataFile <- paste0(dataFolder, "/clusters", dataFileSuffix, ".tab")
        utils::write.table(clustersData, file=clustersDataFile, quote=TRUE, sep="\t", row.names=FALSE, col.names=TRUE)

        # Write out sample/cluster association
        memberData <- cluster.getMemberData (clustersData)
        memberDataFile <- paste0(dataFolder, "/clusterMembers", dataFileSuffix)
        writeSampleData(memberData, memberDataFile)

        # Write out cluster stats
        statsData <- cluster.getClusterStats (ctx, clustersData, memberData)
        statsDataFile <- paste0(dataFolder, "/clusterStats", dataFileSuffix)
        writeSampleData(statsData, statsDataFile)
        
        # Keep the cluster info for storing in the context
        clusterSetInfo <- list(clusterSetName=clusterSetName,
                               sampleSetName=sampleSetName,
                               minIdentity=minIdentity,
                               clusters=clustersData,	# dataframe
                               members=memberData, 	# dataframe
                               stats=statsData)		# dataframe
        clusterSetInfos[[minIdentityLabel]] <- clusterSetInfo
    }
    
    # Reference the cluster data from the context
    sampleSet$clusters[[clusterSetName]] <- clusterSetInfos	#; print(names(sampleSet$clusters))
}
#
# Retrieve the cluster data from the context
#
cluster.getClustersSetFromContext <- function(userCtx, sampleSetName, clusterSetName) {
    sampleSet <- userCtx$sampleSets[[sampleSetName]]		#; print(names(userCtx)); print(names(userCtx$sampleSets))
    if (is.null(sampleSet)) {
        stop(paste("Sample set not found:", sampleSetName))        
    }								#; print(names(sampleSet));    
    clusterSetInfos <- sampleSet$clusters[[clusterSetName]]
    if (is.null(clusterSetInfos)) {
       stop(paste("Invalid cluster set specified:", clusterSetName))        
    }
    clusterSetInfos
}
#
#
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
#
#
cluster.getClusterPalette <- function(ctx, clusterIds) {
    # Get the name of all the clusters and sort them 
    clusterIds <- sort(clusterIds)
    # Put "Other" at the end if it is there
    if ("Other" %in% clusterIds) {
        clusterIds <- c(clusterIds[which(clusterIds != "Other")], "Other")
    }									#; print(clusterIds)
    # Construct a palette using the default palette, recycling it if there are too many clusters
    # and adding white as the last colour for class "Other" (samples not belionging to a cluster)
    colPalette <- graphics.getColourPalette (ctx)
    clusterPalette <- rep_len(colPalette, length.out=(length(clusterIds)-1))
    clusterPalette <- c(clusterPalette,"white")
    names(clusterPalette) <- clusterIds					#; print(clusterPalette)
    clusterPalette
}
#
# #######################################################################################
#
# Descriptive data about the clusters (e.g. prevalence of mutations, etc.)
#
cluster.getClusterStats <- function(ctx, clustersData, clusterMembers) {
    config <- ctx$config

    clNames <- rownames(clustersData)			#; print(head(clustersData))	#; print(clNames)
    
    # Create a table of stats data
    statsData <- NULL
    statsCols <- NULL
    for (clIdx in 1 : length(clNames)) {
        clName <- clNames[clIdx]
        
        # Get the metadata for this cluster
        clSampleNames <- clusterMembers$Sample[which(clusterMembers$Cluster==clName)]
        clSampleMeta <- ctx$unfiltered$meta[clSampleNames,]
        
        # Get the sample count first
        statsCols <- "Count"
        statValues <- length(clSampleNames)

        # Get the drug resistance prevalences for this cluster
        resistanceColumns  <- config$cluster.stats.drugs
        if (!is.null(resistanceColumns)) { 
            clPrevalence <- meta.getResistancePrevalence (ctx, clSampleMeta, resistanceColumns)
            clPrevalence <- format(as.numeric(clPrevalence), digits=2, nsmall=2)
            statsCols <- c(statsCols, resistanceColumns)
            statValues <- c(statValues, clPrevalence)
        }

        # Get the counts for this cluster
        countColumns  <- config$cluster.stats.alleleCounts
        if (!is.null(countColumns)) { 
            clCounts <- meta.getValueCounts (clSampleMeta, countColumns)
            statsCols <- c(statsCols, countColumns)
            statValues <- c(statValues, clCounts)
        }

        # Get the mutation prevalences for this cluster
        mutationColumns  <- config$cluster.stats.mutations
        if (!is.null(mutationColumns)) { 
            clPrevalence <- meta.getMutationPrevalence (ctx, clSampleMeta, mutationColumns)
            clPrevalence <- format(as.numeric(clPrevalence), digits=2, nsmall=2)
            statsCols <- c(statsCols, mutationColumns)
            statValues <- c(statValues, clPrevalence)
        }

        # Stick the joined row data to the Stats table
        names(statValues) <- statsCols
        statsData <- rbind(statsData, statValues)
    }
    statsData <- data.frame(statsData)
    colnames(statsData) <- statsCols
    rownames(statsData) <- clNames
    statsData
}
#
cluster.getClusterStatsText <- function(clusterStats, clustersName) {
    statsNames  <- colnames(clusterStats)
    statsValues <- clusterStats[clustersName,]
    clusterInfoTextLines <- paste(statsNames, statsValues, sep=": ")
    clusterInfoText <- paste(clusterInfoTextLines, collapse = "\n")		#; print(clusterInfoText)
    clusterInfoText
}
#
###############################################################################
# Graph-based Clustering
###############################################################################
cluster.graphCommunityMethods <- c("louvain")
cluster.graphMethods          <- c("allNeighbours", cluster.graphCommunityMethods)

cluster.findClustersFromGraph <- function (ctx, clusterSetName, method, minIdentity, imputeBarcodes, params) {
    datasetName <- ifelse (imputeBarcodes, "imputed", "filtered")
    config <- ctx$config

    # Get a table of pairwise distance/identity values for all pairs of samples that meet the threshold
    distData  <- distance.retrieveDistanceMatrix (ctx, datasetName)		#; print(nrow(distData))

    edgeData <- clusterGraph.getPairwiseIdentityData (distData, minIdentity, params)
    edgeData$weight <- (edgeData$Identity * edgeData$Identity)			#; print(nrow(edgeData))
    
    nodeNames <- unique(c(as.character(edgeData$Sample1),as.character(edgeData$Sample2)))
    nodeCount <- length(nodeNames)						#; print (paste(length(nodeNames),length(unique(edgeData$Sample1)),length(unique(edgeData$Sample2))))
    nodeData  <- data.frame(NodeName=nodeNames, Count=rep(1,nodeCount), NodeType=rep("sample",nodeCount))
    gr <- igraph::graph_from_data_frame(edgeData, directed=FALSE, vertices=nodeData)	#; print(paste("processClusters",minIdentity))

    # Perform clustering from the graph, identifying all clusters of sufficient size
    if (method == "allNeighbours") {
        clustersList <- cluster.findAllNeighbourClusters (gr, params)
    } else if (method %in% cluster.graphCommunityMethods) {
        clustersList <- cluster.findGraphCommunities (gr, method, params)
    }
    
    # Check if we actually found any clusters
    clCount <- length(clustersList)						#; print(clCount)
    if (clCount == 0) {
        print(paste("No clusters of desired minimum side were found at identity threshold", minIdentity))
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
    clustersData <- data.frame(Count=sampleCounts, SampleList=sampleLists, 
                               stringsAsFactors=FALSE)				#; print(head(clustersData))
    clustersData
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
    minCount <- param.getParam ("cluster.minSize", params)
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
    
    minCount <- param.getParam ("cluster.minSize", params)
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
