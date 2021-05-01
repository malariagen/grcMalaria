###############################################################################
# Network Analysis Tasks
# Draw more edges, recompute subgraphs at every frequency level
###############################################################################
default.graph.connectIdentityMin <- 0.4
default.graph.weightFunction <- "identity"
default.graph.layoutAlgorithm <- "fr"

###############################################################################
# Network Analysis 
###############################################################################
graph.execute <- function (analysisContext, analysisName, params) {
    # Set up output folders
    dataFolder <- getOutFolder(analysisName, c("graph", "data"))
    
    # Get the connectivity thresholds, sorted from highest to lowest
    # Analyze connectivity at each thrshold level
    clusterIdentityLevels <- cluster.getIdentityLevels (params)
    for (idx in 1:length(clusterIdentityLevels)) {
    
        # Figure out the genetic identity cutoff, and keep only pairs above this threshold
        thresholdValue <- clusterIdentityLevels[idx]
        thresholdLabel <- cluster.getIdentityLevelLabel (thresholdValue)
        print(thresholdLabel)
        
        # Identify the larger clusters and get the cluster membership
        clustersData <- cluster.findbyIdentity (analysisContext, analysisName, thresholdValue, params)
        clusterMembers <- cluster.getMemberData (clustersData)        		#; print(head(clusterMembers))

        # If all samples are in a single cluster, then there's no point continuing, exit gracefully
        if (nrow(clustersData) == 1) {
            print(paste("All samples in a single subgroup at identity level",thresholdLabel,"- exiting"))
            break
        }

        # Use the cluster memberships to reduce the distance matrix by collapsing the cluster samples and replacing them with the clusters
        # The resulting table is a table of node distance rather than sample distance
        distData <- graph.buildNodeDistances (analysisContext$distance, clusterMembers)


        # Create a new node table, incorporating both newly found clusters and non-cluster samples
        nodeData <- graph.buildNodeData (analysisContext, clusterMembers)

        
        # Trim data to include only the nodes that are connected by edges; but make sure we preserve all clusters
        minIdentity <- analysis.getParam ("graph.connectIdentityMin", params, default.graph.connectIdentityMin)
        edgeData <- graph.getWeightedEdgeData (distData, minIdentity, params)
        nodeNames <- unique(c(as.character(edgeData$Sample1),as.character(edgeData$Sample2)))


        # Create graph from connectivity data
        gr <- igraph::graph_from_data_frame(edgeData, directed=FALSE, vertices=nodeData[nodeNames,])
        grNodeNames <- names(igraph::V(gr))

        # Make sure we preserve all clusters, by adding any missing ones
        clusterNames <- nodeData$Node[which(nodeData$NodeType=="cluster")]
        addClusterNames <- clusterNames[!(clusterNames %in% grNodeNames)]
        if (length(addClusterNames) > 0) {
            addData <- nodeData[addClusterNames,]
            gr <- gr + vertices(addClusterNames, Count=addData$Count, NodeType=addData$NodeType)
            grNodeNames <- names(igraph::V(gr))
        }						        #; print(length(grNodeNames)); print(grNodeNames)
        
        # Write the graph to a file
        grFile <- paste(dataFolder, "/graph-", analysisName, "-", thresholdLabel, ".graphml", sep="")
        igraph::write_graph(gr, grFile, format="graphml")
        # Write cluster membership to file
        if (nrow(clusterMembers) > 0) {
            clusterMembersFile <- paste(dataFolder, "/clusterMembers-", analysisName, "-", thresholdLabel, ".tab", sep="")
            writeSampleData(clusterMembers, clusterMembersFile)
        }
        # Write cluster/sample distance data to file
        distDataFile <- paste(dataFolder, "/clusterDistance-", analysisName, "-", thresholdLabel, ".tab", sep="")
        writeMatrix(distData, distDataFile)
        # Write cluster/sample information to file
        nodeDataFile <- paste(dataFolder, "/nodes-", analysisName, "-", thresholdLabel, ".tab", sep="")
        write.table(nodeData, file=nodeDataFile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    }
}
#
# Given a table of cluster memberships, remove the relevant samples from the distance matrix and 
# replace them with the clusters, using the mean distances as replacement
#
graph.buildNodeDistances <- function (distData, clusterMembers) {
    if (nrow(clusterMembers) == 0) {
        return (distData)
    }
    clNames <- unique(clusterMembers$Cluster)
    for (cldx in 1:length(clNames)) {
        clName <- clNames[cldx]
        clSamples <- clusterMembers$Sample[which(clusterMembers$Cluster == clName)]

        # Get a matrix of distance between cluster samples (cols) and other samples (rows)
        clDist <- distData[,clSamples]

        # Compute the mean distance of cluster samples vs all others, and then stick it a first column and row
        clMeans <- rowMeans(clDist)
        sNames <- c(clName, colnames(distData))
        distData <- rbind(clMeans, distData)
        clMeans <- c(0.0, clMeans)
        distData <- cbind(clMeans, distData)
        rownames(distData) <- sNames
        colnames(distData) <- sNames
       
        # Remove cluster samples from distance matrix
        distData <- distData[!(rownames(distData) %in% clSamples),]
        distData <- distData[,!(colnames(distData) %in% clSamples)]
    }
    distData
}
#
# Build a table of nodes, to distinguish clusters from individual nodes, and giving them a size
#
graph.buildNodeData <- function (analysisContext, clusterMembers) {
    # Create a new node list, including both samples and newly found clusters
    # Start with a list of nodes that are single samples, no clusters
    sampleNames <- rownames(analysisContext$meta)
    nodeData <- data.frame(Node=sampleNames, Count=rep(1,length(sampleNames)), NodeType=rep("sample",length(sampleNames)))
    rownames(nodeData) <- sampleNames
    # If the previous pass produced clusters, we need to create new nodes for them, and remove the samples therein
    if (nrow(clusterMembers) > 0) {
        # Order clusters from the previous pass and get their names
        clusterCounts <- table(clusterMembers$Cluster)
        clusterCounts <- clusterCounts[order(clusterCounts, decreasing=TRUE)]
        clusterNames <- names(clusterCounts)
        #print(clusterNames)
        # For each cluster, remove the samples from the node list, and add the cluster node instead
        for (subIdx in 1:length(clusterNames)) {
            clusterName <- clusterNames[subIdx]
            clusterCount <- clusterCounts[subIdx]
            # Get the list of samples names for this cluster, and remove them from the node list, then add the cluster to replace them
            nodeData <- rbind(nodeData, c(clusterName, clusterCount, "cluster"))
            clusterSamples <- clusterMembers$Sample[which(clusterMembers$Cluster==clusterName)]
            nodeData <- nodeData[!(row.names(nodeData) %in% clusterSamples),]
        }
    }
    rownames(nodeData) <- nodeData$Node
    #print(nodeData)
    nodeData
}
#
###############################################################################
# Network Graph Plotting
###############################################################################
#
graph.executePlots <- function (analysisContext, analysisName, plotList, params) {
    # Set up output folders 
    dataFolder <- getOutFolder(analysisName, c("graph", "data"))
    
    # Get the connectivity thresholds, sorted from highest to lowest
    # Analyze connectivity at each thrshold level
    clusterIdentityLevels <- cluster.getIdentityLevels (params)
    for (idx in 1:length(clusterIdentityLevels)) {
        thresholdValue <- clusterIdentityLevels[idx]
        thresholdLabel <- cluster.getIdentityLevelLabel (thresholdValue)

        # Get the graph node table, incorporating both newly found clusters and non-cluster samples
        nodeDataFile <- paste(dataFolder, "/nodes-", analysisName, "-", thresholdLabel, ".tab", sep="")
        if (!file.exists(nodeDataFile)) {
            print(paste("No data found for identity level",thresholdLabel,"- skipping"))
            next
        }
        nodeData <- read.table(nodeDataFile, as.is=TRUE, header=TRUE, sep="\t")
        rownames(nodeData) <- nodeData$Node
        
        # Get the sample cluster membership table
        subMembersFile <- paste(dataFolder, "/clusterMembers-", analysisName, "-", thresholdLabel, ".tab", sep="")
        clusterMembers <- readSampleData(subMembersFile)

        # Read the graph from file
        grFile <- paste(dataFolder, "/graph-", analysisName, "-", thresholdLabel, ".graphml", sep="")
	gr <- igraph::read_graph(grFile, format="graphml")
	layout <- graph.getLayout(gr, params)
        
	# Perform plots for this onnectivity threshold
        plotFolder <- getOutFolder(analysisName, c("graph", "plots"))	#; print(plotList)

        for (plotIdx in 1:length(plotList)) {
            plotDef <- plotList[[plotIdx]]
            plotName <- plotDef$name
            plotFilenameRoot  <- paste(plotFolder, "/graph-", analysisName, "-", thresholdLabel, "-", plotName, sep="")
            graph.makeSinglePlot (gr, nodeData, clusterMembers, analysisName, thresholdValue, plotDef, layout, analysisContext$meta, plotFilenameRoot)            
        }
        
        # Also plot the clusters-only plot
        #graph.makeClusterPlot (gr, nodeData, distData, analysisName, thresholdValue, plotFolder, params)
        plotFilenameRoot  <- paste(plotFolder, "/graph-", analysisName, "-", thresholdLabel, "--clusters", sep="")
        graph.makeClusterPlot (gr, nodeData, analysisName, thresholdValue, layout, plotFilenameRoot)            
    }
}
#
graph.getLayout <- function (gr, params) {
    layoutAlgo <- analysis.getParam ("graph.layoutAlgorithm", params, default.graph.layoutAlgorithm)
    if (layoutAlgo == "fr") {
        layout <- igraph::layout_with_fr(gr)
    } else if (layoutAlgo == "kk") {
        layout <- igraph::layout_with_kk(gr)
    } else if (layoutAlgo == "lgl") {
        layout <- igraph::layout_with_lgl(gr)
    } else if (layoutAlgo == "drl") {
        layout <- igraph::layout_with_drl(gr)
    } else if (layoutAlgo == "mds") {
        layout <- igraph::layout_with_mds(gr)
    } else if (layoutAlgo == "graphopt") {
        layout <- igraph::layout_with_graphopt(gr)
    }
    layout
}
#
# Get Weighted Graph Edge Data
# Return a table with edge definition, plus weights (based on distance via a transformation function)
#
graph.getWeightedEdgeData <- function (distData, minIdentity, params) {
    # Get a table of pairwise distance/identity values for all pairs of samples that meet the threshold
    pairData <- cluster.getPairwiseIdentityData (distData, minIdentity, params)
    # Compute the link weights
    fnName <- analysis.getParam ("graph.weightFunction", params, default.graph.weightFunction)
    # Default function
    if (fnName == "identity") {
        identityToWeight <- function (identities) identities
    } else if (fnName == "identitySquared") {
        identityToWeight <- function (identities) (identities * identities)
    }
    pairData$weight <- identityToWeight(pairData$Identity)
    pairData
}
#
#
#
graph.makeSinglePlot <- function (gr, nodeData, clusterMembers, analysisName, thresholdValue, plotDef, grLayout, sampleMeta, plotFilenameRoot) {
    thresholdLabel <- cluster.getIdentityLevelLabel (thresholdValue)
    plotName <- plotDef$name
    print (paste("Graph Plot: ", analysisName, thresholdLabel, plotName))
    
    # Set up the graphical attributes for rendering each individual sample
    #print(colnames(sampleMeta))
    attrList <- applyGraphicalAttributes(sampleMeta, plotDef$render)
    sampleMeta <- attrList$meta
    legendData <- attrList$legend

    # The following code needs to construct the following things: 
    #     1) A vector of vertex marker types, which specifies "pie" for vertices representing clusters, and "circle" for vertices 
    #        representing individual non-cluster samples.
    #     2) For the cluster vertices, which are represented as pie charts, a list of vectors of sample count by colour.
    #        The table must contain entries for all vertices, including individual non-cluster samples, which will be set to a dummy default 
    #        vector of counts (could be anything other than all zeroes)
    #     3) For individual non-cluster samples, a vector of colours (e.g. like in PCA).
    #        It will also contain entries for the pie chart vertices, whih will be set to default colour (can be any valid colour)

    # Get the graphical parameters for each sample- only colour used here
    nodeColours <- sampleMeta$plot__colour
    names(nodeColours) <- rownames(sampleMeta)

    allColours <- unique(nodeColours)
    # Put the default colour at the front
    allColours <- c("gray", allColours[allColours != "gray"])
    numColours <- length(allColours)
    
    # Make a dummy table of default colour for the individual samples
    sampleNames <- rownames(sampleMeta)
    nodePieTable <- matrix(c(1, rep(0, (numColours-1))), byrow=TRUE,
                 nrow=length(sampleNames), ncol=length(allColours), dimnames=list(sampleNames,allColours))
   
    # Create a vector of shapes
    nodeShapes <- rep("circle",length(sampleNames))
    names(nodeShapes) <- sampleNames

    # Create a vector of sample counts
    nodeSampleCounts <- rep(1,length(sampleNames))
    names(nodeSampleCounts) <- sampleNames
    
    if (!is.null(clusterMembers)) {
        clusterNames <- unique(clusterMembers$Cluster)
        # Make a table of colour counts for all the clusters
        subPieTable <- matrix(0, nrow=length(clusterNames), ncol=length(allColours), dimnames=list(clusterNames,allColours))
        subSampleCounts <- c()
        for (subIdx in 1 : length(clusterNames)) {
            clusterName <- clusterNames[subIdx]							#; print(clusterName)
            subSamples <- clusterMembers$Sample[which(clusterMembers$Cluster==clusterName)]
            subSampleCounts <- c(subSampleCounts,length(subSamples))
            subColours <- nodeColours[subSamples]
            subColCounts <- table(subColours)
            subPieTable[subIdx, names(subColCounts)] <- subColCounts
        }
        names(subSampleCounts) <- clusterNames
        # Merge with the pie colour table for the sample
        nodePieTable <- rbind(subPieTable, nodePieTable)
        
        # Create a dummy circle colour vector for the clusters
        subColours <- rep("gray", length(clusterNames))
        names(subColours) <- clusterNames
        # Merge the colour vectors
        nodeColours <- c(nodeColours, subColours)

        # Create a vector of shapes
        subShapes <- rep("pie",length(clusterNames))
        names(subShapes) <- clusterNames
        # Merge the shape vectors
        nodeShapes <- c(nodeShapes, subShapes)
    
        # Merge the sample counts vectors
        nodeSampleCounts <- c(subSampleCounts, nodeSampleCounts)
    }
    
    # We now need to order all of the parameters according to the node ordering in the graph
    grNodeNames <- igraph::V(gr)$name
    #print(grNodeNames)
    grShapes  <- nodeShapes[grNodeNames]
    grColours <- nodeColours[grNodeNames]

    grPieTable <- t(nodePieTable[grNodeNames,])
    #print(grPieTable)
    grPieList  <- split(grPieTable, rep(1:ncol(grPieTable), each=nrow(grPieTable)))

    grSampleCounts <- nodeSampleCounts[grNodeNames]
    #print(nodeSampleCounts)
    #print(grSampleCounts)
    markSizes <- sqrt(grSampleCounts)
    #markSizes <- normalizeMarkerSizes(markSizes)
    
    edgeWeights <- as.numeric(igraph::E(gr)$weight)
    grEdgeColours <- graph.getEdgeWeightColours(edgeWeights)
    #grEdgeWidths <- graph.getEdgeWidth(edgeWeights)
    grEdgeWidths <- 1

    initializeGraphics (getGraphicsFilename (plotFilenameRoot), widthInch=24, heightInch=24)
    #print(grShapes)
    #print(markSizes)
    #print(grPieList)
    #print(grColours)
    plot(gr, layout=grLayout,
         edge.color=grEdgeColours, edge.width=grEdgeWidths, vertex.shape=grShapes, vertex.label=NA, vertex.size=markSizes, 
         vertex.pie=grPieList, vertex.pie.color=list(allColours), vertex.color=grColours)
    legend("bottomright", legendData$label, col="black", pt.bg=legendData$colour, pch=21)
    dev.off()
}
#
#
#
#
graph.makeClusterPlot <- function (gr, nodeData, analysisName, thresholdValue, grLayout, plotFilenameRoot) {
    thresholdLabel <- cluster.getIdentityLevelLabel (thresholdValue)
    print (paste("Graph Cluster Plot: ", analysisName, thresholdLabel))
    
    clusterPalette <- cluster.getPalette (params)
    clusterNames <- nodeData$Node[which(nodeData$NodeType=="cluster")]
    numClusters <- min(length(clusterPalette),length(clusterNames))		#; print(numClusters)
    clusterPalette <- clusterPalette[1:numClusters]
    clusterNames <- clusterNames[1:numClusters]					#; print(clusterNames)
    names(clusterPalette) <- clusterNames					#; print(clusterPalette)
    										#; print(nrow(nodeData)); print(rownames(nodeData))
    nodeData$color <- rep("gray", nrow(nodeData))				#; print(nodeData$color)
    nodeData[clusterNames,"color"] <- clusterPalette				#; print(nodeData$color)

    sampleCounts <- as.integer(nodeData$Count)
    nodeData$size <- sqrt(sampleCounts)  					#; print(nodeData$size)
    
    nodeData$label <- rep("", nrow(nodeData))					#; print(nodeData$label)
    nodeData[clusterNames,"label"] <- clusterNames
    
    # We now need to order all of the parameters according to the node ordering in the graph
    grNodeNames <- igraph::V(gr)$name							#; print(grNodeNames)
    nodeData <- nodeData[grNodeNames,]
    
    edgeWeights <- as.numeric(igraph::E(gr)$weight)
    grEdgeColours <- graph.getEdgeWeightColours(edgeWeights)
    #grEdgeWidths <- graph.getEdgeWidth(edgeWeights)
    grEdgeWidths <- 1

    initializeGraphics (getGraphicsFilename (plotFilenameRoot), widthInch=24, heightInch=24)
    plot(gr, layout=grLayout,
         edge.color=grEdgeColours, edge.width=grEdgeWidths, 
         vertex.shape="circle", vertex.label=nodeData$label, vertex.size=nodeData$size, vertex.color=nodeData$color)
    #legend("bottomright", clusterNames, col="black", pt.bg=clusterPalette, pch=21)
    dev.off()
}
#
graph.makeClusterPlot.TO_BE_REMOVED <- function (gr, nodeData, distData, analysisName, thresholdValue, plotFolder, params) {
    thresholdLabel <- cluster.getIdentityLevelLabel (thresholdValue)
    print (paste("Graph Cluster Plot: ", analysisName, thresholdLabel))
    
    # Remove all objects from the graph other than the top clusters (as many as can be represented with the cluster palette)
    clusterPalette <- cluster.getPalette (params)
    maxClusters <- length(clusterPalette)
    clusterNames <- nodeData$Node[which(nodeData$NodeType=="cluster")]
    clusterNames <- clusterNames[order(clusterNames)]				#; print(clusterNames)
    if (length(clusterNames) > maxClusters) {
        clusterNames <- clusterNames[1:maxClusters]				#; print(clusterNames)
    }
    nodeData <- nodeData[clusterNames,]
    clusterCount <- nrow(nodeData)						#; print(clusterCount)
    sampleCounts <- as.integer(nodeData$Count)
    nodeData$size <- 5+sqrt(sampleCounts)  					#; print(nodeData$size)
    nodeData$color <- clusterPalette[1:clusterCount]				#; print(nodeData$color)
    nodeData$label <- nodeData$Node

    distData <- distData[clusterNames, clusterNames]
    edgeData <- graph.getWeightedEdgeData (distData, 0.0, params)
    edgeData$color <- graph.getEdgeWeightColours (edgeData$weight)
    edgeData$width <- graph.getEdgeWidth (edgeData$weight)
    edgeData$weight <- graph.normalizeWeights (edgeData$weight)

    # Create graph from connectivity data
    gr <- igraph::graph_from_data_frame(edgeData, directed=FALSE, vertices=nodeData)
    layout <- graph.getLayout(gr, params)
    
    # Write the graph to a file
    #grFile <- paste(plotFolder, "/graph-clusters-", analysisName, "-", thresholdLabel, ".graphml", sep="")
    #write_graph(gr, grFile, format="graphml")

    plotFilenameRoot  <- paste(plotFolder, "/graph-clusters-", analysisName, "-", thresholdLabel, sep="")
    initializeGraphics (getGraphicsFilename (plotFilenameRoot), widthInch=24, heightInch=24)
    plot(gr, layout=layout, vertex.shape="circle")
    #legend("bottomright", legendData$label, col="black", pt.bg=legendData$colour, pch=21)
    dev.off()
}
#
graph.getEdgeWeightColours <- function(weights) {
    minWeight <- min(weights)
    maxWeight <- max(weights)
    steps <- 50
    stepSize <- (maxWeight - minWeight)/steps
    idx <- (1 + trunc((maxWeight - weights) / stepSize))
    grays <- gray.colors(steps, start=0.0, end=0.95)		#; grays <- rainbow(steps)
    grays <- c(grays,grays[steps])
    colours <- grays[idx]					#; colours <- rep("red",length(weights))
    colours
}
#
graph.getEdgeWidth <- function(weights) {
    minWeight <- min(weights)
    maxWeight <- max(weights)
    prop <- ((weights - minWeight) / (maxWeight - minWeight))
    minWidth <- 0.5
    maxWidth <- 10
    widths <- minWidth + (prop * (maxWidth - minWidth))		#; widths <- rep(10,length(weights))
    widths
}
#
graph.normalizeWeights <- function(weights) {
    minWeight <- min(weights)
    maxWeight <- max(weights)
    prop <- ((weights - minWeight) / (maxWeight - minWeight))
    prop
}
#
###############################################################################
# Markers / igraph shapes
###############################################################################
graph.pchVertex <- function(coords, v = NULL, params) {
    v.pch    <- params("vertex", "pch")
    v.color  <- params("vertex", "color")
    v.lcolor <- params("vertex", "lcolor")
    v.lwd    <- params("vertex", "lwd")
    v.size   <- params("vertex", "size")
    points(x = coords[, 1], y = coords[, 2], pch=v.pch, bg=v.color, col=v.lcolor, lwd=v.lwd, cex=v.size)
}
igraph::add_shape("pch", plot=graph.pchVertex)
#
#add_shape("pch", clip=shapes("circle")$clip, plot=graph.pchVertex)
#
# Test code
#g <- graph.ring(10, dir = FALSE)
#plot(g, vertex.shape="pch", vertex.label=NA, vertex.pch=25, vertex.color=c("red","yellow"), vertex.lcolor="green", vertex.size=3, vertex.lwd=2.5)
#
