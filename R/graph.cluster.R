###############################################################################
# Network Analysis Tasks
###############################################################################
clusterGraph.execute <- function (userCtx, sampleSetName, plotList, params) {
    sampleSet <- userCtx$sampleSets[[sampleSetName]]
    ctx     <- sampleSet$ctx
    dataset <- ctx$imputed
    distance <- dataset$distance
    meta     <- dataset$meta
    
    # Set up output folders
    dataFolder <- getOutFolder(ctx, sampleSetName, c("graph", "data"))

    # Get the cluster data
    clusterSetName  <- analysis.getParam ("cluster.clusterSet.name", params)	#; print(clusterSetName)
    clusterSetInfos <- cluster.getClustersSetFromContext (userCtx, sampleSetName, clusterSetName)
    
    # Process the clusters for different thresholds of min identity 
    setCount <- length(clusterSetInfos)
    for (idIdx in 1:setCount) {
        clusterSetInfo <- clusterSetInfos[[idIdx]]
        if (is.null(clusterSetInfo$clusters)) {		# No clusters present at this identity level
            next
        }
        
        # Build a graph from this distance matrix and clusters set
        graphInfo <- clusterGraph.buildGraph (distance, clusterSetInfo, clustersOnly=TRUE, params)
        
        # Save the graph data to file
        clusterGraph.saveGraphData (graphInfo, dataFolder)

        # Plot the clusters-only plot
        clusterGraph.makeClusterPlot (ctx, graphInfo, params)
        
        # TODO - Plot Lists not yet implemented
        #for (plotIdx in 1:length(plotList)) {
        #    plotDef <- plotList[[plotIdx]]
        #    plotName <- plotDef$name
        #    plotFilenameRoot  <- paste(plotFolder, "/graph-", sampleSetName, "-", minIdentityLabel, "-", plotName, sep="")
        #    clusterGraph.makeSinglePlot (gr, nodeData, clusterMembers, sampleSetName, minIdentityLabel, plotDef, layout, meta, plotFilenameRoot)            
        #}
    }
}
#
# 
#
clusterGraph.saveGraphData <- function (graphInfo, dataFolder) {
    clusterSetInfo <- graphInfo$clusterSetInfo
    filenameSuffix <- paste(clusterSetInfo$sampleSetName, 
                            clusterSetInfo$clusterSetName, 
                            getMinIdentityLabel(clusterSetInfo$minIdentity), sep="-")

    # Write the graph to a file
    igraph::write_graph(graphInfo$graph, 
                        paste0(dataFolder, "/", "graph", "-", filenameSuffix, ".graphml"), format="graphml")
    # Write cluster/sample distance data to file
    writeMatrix(graphInfo$distData, 
                        paste0(dataFolder, "/", "nodeDistance", "-", filenameSuffix, ".tab"))
    # Write cluster/sample information to file
    utils::write.table(graphInfo$nodeData, 
                       file=paste0(dataFolder, "/", "nodes", "-", filenameSuffix, ".tab"), 
                       quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}
#
###############################################################################
# Graph Plotting
###############################################################################
#
#
#
clusterGraph.makeClusterPlot <- function (ctx, graphInfo, params) {
    config         <- ctx$config
    gr             <- graphInfo$graph
    nodeData       <- graphInfo$nodeData
    clusterSetInfo <- graphInfo$clusterSetInfo
    sampleSetName  <- clusterSetInfo$sampleSetName
    clusterSetName <- clusterSetInfo$clusterSetName
    clustersOnly   <- graphInfo$clustersOnly

    minIdentityLabel <- getMinIdentityLabel (clusterSetInfo$minIdentity)	#; print(minIdentityLabel)
    print (paste("Graph Cluster Plot:", sampleSetName, clusterSetName, minIdentityLabel))
    
    # Construct a palette using the default palette, recycling it if there are too many clusters
    clusterNames <- nodeData$Node[which(nodeData$NodeType=="cluster")]		#; print(clusterNames)
    clusterPalette <- cluster.getClusterPalette (ctx, clusterNames)		#; print(clusterPalette)
    
    nodeData$color <- rep("gray", nrow(nodeData))				#; print(nodeData$color)
    nodeData[clusterNames,"color"] <- clusterPalette				#; print(nodeData$color)

    sampleCounts <- as.integer(nodeData$Count)
    nodeData$size <- sqrt(sampleCounts)  					#; print(nodeData$size)
    
    nodeData$label <- rep("", nrow(nodeData))					#; print(nodeData$label)
    nodeData[clusterNames,"label"] <- clusterNames
    
    # We now need to order all of the parameters according to the node ordering in the graph
    grNodeNames <- igraph::V(gr)$name						#; print(grNodeNames)
    nodeData <- nodeData[grNodeNames,]
    
    edgeIdentities <- as.numeric(igraph::E(gr)$Identity)			#; print(edgeIdentities)
    grEdgeColours <- clusterGraph.getEdgeWeightColours(edgeIdentities)		#; print(grEdgeColours)
    grEdgeWidths  <- clusterGraph.getEdgeWidth(edgeIdentities)
    #print(igraph::edge_attr(gr))

    # Set the layout for drawing the graph
    grLayout <- clusterGraph.getLayout(gr, params)
    
    # Draw the plot and save to file
    filenameSuffix <- paste(sampleSetName, clusterSetName, minIdentityLabel, sep="-")
    clGraphFilename <- paste("graph", filenameSuffix, "clusters", sep="-")
    plotFolder <- getOutFolder(ctx, sampleSetName, c("graph", "plots"))
    clGraphFile <- paste(plotFolder, clGraphFilename, sep="/")
    initializeGraphics (getGraphicsFilename (clGraphFile), widthInch=24, heightInch=24)
    plot(gr, layout=grLayout,
         edge.color=grEdgeColours, edge.width=grEdgeWidths, 
         vertex.shape="circle", vertex.label=nodeData$label, vertex.size=nodeData$size, vertex.color=nodeData$color)
    grDevices::dev.off()
}
#
#
#


#
#
#
clusterGraph.getEdgeWeightColours <- function(weights) {
    minWeight <- min(weights)
    maxWeight <- max(weights)
    steps <- 50
    stepSize <- (maxWeight - minWeight)/steps
    idx <- (1 + trunc((maxWeight - weights) / stepSize))
    grays <- grDevices::gray.colors(steps, start=0.0, end=0.95)
    grays <- c(grays,grays[steps])
    colours <- grays[idx]
    colours
}
#
#
#
clusterGraph.getEdgeWidth <- function(weights) {
    minWeight <- min(weights)
    maxWeight <- max(weights)
    prop <- ((weights - minWeight) / (maxWeight - minWeight))
    minWidth <- 1
    maxWidth <- 10
    widths <- minWidth + (prop * (maxWidth - minWidth))		#; widths <- rep(10,length(weights))
    widths
}
#
#
#
clusterGraph.getLayout <- function (gr, params) {
    layoutAlgo <- analysis.getParam ("graph.layoutAlgorithm", params)
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
###############################################################################
# Graph Construction
###############################################################################
clusterGraph.buildGraph <- function (distance, clusterSetInfo, clustersOnly, params) {

    # Analyze connectivity at a given identity threshold
    clusterData <- clusterSetInfo$clusters
    clusterMembers <- clusterSetInfo$members
    minIdentity <- clusterSetInfo$minIdentity				#; print(minIdentity)
    minIdentityLabel <- getMinIdentityLabel (minIdentity)		#; print(minIdentityLabel)
        
    # Use the cluster memberships to reduce the distance matrix by collapsing the cluster samples and replacing them with the clusters
    # The resulting table is a table of node distance rather than sample distance
    distData <- clusterGraph.buildNodeDistances (distance, clusterMembers)

    # If we're only doing the clusters, trim out the loose samples
    if (clustersOnly) {
        clusterIds <- clusterData$ClusterId
        distData <- distData[clusterIds,clusterIds]
    }

    # Create a node table, incorporating both newly found clusters and non-cluster samples
    nodeNames <- rownames(distData)
    nodeData <- clusterGraph.buildNodeData (nodeNames, clusterData)

    # Trim data to include only the nodes that are connected by edges; but make sure we preserve all clusters
    #minIdentity <- analysis.getParam ("graph.connectIdentityMin", params)
    minConnectIdentity <- minIdentity / 2
    #minConnectIdentity <- 0
    edgeData <- clusterGraph.getWeightedEdgeData (distData, minConnectIdentity, params)
    nodeNames <- unique(c(as.character(edgeData$Sample1),as.character(edgeData$Sample2)))

    # Create graph from connectivity data
    gr <- igraph::graph_from_data_frame(edgeData, directed=FALSE, vertices=nodeData[nodeNames,])
    grNodeNames <- names(igraph::V(gr))

    # Make sure we preserve all clusters, by adding any missing ones
    clusterNames <- nodeData$Node[which(nodeData$NodeType=="cluster")]
    addClusterNames <- clusterNames[!(clusterNames %in% grNodeNames)]
    if (length(addClusterNames) > 0) {
        addData <- nodeData[addClusterNames,]
        gr <- gr + igraph::vertices(addClusterNames, Count=addData$Count, NodeType=addData$NodeType)
    }						        #; print(names(igraph::V(gr)))

    graphInfo <- list (
        graph = gr,
        clusterSetInfo = clusterSetInfo,
        nodeData = nodeData,
        edgeData = edgeData,
        distData = distData,
        clustersOnly = clustersOnly
    )
    graphInfo
}
#
# Given a table of cluster memberships, remove the relevant samples from the distance matrix and 
# replace them with the clusters, using the mean distances as replacement
#
clusterGraph.buildNodeDistances <- function (distData, clusterMembers) {
    if (nrow(clusterMembers) == 0) {
        return (distData)
    }
    clNames <- unique(clusterMembers$Cluster)
    for (cldx in 1:length(clNames)) {
        clName <- clNames[cldx]
        clSamples <- clusterMembers$Sample[which(clusterMembers$Cluster == clName)]        #; print(clSamples)

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
clusterGraph.buildNodeData <- function (nodeNames, clusterData) {

    clusterIds    <- clusterData$ClusterId
    clusterCounts <- clusterData$Count
    
    nodeCounts <- rep(1, length(nodeNames))
    names(nodeCounts) <- nodeNames
    for (idx in 1:length(clusterIds)) {
        clusterId <- clusterIds[idx]
        nodeCounts[clusterId] <- clusterCounts[idx]
    }
    nodeTypes <- rep("sample", length(nodeNames))
    nodeTypes[which(nodeNames %in% clusterIds)] <- "cluster"

    # Create a new node list, including both samples and newly found clusters
    nodeData <- data.frame(Node=nodeNames, Count=nodeCounts, NodeType=nodeTypes)
    rownames(nodeData) <- nodeNames			    #print(nodeData)
    nodeData
}
#
# Get Weighted Graph Edge Data
# Return a table with edge definition, plus weights (based on distance via a transformation function)
#
clusterGraph.getWeightedEdgeData <- function (distData, minIdentity, params) {
    # Get a table of pairwise distance/identity values for all pairs of samples that meet the threshold
    pairData <- clusterGraph.getPairwiseIdentityData (distData, minIdentity, params)
    # Compute the link weights
    power <- analysis.getParam ("graph.weightPower", params)
    pairData$weight <- 100 * (pairData$Identity^power)
    pairData
}
#
# Get Unweighted Graph Edge Data
# Turn the distance matrix into a table of pairwise identity with four columns: "Sample1", "Sample2", "Distance", "Identity"
#
clusterGraph.getPairwiseIdentityData <- function (distData, minIdentity, params) {
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
#
#
#
###############################################################################
###############################################################################
# TODO - Graph plotting with color scheme- plot lists are to be re-introduced
###############################################################################
###############################################################################
#
#
#
clusterGraph.TODO.makeSinglePlot <- function (gr, nodeData, clusterMembers, sampleSetName, thresholdLabel, plotDef, grLayout, sampleMeta, plotFilenameRoot) {
    plotName <- plotDef$name
    print (paste("Graph Plot: ", sampleSetName, thresholdLabel, plotName))
    
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
    grEdgeColours <- clusterGraph.getEdgeWeightColours(edgeWeights)
    #grEdgeWidths <- clusterGraph.getEdgeWidth(edgeWeights)
    grEdgeWidths <- 1

    initializeGraphics (getGraphicsFilename (plotFilenameRoot), widthInch=24, heightInch=24)
    #print(grShapes)
    #print(markSizes)
    #print(grPieList)
    #print(grColours)
    plot(gr, layout=grLayout,
         edge.color=grEdgeColours, edge.width=grEdgeWidths, vertex.shape=grShapes, vertex.label=NA, vertex.size=markSizes, 
         vertex.pie=grPieList, vertex.pie.color=list(allColours), vertex.color=grColours)
    graphics::legend("bottomright", legendData$label, col="black", pt.bg=legendData$colour, pch=21)
    grDevices::dev.off()
}
#
###############################################################################
# Markers / igraph shapes
###############################################################################
clusterGraph.TODO.pchVertex <- function(coords, v = NULL, params) {
    v.pch    <- params("vertex", "pch")
    v.color  <- params("vertex", "color")
    v.lcolor <- params("vertex", "lcolor")
    v.lwd    <- params("vertex", "lwd")
    v.size   <- params("vertex", "size")
    graphics::points(x = coords[, 1], y = coords[, 2], pch=v.pch, bg=v.color, col=v.lcolor, lwd=v.lwd, cex=v.size)
}
igraph::add_shape("pch", plot=clusterGraph.TODO.pchVertex)
#
#add_shape("pch", clip=shapes("circle")$clip, plot=clusterGraph.TODO.pchVertex)
#
