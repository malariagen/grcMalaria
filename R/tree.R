#
# NJ Tree plot settings - TODO
#
tree.colouredBranches <- TRUE
tree.showLeaves       <- TRUE
tree.nonLeafColour    <- "darkgrey"
tree.branchThickness  <- 3
tree.popLabels        <- FALSE
tree.sampleLabels     <- FALSE
#
###############################################################################
# Tree Analysis Routines
# Computes and plots trees (e.g. NJ trees)
# Multiple visualization (e.g. different colouring schemes) can be performed in a single run.
###############################################################################
#
tree.execute <- function(userCtx, sampleSetName, method, params) {
    sampleSet <- userCtx$sampleSets[[sampleSetName]]
    ctx <- sampleSet$ctx
    #
    # Get the output folders
    #
    dataFolder      <- getOutFolder(ctx$config, sampleSetName, c(method, "data"))
    plotsRootFolder <- getOutFolder(ctx$config, sampleSetName, c(method, "plots"))
    #
    # Get the metadata and distance and genotypes
    #
    useImputed <- param.getParam ("phylo.impute", params)	#; print(useImputed)
    datasetName <- ifelse (useImputed, "imputed", "filtered")
    dataset <- ctx[[datasetName]]
    sampleMeta <- dataset$meta
    distData  <- distance.retrieveDistanceMatrix (ctx, datasetName)
    #
    # Build a Tree (currently only NJ trees supported)
    #
    tree <- NULL
    if (method == "njt") {
        tree <- ape::nj(as.matrix(distData));
    } else  {
        stop(paste("Invalid tree type specified:",method))
    }
    #
    # Write tree to file, in case we want to sue it later
    #
    treeFilename  <- paste0(dataFolder, "/", method, "-", sampleSetName, ".newick")
    ape::write.tree(tree, file=treeFilename)
    #
    # Order the metadata according to the tree results
    #
    sampleNames <- tree$tip.label
    sampleMeta <- sampleMeta[sampleNames,]
    #
    # Execute the plots
    #
    plotList <- param.getParam ("plot.plotList", params)		#; print(plotList)
    for (plotIdx in 1:length(plotList)) {
        plotDef <- plotList[[plotIdx]]
        plotDefName <- plotDef$name;
        plotName <- paste(sampleSetName, plotDefName, method, sep="-")
        print (paste("Tree Plot: ", plotName))
        #
        # Set up the graphical attributes for rendering
        #
        ga <- graphics.getGraphicalAttributes(userCtx, sampleMeta, plotDef$attributes)
        gaData <- ga$sampleAttrData
        legendData <- ga$legendData
        #
        # Merge graphics attributes and metadata
        #
        plotMetadata <- cbind(sampleMeta, gaData)
        sampleNames <- rownames(plotMetadata)
        sampleCount <- length(sampleNames)
        #
        # Write out the plot data (facilitate debug)
        #
        plotMetadataFilename  <- paste0(dataFolder, "/", "plotData-", plotName)
        writeSampleData(plotMetadata, plotMetadataFilename)
        #
        # Do the plots
        #
        plotsFolder  <- getSubFolder (plotsRootFolder, plotDefName)
        #
        # Graphics elements
        #
        tipLabels <- NULL
        tree.showLeafSymbol <- FALSE
        tree.showLeafText <- FALSE
        if (tree.showLeaves) {
            if (tree.sampleLabels) {
                tipLabels <- sampleNames
                tipSizeFactor <- 0.5
                tree.showLeafText <- TRUE
            } else {
                tipLabels <- rep("",sampleCount)
                tipSizeFactor <- 1
                tree.showLeafSymbol <- TRUE
            }
        }
        tree$tip.label <- tipLabels

        tipColours <- plotMetadata$plot__colour
        tipPch     <- plotMetadata$plot__pch
        tipSize    <- plotMetadata$plot__size
        tipLcolours<- plotMetadata$plot__lcolour
        tipLwd     <- plotMetadata$plot__lwd
        #
        # Colour the tree edges
        #
        edgeColours <- tree.colourEdges (tree, plotMetadata)
        #
        # Plot tree
        #
        treeName  <- paste("tree", plotName, sep="-")
        graphicFilenameRoot  <- paste(plotsFolder, treeName, sep="/")
        initializeGraphics (getGraphicsFilename (graphicFilenameRoot), params)
        graphics::par(mar=c(3.1, 3.1, 3.1, 3.1))
        treePlot <- NULL
        treePlot <- ape::plot.phylo(tree, type="unrooted", root.edge=FALSE, edge.color=edgeColours, edge.width=tree.branchThickness, 
                                    show.tip.label=tree.showLeafText, tip.color=tipColours, label.offset=2, cex=tipSizeFactor, font=2)
    
        sizes <- 2.0 * as.numeric(tipSize)
        lineWidths <- 2.0 * as.numeric(tipLwd)
        if (tree.showLeafSymbol) {
            ape::tiplabels(NULL, pch=as.numeric(tipPch), cex=sizes, 
                       frame="none", col=tipLcolours, bg=tipColours, lwd=lineWidths)
        }
        grDevices::dev.off()
        #
        # Create a separate legend
        #
        graphicFilenameRoot  <- paste0(plotsFolder, "/", treeName, "-legend")
        initializeGraphics (getGraphicsFilename (graphicFilenameRoot), params)
        graphics::plot.new()
        graphics::par(mar=c(0,0,0,0))
        graphics::legend("topleft", ncol=1, inset=0.05, cex=1.0, 
                         legendData$label, pt.bg=legendData$colour, col=legendData$lcolour, 
                         pch=as.numeric(legendData$pch), pt.lwd=as.numeric(legendData$lwd), 
                         pt.cex=as.numeric(legendData$size))

        grDevices::dev.off()        
    }
}


tree.colourEdges <- function(tree, sampleMeta) {
    # print(tree$edge)
    # tree$edge is a two-column matrix of mode numeric where each row represents an edge of the
    # tree; the nodes and the tips are symbolized with numbers; the tips are numbered 1, 2, . . . , 
    # and the nodes are numbered after the tips. For each row, the first column gives the ancestor

    numSamples <- nrow(sampleMeta)
    numEdges <- dim(tree$edge)[1]
    edgeColours <- rep(tree.nonLeafColour, numEdges);
    if (tree.colouredBranches) {
        parentEdgeIds <- c()
        for (leafNodeId in 1:nrow(sampleMeta)) {
            leafNodeIdx <- which(tree$edge[,2] == leafNodeId)
            edgeColours[leafNodeIdx] <- sampleMeta$plot__colour[leafNodeId]
            parentEdgeIds <- c(parentEdgeIds, tree$edge[leafNodeIdx,1])    # Record the node
        }
        #print(parentEdgeIds)
        edgeColours <- tree.colourInnerEdges (tree, edgeColours, parentEdgeIds)
    }
    return (edgeColours)
}

tree.colourInnerEdges <- function(tree, edgeColours, innerEdgeIds) {
    if (length(innerEdgeIds) == 0) {
        return (edgeColours);
    }
    innerEdgeIds <- unique(innerEdgeIds)
    parentEdgeIds <- c()
    for (idx in 1:length(innerEdgeIds)) {
        innerEdgeId <- innerEdgeIds[idx]
        innerEdgeIdx <- which(tree$edge[,2] == innerEdgeId)
        #print(paste(idx, innerEdgeId))
        childIndexes <- which(tree$edge[,1] == innerEdgeId)
        childColours <- unique(edgeColours[childIndexes])
        #print(childColours)
        if (length(childColours) == 1) { # All the same colour
            edgeColours[innerEdgeIdx] <- childColours[1]
            parentEdgeIds <- c(parentEdgeIds, tree$edge[innerEdgeIdx,1])
        }
    }
    edgeColours <- tree.colourInnerEdges (tree, edgeColours, parentEdgeIds)  # Recurse to deal with parent edges
    return (edgeColours)
}

