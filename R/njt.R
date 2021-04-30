#
# NJ Tree plot settings
#
nj.colouredBranches <- TRUE
nj.showLeaves       <- TRUE
nj.nonLeafColour    <- "darkgrey"
nj.branchThickness  <- 3
nj.popLabels        <- FALSE
nj.sampleLabels     <- FALSE
#
###############################################################################
# NJT Analysis Routines
# Computes and plots NJ trees.
# Multiple analyses, based on different subsets of samples, can be performed in a single run.
###############################################################################
#
njt.execute <- function(analysisContext, analysisName) {
    dataFolder <- getOutFolder(analysisName, c("njt", "data"))
    distData  <- analysisContext$distance

    # Build a NJ Tree
    tree <- nj(as.matrix(distData));
    treeFilename  <- paste(dataFolder, "/njt-", analysisName, ".newick", sep="")
    write.tree(tree, file=treeFilename)
}

njt.executePlots <- function(analysisContext, analysisName, plotList) {
    dataFolder <- getOutFolder(analysisName, c("njt", "data"))
    sampleMeta <- analysisContext$meta

    # Retrieve a NJ Tree
    treeFilename  <- paste(dataFolder, "/njt-", analysisName, ".newick", sep="")
    tree <- read.tree(file=treeFilename)
    sampleNames <- tree$tip.label
    sampleMeta <- sampleMeta[sampleNames,]
  
    # Execute the plots
    #print(plotList)
    for (pIdx in 1:length(plotList)) {
        plotDef <- plotList[[pIdx]]
        njt.plotNjt (analysisName, tree, sampleMeta, plotDef)
    }
}

njt.plotNjt <- function(analysisName, tree, sampleMeta, plotDef) {
    plotName <- plotDef$name;
    print (paste("Plot: ",plotName))
  
    # Set up the graphical attributes for rendering
    dataList <- applyGraphicalAttributes(sampleMeta, plotDef$render)
    plotMetadata <- dataList$meta
    legendData <- dataList$legend
    sampleCount <- nrow(plotMetadata)
  
    # Do the plots
    plotsFolder <- getOutFolder(analysisName, c("njt", "plots"))
    
    # Graphics elements
    tipLabels <- NULL
    nj.showLeafSymbol <- FALSE
    nj.showLeafText <- FALSE
    if (nj.showLeaves) {
        if (nj.sampleLabels) {
            tipLabels <- plotMetadata[,sampleMetaSampleColumnName]
            tipSizeFactor <- 0.5
            nj.showLeafText <- TRUE
        } else {
            tipLabels <- rep("",sampleCount)
            tipSizeFactor <- 1
            nj.showLeafSymbol <- TRUE
        }
    }
    tree$tip.label <- tipLabels
    
    tipColours <- plotMetadata$plot__colour
    tipPch     <- plotMetadata$plot__pch
    tipSize    <- plotMetadata$plot__size
    tipLcolours<- plotMetadata$plot__lcolour
    tipLwd     <- plotMetadata$plot__lwd

    # Colour the tree edges
    edgeColours <- njt.colourEdges (tree, plotMetadata)
    
    # Plot tree
    treeName  <- paste("njt", analysisName, plotName, sep="-")
    graphicFilenameRoot  <- paste(plotsFolder, treeName, sep="/")
    initializeGraphics (getGraphicsFilename (graphicFilenameRoot), widthInch=24, heightInch=24)
    par(mar=c(3.1, 3.1, 3.1, 3.1))
    treePlot <- NULL
    treePlot <- plot.phylo(tree, type="unrooted", root.edge=FALSE, edge.color=edgeColours, edge.width=nj.branchThickness, 
                show.tip.label=nj.showLeafText, tip.color=tipColours, label.offset=2, cex=tipSizeFactor, font=2)
    #print(treePlot)
    
    sizes <- 2.0 * as.numeric(tipSize)
    lineWidths <- 2.0 * as.numeric(tipLwd)
    if (nj.showLeafSymbol) {
        tiplabels(NULL, pch=as.numeric(tipPch), cex=sizes, 
                  frame="none", col=tipLcolours, bg=tipColours, lwd=lineWidths)
    }
    dev.off()
}

njt.colourEdges <- function(tree, sampleMeta) {
    # print(tree$edge)
    # tree$edge is a two-column matrix of mode numeric where each row represents an edge of the
    # tree; the nodes and the tips are symbolized with numbers; the tips are numbered 1, 2, . . . , 
    # and the nodes are numbered after the tips. For each row, the first column gives the ancestor

    numSamples <- nrow(sampleMeta)
    numEdges <- dim(tree$edge)[1]
    edgeColours <- rep(nj.nonLeafColour, numEdges);
    if (nj.colouredBranches) {
        parentEdgeIds <- c()
        for (leafNodeId in 1:nrow(sampleMeta)) {
            leafNodeIdx <- which(tree$edge[,2] == leafNodeId)
            edgeColours[leafNodeIdx] <- sampleMeta$plot__colour[leafNodeId]
            parentEdgeIds <- c(parentEdgeIds, tree$edge[leafNodeIdx,1])    # Record the node
        }
        #print(parentEdgeIds)
        edgeColours <- njt.colourInnerEdges (tree, edgeColours, parentEdgeIds)
    }
    return (edgeColours)
}

njt.colourInnerEdges <- function(tree, edgeColours, innerEdgeIds) {
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
    edgeColours <- njt.colourInnerEdges (tree, edgeColours, parentEdgeIds)  # Recurse to deal with parent edges
    return (edgeColours)
}

