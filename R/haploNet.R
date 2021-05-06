###############################################################################
# Haplotype Network Analysis.
# Computes and plots PCA, using a pairwise distance matrix as an input.
# Multiple analyses, based on different subsets of samples, can be performed in a single run.
#
################################################################################
#
haploNet.execute <- function(ctx, datasetName, analysisName, plotList, params) {
    dataset <- ctx[[datasetName]]
    # Set up output folders
    dataFolder <- getOutFolder(ctx, analysisName, c("haploNet", "data"))
    plotsFolder <- getOutFolder(ctx, analysisName, c("haploNet", "plots"))
    
    sampleMeta <- dataset$meta
    sampleNames <- rownames(sampleMeta)
    barcodeData <- dataset$barcodes[sampleNames,]
    
    #barcodeData <- barcodeData[,1:20]
    barcodeSeqs <- ape::as.DNAbin(as.matrix(barcodeData))
    h <- pegas::haplotype(barcodeSeqs)
    
    # Keep only haplotypes that are common (i.e. in a number of parasites)
    minHaploCount <- analysis.getParam ("haploNet.minHaploCount", params)
    h <- subset(h, minfreq=minHaploCount)
    
    net <- pegas::haploNet(h)
    
    # Execute the plots
    #print(plotList)
    for (pIdx in 1:length(plotList)) {
        plotDef <- plotList[[pIdx]]
        haploNet.plotNet (net, h, analysisName, sampleMeta, plotDef, plotsFolder)
    }
}
#maxRecords <- 100

haploNet.plotNet <- function(net, h, analysisName, sampleMeta, plotDef, plotsFolder) {
    
    plotName <- plotDef$name
    print (paste("Haplotype Network Plot:",plotName))
        
    # Set up the graphical attributes for rendering each individual sample
    attrList <- applyGraphicalAttributes(sampleMeta, plotDef$render)
    sampleMeta <- attrList$meta
    legendData <- attrList$legend

    # Get the graphical parameters for each sample- only colour used here
    sampleColours <- sampleMeta$plot__colour
    names(sampleColours) <- rownames(sampleMeta)
    
    # Create a colour table for the pie charts
    pieColTable<-with(
        utils::stack(stats::setNames(attr(h, "index"), rownames(h))), 
        table(hap=ind, pop=sampleColours[values])
    )
    pieCols <- colnames(pieColTable)

    # Do the plot
    netName  <- paste("haploNet", analysisName, plotName, sep="-")
    plotFilenameRoot  <- paste(plotsFolder, netName, sep="/")
    initializeGraphics (getGraphicsFilename (plotFilenameRoot), widthInch=24, heightInch=24)
    
    plot(net, size=attr(net, "freq"), scale.ratio=2, bg=pieCols, cex=0.8, pie=pieColTable, labels=FALSE, fast=TRUE, show.mutation=2)
    graphics::legend("bottomright", legendData$label, col="black", pt.bg=legendData$colour, pch=21)
    
    grDevices::dev.off()
}
