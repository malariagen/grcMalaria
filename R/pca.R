phylo.pcCount <- 10
###############################################################################
# PCA Analysis.
# Computes and plots PCA, using a pairwise distance matrix as an input.
# Multiple analyses, based on different subsets of samples, can be performed in a single run.
#
# Supported methods:
#     pca/PCoA		Analyzes a distance matrix
#     pca/nipals	Analyzes each barcode variant as a variable
#     pca/bpca		Analyzes each barcode variant as a variable
#     pca/ppca		Analyzes each barcode variant as a variable
################################################################################
#
pca.execute <- function(userCtx, sampleSetName, pcaMethod) {
    sampleSet <- userCtx$sampleSets[[sampleSetName]]
    ctx <- sampleSet$ctx
    dataset <- ctx$imputed

    dataFolder <- getOutFolder(ctx, sampleSetName, c(pcaMethod, "data"))
    sampleMeta <- dataset$meta
    distData  <- dataset$distance
    genosData <- dataset$genos
    
    # Compute the Principal components, and join them to the metadata
    pcNames <- paste("PC", seq(1:phylo.pcCount), sep="")
    sampleNames <- rownames(sampleMeta)
    
    pcScores <- NULL
    varExplained <- NULL
        
    if (pcaMethod == "umap") {
        #genos <- as.matrix(genosData)
        #pcaResults <- umap(genos)
        #pcScores <- pcaResults$layout
        
        
        
        #pcScores <- pcaResults@scores
        #varExplained <- pcaResults@R2
    } else if (pcaMethod == "PCoA") {
        pcaResults <- stats::cmdscale(as.matrix(distData), eig=TRUE, k=phylo.pcCount)
        pcScores <- pcaResults$points
        varExplained <- pcaResults$eig / sum(abs(pcaResults$eig))
    } else  {
        genos <- as.matrix(genosData)
        pcaResults <- NULL
        if (pcaMethod == "bpca") {
            pcaResults <- pcaMethods::pca(genos, method=pcaMethod, center=TRUE, scale="uv", nPcs=phylo.pcCount, maxSteps=500)
        } else {
            pcaResults <- pcaMethods::pca(genos, method=pcaMethod, center=TRUE, scale="uv", nPcs=phylo.pcCount)
        }
        pcScores <- pcaResults@scores
        varExplained <- pcaResults@R2
    }
        
    rownames(pcScores) <- sampleNames
    colnames(pcScores) <- pcNames
    pcaScoresFilename  <- paste(dataFolder, pca.getDataFileName("", sampleSetName, pcaMethod), sep="/")
    writeSampleData(pcScores, pcaScoresFilename)

    pcaVarData <- data.frame(pcNames, varExplained[1:length(pcNames)])
    colnames(pcaVarData) <- c("PC","VarianceExplained")
    pcaVarFilename  <- paste(dataFolder, pca.getDataFileName("-varExplained", sampleSetName, pcaMethod), sep="/")
    utils::write.table(pcaVarData, file=pcaVarFilename, sep="\t", quote=FALSE, row.names=FALSE)
}

pca.getDataFileName <- function(suffix, sampleSetName, pcaMethod) {
    fn <- paste("pca", suffix, "-", sampleSetName, "-", pcaMethod, ".tab", sep="")
    fn
}

pca.executePlots <- function(userCtx, sampleSetName, pcaMethod, plotList) {
    sampleSet <- userCtx$sampleSets[[sampleSetName]]
    ctx <- sampleSet$ctx
    dataset <- ctx$imputed

    sampleMeta <- dataset$meta

    dataFolder      <- getOutFolder(ctx, sampleSetName, c(pcaMethod, "data"))
    plotsRootFolder <- getOutFolder(ctx, sampleSetName, c(pcaMethod, "plots"))

    # Read in the PCA results and attach them to the metadata
    pcaScoresFilename  <- paste(dataFolder, pca.getDataFileName("", sampleSetName, pcaMethod), sep="/")
    pcaData <- readSampleData(pcaScoresFilename)
    sampleMeta <- cbind(sampleMeta, pcaData)
    
    # Execute the plots
    #print(plotList)
    for (plotIdx in 1:length(plotList)) {
        plotDef <- plotList[[plotIdx]]
        plotDefName <- plotDef$name
        plotName <- paste(sampleSetName, plotDefName, pcaMethod, sep="-")
        print (paste("PCA Plot: ",plotName))
  
        # Set up the graphical attributes for rendering
        #print(colnames(sampleMeta))
        #print(plotDef$render)
        dataList <- applyGraphicalAttributes(sampleMeta, plotDef$render)
        plotMetadata <- dataList$meta
        legendData <- dataList$legend
        
        # Order the samples so they are plotted in the correct stack order
        #print(plotMetadata$plot__order)
        excludedIdx <- which(as.numeric(plotMetadata$plot__order) < 0)
        if (length(excludedIdx) > 0) {
            print(paste("Excluded samples:", length(excludedIdx)))
            plotMetadata <- plotMetadata[-excludedIdx,]
        }
        plotMetadata <- plotMetadata[order(-plotMetadata$plot__order),]

        # Do the plots
        plotsFolder  <- getSubFolder (plotsRootFolder, plotDefName)
        pca.plotThreeComponents (plotName, plotMetadata, legendData, "PC1", "PC2", "PC3", plotsFolder);
        pca.plotThreeComponents (plotName, plotMetadata, legendData, "PC1", "PC4", "PC5", plotsFolder);
        pca.plotThreeComponents (plotName, plotMetadata, legendData, "PC1", "PC6", "PC7", plotsFolder);
        pca.plotThreeComponents (plotName, plotMetadata, legendData, "PC1", "PC8", "PC9", plotsFolder);
    }
}

###############################################################################
# Principal Component plotting
###############################################################################
pca.plotPrincipalComponents <- function(plotName, sampleMeta, legendData, 
                                    pcAIdx, pcBIdx, pcCIdx=NULL, plotsFolder) {
    plotThree <- !is.null(pcCIdx)
    sampleCount <- nrow(sampleMeta)
  
    #Get the graphical parameters for each sample
    sampleColours <- sampleMeta$plot__colour
    samplePch     <- sampleMeta$plot__pch
    sampleSize    <- sampleMeta$plot__size
    sampleLwd     <- sampleMeta$plot__lwd
    sampleLcolours<- sampleMeta$plot__lcolour
  
    # Plot
    plotFilename  <- paste("pca-",plotName,"-",pcAIdx,"_",pcBIdx, sep="")
    if (plotThree) {
        plotFilename  <- paste(plotFilename, pcCIdx, sep="_")
    }
    print(paste("Plotting ", plotFilename))
    graphicFilenameRoot  <- paste(plotsFolder, plotFilename, sep="/")
    initializeGraphics (getGraphicsFilename (graphicFilenameRoot), widthInch=12, heightInch=8, resolution=300)
    #print(colnames(sampleMeta))
  
    if (plotThree) {
        graphics::par(mfrow=c(2,2))
    } else {
#        graphics::par(mfrow=c(1,2))
    }
    graphics::par(mar=c(4,4,1,1) + 0.1)
    x <- sampleMeta[,pcAIdx]
    y <- sampleMeta[,pcBIdx]
    plot(x,y, xlab=pcAIdx, ylab=pcBIdx, cex.lab=0.8, cex.axis=0.8, bg=sampleColours, col=sampleLcolours, pch=samplePch, lwd=sampleLwd, cex=sampleSize)
    if (plotThree) {
        x <- sampleMeta[,pcCIdx]
        y <- sampleMeta[,pcBIdx]
        plot(x,y, xlab=pcCIdx, ylab=pcBIdx, cex.lab=0.8, cex.axis=0.8, bg=sampleColours, col=sampleLcolours, pch=samplePch, lwd=sampleLwd, cex=sampleSize)
        x <- sampleMeta[,pcAIdx]
        y <- sampleMeta[,pcCIdx]
        plot(x,y, xlab=pcAIdx, ylab=pcCIdx, cex.lab=0.8, cex.axis=0.8, bg=sampleColours, col=sampleLcolours, pch=samplePch, lwd=sampleLwd, cex=sampleSize)
        plot(x,y,xaxt="n",xlab="", ylab="",yaxt="n",type="n",bty="n")
    }
    grDevices::dev.off()
    
    # Do a separate legend- useful if the names are long
    plotFilename  <- paste("pca-",plotName,"-legend", sep="")
    graphicFilenameRoot  <- paste(plotsFolder, plotFilename, sep="/")
    initializeGraphics (getGraphicsFilename (graphicFilenameRoot), widthInch=8, heightInch=12, resolution=300)
    #graphics::par(mar=c(1,1,1,1) + 0.1)
    #plot(x,y,xaxt="n",xlab="", ylab="",yaxt="n",type="n",bty="n")
    graphics::plot.new()
    graphics::par(mar=c(0,0,0,0))
    graphics::legend("topleft", ncol=1, inset=0.05, cex=1.0, 
             legendData$label, pt.bg=legendData$colour, col=legendData$lcolour, 
             pch=as.numeric(legendData$pch), pt.lwd=as.numeric(legendData$lwd), 
             pt.cex=as.numeric(legendData$size))
    grDevices::dev.off()
}

pca.plotThreeComponents <- function(plotName, sampleMeta, legendData, pcAIdx, pcBIdx, pcCIdx, plotsFolder) {
    pca.plotPrincipalComponents (plotName, sampleMeta, legendData, pcAIdx, pcBIdx, pcCIdx, plotsFolder)
}

pca.plotTwoComponents <- function(plotName, sampleMeta, legendData, pcAIdx, pcBIdx, plotsFolder) {
    pca.plotPrincipalComponents (plotName, sampleMeta, legendData, pcAIdx, pcBIdx, NULL, plotsFolder)
}

