###################################################################
# Imputed Data, with filled-in missingness in filtered barcodes 
###################################################################
#
impute.createImputedDataset <- function (ctx) {
    print("Initializing Imputed Dataset")
    imputedDs <- context.createDataset (ctx, "imputed")
    #
    filteredDs <- ctx$filtered
    config     <- ctx$config
    #
    imputedMetaFile    <- meta.getMetaDataFile(ctx, "imputed")
    imputedBarcodeFile <- barcode.getBarcodeDataFile (ctx, "imputed")
    #
    if (rdaFileExists(imputedMetaFile) & rdaFileExists(imputedBarcodeFile)) {
        meta <- readRdaSampleData (imputedMetaFile)		#; print(colnames(meta))
        meta.setDatasetMeta (ctx, "imputed", meta, store=FALSE)
        barcodeData <- readRdaSampleData (imputedBarcodeFile)
        barcode.setDatasetBarcodes (ctx, "imputed", barcodeData, store=FALSE)
        print(paste("Loaded imputed barcodes - Samples:", nrow(barcodeData), "x SNPs:", ncol(barcodeData)))
    } else {
        meta  <- impute.filterSampleData (ctx)
        meta.setDatasetMeta (ctx, "imputed", meta, store=TRUE)
        sampleNames <- rownames(meta)
        barcodeData <- filteredDs$barcodes
        barcodeData <- barcodeData[sampleNames,]
        impBarcodeData <- impute.imputeBarcodes (ctx, barcodeData)			#; print(head(impBarcodeData))
        barcode.setDatasetBarcodes (ctx, "imputed", impBarcodeData, store=TRUE)
        print(paste("Computed imputed barcodes - Samples:", nrow(barcodeData), "x SNPs:", ncol(barcodeData)))
    }
    #
    # Get the genotypes and distance matrix from imputed data
    #
    geno.initialize(ctx, "imputed")
    distance.initializeDistanceMatrix (ctx, "imputed")
}
#
impute.filterSampleData <- function (ctx) {
    srcMeta  <- ctx$filtered$meta
    srcBarcodes <- ctx$filtered$barcodes
    maxImputeCounts <- ctx$config$maxImputedProportion * ncol(srcBarcodes)    #; print(maxImputeCounts)
    imputeCounts <- apply(srcBarcodes, 1, function(x) length(which((x=="N")|(x=="X"))))   #; print(imputeCounts)
    keepIdx <- which (imputeCounts <= maxImputeCounts)
    imputedMeta <- srcMeta[keepIdx,]
    imputedMeta 
}
#
impute.imputeBarcodes <- function (ctx, barcodeData) {
    inSampleNames <- rownames(barcodeData)
    #
    barcodeMeta <- barcode.getMetadata(ctx, barcodeData)
    refs  <- barcodeMeta$Ref; nrefs <- barcodeMeta$Nonref
    #
    ms <- distance.findMostSimilarSamples (ctx, "filtered")
    msIndexData <- ms$indexes; msDistanceData <- ms$distances
    #
    # Set the sample order to be the same as for the "most similar sample" table,
    # or else we will not be able to use the sample indexes in that table.
    #
    sampleNames <-  colnames(msIndexData)
    barcodeData <- barcodeData[sampleNames,]
    #
    sampleCount <- length(sampleNames)
    snpCount <- ncol(barcodeData)
    impMat <- as.matrix(barcodeData)
    for (sIdx in 1:sampleCount) {
        alleles <- barcodeData[sIdx,]
        impIndexes <- which(alleles %in% c("X","N"))
        if (length(impIndexes) == 0) {
            next
        }
        msIndexes <- msIndexData[,sIdx]
        msScores <- (1 - msDistanceData[,sIdx])
        for (vIdx in impIndexes) {
            msAlleles <- barcodeData[msIndexes,vIdx]
            refScore  <- sum(msScores[which(msAlleles==refs[vIdx])])
            nrefScore <- sum(msScores[which(msAlleles==nrefs[vIdx])])
            impMat[sIdx,vIdx] <- ifelse(refScore >= nrefScore, refs[vIdx], nrefs[vIdx])
        }
    }
    impBarcodeData <- data.frame(impMat)
    rownames(impBarcodeData) <- sampleNames
    colnames(impBarcodeData) <- colnames(barcodeData)		#; print(barcodeData[1:10,1:30]); print(impBarcodeData[1:10,1:30])
    #
    # Restore initial sample order
    #
    impBarcodeData <- impBarcodeData[inSampleNames,]
    impBarcodeData
}
