##############################################################################
# Basic dataset functions
################################################################################
#
# Creates an empty (data-free) dataset.
#
dataset.createDataset <- function (datasetName) {
    ds <- new.env()
    ds$name <- datasetName
    ds
}
#
# Creates an empty (data-free) dataset and attaches it to the context, under the name specified.
#
dataset.addDatasetToContext <- function (ctx, dataset) {
    datasetName <- dataset$name
    ctx[[datasetName]] <- dataset
}
#
# Returns a dataset that retains only the selected samples of the source dataset.
#
dataset.trimDatasetBySample <- function (srcDataset, selSamples) {
    # Some of the selected samples may be filtered out of the doriginal dataset
    trimSamples <- selSamples[which(selSamples %in% srcDataset$samples)]
    dataset <- dataset.createDataset (srcDataset$name)
    if (!is.null(srcDataset$barcodeGenoData)) {
        dataset$barcodeGenoData <- genotype.filterGenotypeDataBySample (srcDataset$barcodeGenoData, trimSamples)
    }
    dataset$samples <- trimSamples;
    dataset
}
#
##############################################################################
# Dataset creation functions
################################################################################
#
#
#
dataset.createUnfilteredDataset <- function (ctx) {
    print("Initializing Unfiltered Dataset")
    config <- context.getConfig(ctx)
    datasetName <- "unfiltered"
    dataset <- dataset.createDataset (datasetName)
    dataset.addDatasetToContext (ctx, dataset)
    #
    meta <- context.getMeta(ctx)
    dataset$samples <- rownames(meta)
    print(paste("Loaded sample data - Samples:", length(dataset$samples)))
    #
    #
    # If there is cached ready-parsed barcode genotype data, use it, otherwise parse genotypes from scratch
    #
    print("Processing barcoding genotypes")
    barcodeGenoCacheFile <- getContextCacheFile(ctx$rootCtx, datasetName, "barcodes", "barcodeGenotypes")
    if (rdaFileExists(barcodeGenoCacheFile)) {
        print("Retrieving barcode genotypes from cache")
        barcodeGenoData <- readRdaSampleData (barcodeGenoCacheFile)
    } else {
        print("Parsing barcode genotypes from GRC data file")
        barcodeGenoData  <- genotype.processGenotypes (meta, config$barcodeColumns)
        writeRdaSampleData(barcodeGenoData, barcodeGenoCacheFile) 
    }
    dataset$barcodeGenoData <- barcodeGenoData;
}
#
#
# Filtered Data, with high-missingness barcodes removed
#
dataset.createFilteredDataset <- function (ctx) {
    print("Initializing Filtered Dataset")
    config <- context.getConfig(ctx)
    datasetName <- "filtered"
    dataset <- dataset.createDataset (datasetName)
    dataset.addDatasetToContext (ctx, dataset)
    #
    # If there is cached ready-trimmed barcode genotype data, use it, otherwise filter genotypes
    #
    print("Processing samples for missingness")
    barcodeGenoCacheFile <- getContextCacheFile(ctx$rootCtx, datasetName, "barcodes", "barcodeGenotypes")
    if (rdaFileExists(barcodeGenoCacheFile)) {
        print("Retrieving filtered barcode genotypes from cache")
        barcodeGenoData <- readRdaMatrix (barcodeGenoCacheFile)
        
    } else {
        print("Filtering barcode genotypes")
        #
        # Remove "hopeless" samples with missingness > 0.5
        #
        print("Removing samples missing more than half their genotypes")
        barcodeGenoData <- ctx$rootCtx$unfiltered$barcodeGenoData
        barcodeGenoData <- dataset.filterBarcodesBySampleMissingProp (barcodeGenoData, 0.5)
        #
        # Remove SNPs with too much missingness
        #
        print("Filtering barcode genotypes by column missingness")
        maxMissingProp <- 1 - config$minSnpTypability
        barcodeGenoData <- dataset.filterBarcodesByColumnMissingProp (barcodeGenoData, maxMissingProp)
        #
        # Remove samples with too much missingness
        #
        print("Filtering barcode genotypes by sample missingness")
        maxMissingProp <- 1 - config$minSampleTypability
        barcodeGenoData <- dataset.filterBarcodesBySampleMissingProp (barcodeGenoData, maxMissingProp)
        #
        print(paste("After filtering, barcoding uses", length(barcodeGenoData$columns), "variants for", 
                     length(barcodeGenoData$samples), "samples."))
        writeRdaSampleData(barcodeGenoData, barcodeGenoCacheFile) 
    }
    dataset$barcodeGenoData <- barcodeGenoData;
    dataset$samples <- barcodeGenoData$samples;
    #
    # Get the sample pairwise distance matrix
    #
    print("Getting pairwise genetic distances")
    distanceMatrixCacheFile <- getContextCacheFile(ctx$rootCtx, datasetName, "distance", "sampleDistance")
    if (rdaFileExists(distanceMatrixCacheFile)) {
        print("Retrieving pairwise genetic distances from cache")
        distData <- readRdaMatrix (distanceMatrixCacheFile)
    } else {
        print("Estimating pairwise genetic distances - please wait")
        distData <- NULL
        distData <- dist_calculateDistanceMatrix (barcodeGenoData)
        writeRdaMatrix(distData, distanceMatrixCacheFile) 
    }
    dataset$distData <- distData;
    ctx
}
#
# Remove samples with too much missingness
#
dataset.filterBarcodesBySampleMissingProp <- function (barcodeGenoData, maxMissingProp) {
        colCount      <- barcodeGenoData$columnCount
        missingCounts <- barcodeGenoData$sampleMissingCounts
        missingProps  <- missingCounts / colCount
        samples        <- barcodeGenoData$samples
        selSamples     <- samples[which(missingProps <= maxMissingProp)]
        barcodeGenoData <- genotype.filterGenotypeDataBySample (barcodeGenoData, selSamples)
        barcodeGenoData
}
#
# Remove SNPs with too much missingness
#
dataset.filterBarcodesByColumnMissingProp <- function (barcodeGenoData, maxMissingProp) {
        sampleCount   <- barcodeGenoData$sampleCount
        missingCounts <- barcodeGenoData$columnMissingCounts
        missingProps  <- missingCounts / sampleCount
        columnNames    <- names(missingCounts)
        selColumns  <- columnNames[which(missingProps <= maxMissingProp)]
        barcodeGenoData <- genotype.filterGenotypeDataByColumn (barcodeGenoData, selColumns)
        barcodeGenoData
}
#
#
#
dataset.createImputedDataset <- function (ctx) {
    print("Initializing Imputed Dataset")
    config <- context.getConfig(ctx)
    datasetName <- "imputed"
    dataset <- dataset.createDataset (datasetName)
    dataset.addDatasetToContext (ctx, dataset)
    #
    # If there is cached ready-imputed barcode genotype data, use it, otherwise filter genotypes
    #
    barcodeGenoCacheFile <- getContextCacheFile(ctx$rootCtx, datasetName, "barcodes", "barcodeGenotypes")
    if (rdaFileExists(barcodeGenoCacheFile)) {
        print("Retrieving imputed barcode genotypes from cache")
        barcodeGenoData <- readRdaSampleData (barcodeGenoCacheFile)
    } else {
        print("Imputing barcode genotypes")
        #
        # We have to first determine which samples require lower levels of imputation 
        #
        barcodeGenoData <- ctx$rootCtx$filtered$barcodeGenoData
        #
        print("Filtering barcode genotypes by sample missingness")
        selSamples      <- impute.selectImputableSamples (ctx, barcodeGenoData)
        barcodeGenoData <- genotype.filterGenotypeDataBySample (barcodeGenoData, selSamples)
        barcodeGenoData <- impute.imputeBarcodeGenos (ctx, barcodeGenoData)
        print(paste("After imputation, barcoding uses", length(barcodeGenoData$samples), "samples."))
        writeRdaSampleData(barcodeGenoData, barcodeGenoCacheFile) 
    }
    dataset$barcodeGenoData <- barcodeGenoData;
    dataset$samples <- barcodeGenoData$samples;
    #
    # Get barcode table and sample barcode strings
    #
    dataset$barcodeGenoTable <- impute.buildBarcodeGenotypeTable(barcodeGenoData);
    dataset$barcodes <- impute.buildBarcodeStrings(dataset$barcodeGenoTable);
    #
    # Make a distance matrix from imputed data
    #
    print("Getting pairwise genetic distances")
    distanceMatrixCacheFile <- getContextCacheFile(ctx$rootCtx, datasetName, "distance", "sampleDistance")
    if (rdaFileExists(distanceMatrixCacheFile)) {
        print("Retrieving pairwise genetic distances from cache")
        distData <- readRdaMatrix (distanceMatrixCacheFile)
    } else {
        print("Estimating pairwise genetic distances - please wait")
        #distData <- NULL
        distData <- dist_calculateDistanceMatrix (barcodeGenoData)
        writeRdaMatrix(distData, distanceMatrixCacheFile)
    }
    dataset$distData <- distData;
    ctx
}
#
#
#
dataset.getSampleBarcodeGenos <- function (ctx, sampleName, useImputed=FALSE) {
    dataset <- ctx$filtered
    if (useImputed) dataset <- ctx$imputed
    colNames <- as.character(names(dataset$barcodeGenoData$columnGenoData))
    colCount <- length(colNames)
    genos <- vector(mode="integer",length=colCount)
    for (cIdx in 1:colCount) {
        cName <- colNames[cIdx]
        cGenos <- dataset$barcodeGenoData$columnGenoData[[cName]]$sampleGenotypes
        genos[cIdx] <- cGenos[sampleName]
    }
    names(genos) <- colNames
    genos
}
