###################################################################
# Root Context - Contains all data needed from analysis
###################################################################
#
context.createRootContext <- function (grcData, config, clearCacheData=FALSE) {
    newCtx <- new.env()
    #
    # Clean out the data cache if specified
    #
    newCtx$id         <- digest::digest(grcData, algo="md5")
    if (clearCacheData) {
        print("Clearing data cache")
        cacheRootFolder <- paste(config$folder.data, newCtx$id, sep="/")
        if (file.exists(cacheRootFolder)) {
            unlink(cacheRootFolder, recursive=TRUE)
        }
    }
    #
    newCtx$rootCtx    <- newCtx			# This is a root context (user-created)
    newCtx$meta       <- grcData		# TODO: We should strip the data columns
    newCtx$config     <- config
    #
    # Parse all the data in the genotype columns (phenotype-relevant and barcoding)
    #
    print("Processing phenotype-relevant genotypes")
    alleleGenoCacheFile <- getContextCacheFile(newCtx, "root", "alleles", "alleleGenotypes")
    if (rdaFileExists(alleleGenoCacheFile)) {
        print("Retrieving allele genotypes from cache")
        alleleGenoData <- readRdaSampleData (alleleGenoCacheFile)
    } else {
        print("Parsing allele genotypes from GRC data file")
        #
        # There is no ready-made cached data, so have to parse the genotypes from scratch
        #
        alleleGenoData  <- genotype.processGenotypes (newCtx$meta, config$alleleColumns)
        writeRdaSampleData(alleleGenoData, alleleGenoCacheFile) 
    }
    newCtx$alleleGenoData <- alleleGenoData;
    #
    newCtx$sampleSets <- new.env()
    newCtx$userGraphicAttributes <- new.env()
    #
    #
    dataset.createUnfilteredDataset (newCtx)
    dataset.createFilteredDataset(newCtx)
    dataset.createImputedDataset(newCtx)
    newCtx
}
#
###################################################################
#
###################################################################
#
context.getRootContext <- function (ctx) {
    ctx$rootCtx
}
#
context.getConfig <- function (ctx) {
    ctx$rootCtx$config
}
#
context.getSampleSet <- function (ctx, sampleSetName) {
    ctx$rootCtx$sampleSets[[sampleSetName]]
}
#
context.getDistanceMatrix <- function (ctx, sampleSetName=NULL, useImputation=TRUE) {
    if (useImputation) {							#; print(useImputation)
        datasetName <- "imputed"
    } else {
        datasetName <- "filtered"
    }
    sampleNames <- context.getSamples (ctx, datasetName, sampleSetName)		#; print(length(sampleNames)); print(head(sampleNames))
    dataset  <- ctx$rootCtx[[datasetName]]
    distData <- dataset$distData    						#; print(distData[1:10,1:10]); print(nrow(distData)); print(rownames(distData))
    #print(sampleNames[which(!(sampleNames %in% colnames(distData)))])
    distData <- distData[sampleNames,sampleNames]
    distData
}
#
context.getSamples <- function (ctx, datasetName="unfiltered", sampleSetName=NULL) {
    dataset <- ctx$rootCtx[[datasetName]]
    sampleNames <- dataset$samples
    if (!is.null(sampleSetName)) {
        sampleSet <- ctx$rootCtx$sampleSets[[sampleSetName]]
        sampleSetNames <- sampleSet$samples
        sampleNames <- sampleNames[which(sampleNames %in% sampleSetNames)]
    }
    sampleNames
}
#
context.getMeta <- function (ctx, datasetName=NULL) {
    sampleNames <- NULL
    if (!is.null(datasetName)) {
        dataset <- ctx[[datasetName]]
        sampleNames <- dataset$samples
    }
    meta <- context.getSampleMeta (ctx, sampleNames) 
    meta
}
#
context.getSampleMeta <- function (ctx, sampleNames=NULL) {
    meta <- ctx$rootCtx$meta
    if (!is.null(sampleNames)) {
        meta <- meta[sampleNames,]
    }
    meta
}
#
###################################################################
# Create a new context, consisting of the same data as the source 
# context, but filtered to a given sample list
###################################################################
#
context.trimContext <- function (ctx, sampleNames) {
    trimCtx <- new.env()
    trimCtx$rootCtx <- ctx$rootCtx
    trimDataset <- dataset.trimDatasetBySample (ctx$unfiltered, sampleNames)
    dataset.addDatasetToContext (trimCtx, trimDataset)
    if (!is.null(ctx$filtered)) {
        trimDataset <- dataset.trimDatasetBySample (ctx$filtered, sampleNames)
        dataset.addDatasetToContext (trimCtx, trimDataset)
    }
    if (!is.null(ctx$imputed)) {
        trimDataset <- dataset.trimDatasetBySample (ctx$imputed, sampleNames)
        dataset.addDatasetToContext (trimCtx, trimDataset)
    }
    trimCtx$config     <- ctx$config
    trimCtx$sampleSets <- ctx$sampleSets
    trimCtx
}
#
#
#
context.trimContextByTimeInterval <- function (ctx, datasetName, interval) {	#; print(datasetName); print(interval)
    # Get the sample metadata and filter it by time interval
    sampleMeta <- context.getMeta (ctx, datasetName)				#; print(nrow(sampleMeta))
    filterMeta <- meta.filterByDate (sampleMeta, interval$start, interval$end)	#; print(nrow(filterMeta))
    if (is.null(filterMeta)) {
        return(NULL)
    }
    filterSamples <- rownames(filterMeta)					#; print(filterSamples[1:10])
    trimCtx <- context.trimContext (ctx, filterSamples)
    trimCtx
}
