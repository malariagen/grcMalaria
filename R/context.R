###################################################################
# Root Context - Contains all data needed from analysis
###################################################################
#
context.createRootContext <- function (grcData, config, clearCacheData=FALSE) {
    newCtx <- new.env()

    newCtx$id         <- digest::digest(grcData, algo="md5")
    if (clearCacheData) {
        print("Clearing data cache")
        cacheRootFolder <- paste(config$folder.data, newCtx$id, sep="/")
        if (file.exists(cacheRootFolder)) {
            unlink(cacheRootFolder, recursive=TRUE)
        }
    }
    newCtx$rootCtx    <- newCtx			# This is a root context (user-created)
    newCtx$config     <- config
    newCtx$sampleSets <- new.env()
    newCtx$userGraphicAttributes <- new.env()
    
    context.createUnfilteredDataset (newCtx, grcData)
    filter.createFilteredDataset(newCtx)
    impute.createImputedDataset(newCtx)
    newCtx
}
#
context.createUnfilteredDataset <- function (ctx, grcData) {
    print("Initializing Unfiltered Dataset")
    unfilteredDs <- context.createDataset (ctx, "unfiltered")
    meta.setDatasetMeta(ctx, "unfiltered", grcData, store=TRUE)
    print(paste("Loaded metadata - Samples:", nrow(ctx$unfiltered$meta)))
}
#
context.createDataset <- function (ctx, datasetName) {
    ds <- new.env()
    ds$name <- datasetName
    ctx[[datasetName]] <- ds
    ds
}
#
###################################################################
# Create a new context, consisting of the same data as the source context, but filtered to a given sample list
###################################################################
#
context.trimContext <- function (ctx, sampleNames) {
    #trimCtx <- list()
    trimCtx <- new.env()
    trimCtx$rootCtx <- ctx$rootCtx 

    context.addTrimmedDatasetToContext (ctx, "unfiltered", sampleNames, trimCtx)
    if (!is.null(ctx$filtered)) {
        context.addTrimmedDatasetToContext (ctx, "filtered", sampleNames, trimCtx)
    }
    if (!is.null(ctx$imputed)) {
        context.addTrimmedDatasetToContext (ctx, "imputed", sampleNames, trimCtx)
    }
    trimCtx$config     <- ctx$config
    trimCtx$sampleSets <- ctx$sampleSets
    trimCtx
}
#
#
#
context.trimContextByTimeInterval <- function (ctx, interval) {		#; print(interval)
    # Get the sample metadata and filter it by time interval
    dataset <- ctx$unfiltered
    sampleMeta <- dataset$meta							#; print(nrow(sampleMeta))
    filterMeta <- meta.filterByDate (sampleMeta, interval$start, interval$end)	#; print(nrow(filterMeta))
    if (is.null(filterMeta)) {
        return(NULL)
    }
    filterSamples <- rownames(filterMeta)					#; print(filterSamples[1:10])
    trimCtx <- context.trimContext (ctx, filterSamples)
    trimCtx
}
#
#
#
context.addTrimmedDatasetToContext <- function (ctx, datasetName, sampleNames, trimCtx) {
    # Use only samples that are shared with the original context
    dataset <- ctx[[datasetName]]
    dsSampleNames <- rownames(dataset$meta)
    sampleNames <- sampleNames[which(sampleNames %in% dsSampleNames)]
    #
    #trimCtx[[datasetName]] <- list(name=dataset$name)
    context.createDataset (trimCtx, dataset$name)

    meta.setDatasetMeta (trimCtx, datasetName, dataset$meta[sampleNames,], store=FALSE)
    if (!is.null(dataset$barcodes)) {
        barcode.setDatasetBarcodes (trimCtx, datasetName, dataset$barcodes[sampleNames,], store=FALSE)
    }
    if (!is.null(dataset$genos)) {
        geno.setDatasetGenotypes (trimCtx, datasetName, dataset$genos[sampleNames,], store=FALSE)
    }
    trimCtx
}
