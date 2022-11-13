###################################################################
# Root Context - Contains all data needed from analysis
###################################################################
#
analysis.createContext <- function (grcData, config) {
    newCtx <- new.env()
    newCtx$id         <- digest::digest(grcData, algo="md5")
    newCtx$config     <- config
    newCtx$sampleSets <- new.env()
    newCtx$userGraphicAttributes <- new.env()
    
    #newCtx <- list(
    #    id         = digest::digest(grcData, algo="md5"),
    #    config     = config,
    #    sampleSets = list()
    #)
    #
    # Include basic (Unfiltered) Data - complete but raw data - built upon initialization
    # print("Initializing Unfiltered Dataset")
    #newCtx$unfiltered <- list(name="unfiltered")
    
    print("Initializing Unfiltered Dataset")
    unfilteredDs <- analysis.createContextDataset (newCtx, "unfiltered")

    meta.setDatasetMeta(newCtx, "unfiltered", grcData, store=TRUE)
    print(paste("Loaded metadata - Samples:", nrow(newCtx$unfiltered$meta)))
    
    #newCtx <- analysis.createFilteredDataset(newCtx)
    #newCtx <- impute.createImputedDataset(newCtx)
    
    analysis.createFilteredDataset(newCtx)
    impute.createImputedDataset(newCtx)
    newCtx
}
#
analysis.createContextDataset <- function (ctx, datasetName) {
    ds <- new.env()
    ds$name <- datasetName
    ctx[[datasetName]] <- ds
    ds
}
#
###################################################################
# Filtered Data, with high-missingness barcodes removed
###################################################################
#
analysis.createFilteredDataset <- function (ctx, loadFromCache=TRUE) {

    print("Initializing Filtered Dataset")
    filteredDs <- analysis.createContextDataset (ctx, "filtered")
    #ctx$filtered <- list(name="filtered")
    #filteredDs   <- ctx$filtered
    
    unfilteredDs <- ctx$unfiltered
    config       <- ctx$config

    filteredMetaFile        <- meta.getMetaDataFile(ctx, "filtered")
    filteredBarcodeFile     <- barcode.getBarcodeDataFile(ctx, "filtered")

    if (loadFromCache & file.exists(filteredMetaFile) & file.exists(filteredBarcodeFile)) {
        meta <- readSampleData (filteredMetaFile)		#; print(colnames(meta))
        meta.setDatasetMeta (ctx, "filtered", meta, store=FALSE)
        barcodeData <- readSampleData (filteredBarcodeFile)
        barcode.setDatasetBarcodes (ctx, "filtered", barcodeData, store=FALSE)
        print(paste("Loaded filtered barcodes - Samples:", nrow(barcodeData), "x SNPs:", ncol(barcodeData)))
    } else {
        meta.setDatasetMeta (ctx, "filtered", ctx$unfiltered$meta, store=FALSE)
        barcode.initializeBarcodes (ctx, "filtered")
        # Trim the metadata to cover the barcodes selected
        sampleNames   <- rownames(ctx$filtered$barcodes)
        filteredMeta  <- filterSampleData (ctx$unfiltered$meta, sampleNames)
        meta.setDatasetMeta (ctx, "filtered", filteredMeta)
    }

    # Get the genotypes, distance matrix and execute the pop structure analysis
    geno.initialize(ctx, "filtered")
    distance.initialize(ctx, "filtered")
}
#
###################################################################
# Create a new context, consisting of the same data as the source context, but filtered to a given sample list
###################################################################
#
analysis.trimContext <- function (ctx, sampleNames) {
    #trimCtx <- list()
    trimCtx <- new.env()
    analysis.addTrimmedDatasetToContext (ctx, "unfiltered", sampleNames, trimCtx)
    if (!is.null(ctx$filtered)) {
        analysis.addTrimmedDatasetToContext (ctx, "filtered", sampleNames, trimCtx)
    }
    if (!is.null(ctx$imputed)) {
        analysis.addTrimmedDatasetToContext (ctx, "imputed", sampleNames, trimCtx)
    }
    trimCtx$config     <- ctx$config
    trimCtx$sampleSets <- ctx$sampleSets
    trimCtx
}
#
#
#
analysis.trimContextByTimeInterval <- function (ctx, interval) {		#; print(interval)
    # Get the sample metadata and filter it by time interval
    dataset <- ctx$unfiltered
    sampleMeta <- dataset$meta							#; print(nrow(sampleMeta))
    filterMeta <- meta.filterByDate (sampleMeta, interval$start, interval$end)	#; print(nrow(filterMeta))
    if (is.null(filterMeta)) {
        return(NULL)
    }
    filterSamples <- rownames(filterMeta)					#; print(filterSamples[1:10])
    trimCtx <- analysis.trimContext (ctx, filterSamples)
    trimCtx
}
#
#
#
analysis.addTrimmedDatasetToContext <- function (ctx, datasetName, sampleNames, trimCtx) {
    # Use only samples that are shared with the original context
    dataset <- ctx[[datasetName]]
    dsSampleNames <- rownames(dataset$meta)
    sampleNames <- sampleNames[which(sampleNames %in% dsSampleNames)]
    #
    #trimCtx[[datasetName]] <- list(name=dataset$name)
    analysis.createContextDataset (trimCtx, dataset$name)

    meta.setDatasetMeta (trimCtx, datasetName, dataset$meta[sampleNames,], store=FALSE)
    if (!is.null(dataset$barcodes)) {
        barcode.setDatasetBarcodes (trimCtx, datasetName, dataset$barcodes[sampleNames,], store=FALSE)
    }
    if (!is.null(dataset$genos)) {
        geno.setDatasetGenotypes (trimCtx, datasetName, dataset$genos[sampleNames,], store=FALSE)
    }
    if (!is.null(dataset$distance)) {
        distance.setDatasetDistance (trimCtx, datasetName, dataset$distance[sampleNames,sampleNames], store=FALSE)
    }
    trimCtx
}
#
###################################################################
# 
###################################################################
#
analysis.selectSampleSet <- function (userCtx, sampleSetName, select) {
    sampleMeta <- userCtx$unfiltered$meta

    # Select the samples to be analyzed
    sampleMeta <- meta.select(sampleMeta, select)

    # Create a trimmed analysis context containing only data pertaining to the selected samples
    sampleNames <- rownames(sampleMeta)
    trimCtx <- analysis.trimContext (userCtx, sampleNames)
    trimCtx$sampleSets <- NULL
    
    sampleSet <- new.env()
    sampleSet$name=sampleSetName
    sampleSet$select=select
    sampleSet$samples=sampleNames
    sampleSet$ctx=trimCtx
    sampleSet$clusters=new.env()
    sampleSet$valuePalettes=new.env()
    
    userCtx$sampleSets[[sampleSetName]] <- sampleSet
    
    metaOutFolder <- getOutFolder(userCtx$config, c(sampleSetName, "metadata"))
    metaFilename  <- paste(metaOutFolder, "/meta-", sampleSetName, "-unfiltered.tab", sep="")
    utils::write.table(trimCtx$unfiltered$meta, file=metaFilename, sep="\t", quote=FALSE, row.names=FALSE)
    metaFilename  <- paste(metaOutFolder, "/meta-", sampleSetName, "-filtered.tab", sep="")
    utils::write.table(trimCtx$filtered$meta, file=metaFilename, sep="\t", quote=FALSE, row.names=FALSE)

    unfilteredCount <- length(sampleNames)
    filteredCount <- nrow(trimCtx$filtered$meta)
    print(paste0("Selected ", unfilteredCount, " samples for dataset '", sampleSetName, "', including ", filteredCount, " quality filtered samples"))
    userCtx
}
#
###################################################################
# Execute
###################################################################
#
analysis.executeOnSampleSet <- function(userCtx, sampleSetName, tasks, params) {
    if (!sampleSetName %in% names(userCtx$sampleSets)) {
        stop(paste("Sample set not initialized:", sampleSetName))
    }
    print(paste("Analyzing sample set:", sampleSetName))
    measures    <- analysis.getParamIfDefined ("analysis.measures", params)
    aggregation <- analysis.getParamIfDefined ("aggegation.levels", params)

    # Execute computations and plots with different tasks
    for (tIdx in 1:length(tasks)) {
        # Split the task if there is a method part in it
        tParts <- unlist(strsplit(tasks[tIdx], "/"))
        task <- tParts[1]
        method <- NULL
        if (length(tParts)>1) {
            method <- tParts[2]
        }
        print(paste(sampleSetName, "-", ifelse(is.null(method), task, paste(task, method, sep="/"))))

        # Execute the task
        if (task == "pca") {
            pca.execute (userCtx, sampleSetName, method, params)

        } else if (task == "tree") {
            tree.execute (userCtx, sampleSetName, method, params)

        } else if (task == "graph") {
            clusterGraph.execute (userCtx, sampleSetName, params)

        } else if (task == "map") {
            interval <- NULL
            if (method %in% c("drug", "mutation", "alleleProp", "diversity", "sampleCount", "location")) {
                intervals <- params$analysis.timeIntervals
                for (idx in 1:length(intervals)) {
                    interval <- intervals[[idx]]
                    map.execute(userCtx, sampleSetName, interval, method, aggregation, measures, params)
                }
            } else {
                map.execute(userCtx, sampleSetName, interval, method, aggregation, measures, params)
            }

        } else {
            stop(paste("Invalid analysis task:", task))
        }
    }
    print(paste("Analysis", sampleSetName, "completed"))
}
