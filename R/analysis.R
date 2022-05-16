###################################################################
# Root Context - Contains all data needed from analysis
###################################################################
#
analysis.createContext <- function (grcData, config) {
    ctx <- list(
        id         = digest::digest(grcData, algo="md5"),
        config     = config,
        sampleSets = list()
    )
    #
    # Include basic (Unfiltered) Data - complete but raw data - built upon initialization
    #
    print("Initializing Unfiltered Dataset")
    ctx$unfiltered <- list(name="unfiltered")
    ctx <- meta.setDatasetMeta(ctx, "unfiltered", grcData, store=TRUE)
    print(paste("Loaded metadata - Samples:", nrow(ctx$unfiltered$meta)))
    
    ctx <- analysis.createFilteredDataset(ctx)
    ctx <- impute.createImputedDataset(ctx)
    ctx
}
#
###################################################################
# Filtered Data, with high-missingness barcodes removed
###################################################################
#
analysis.createFilteredDataset <- function (ctx, loadFromCache=TRUE) {
    print("Initializing Filtered Dataset")
    ctx$filtered <- list(name="filtered")
    
    filteredDs   <- ctx$filtered
    unfilteredDs <- ctx$unfiltered
    config       <- ctx$config

    filteredMetaFile        <- meta.getMetaDataFile(ctx, "filtered")
    filteredBarcodeFile     <- barcode.getBarcodeDataFile(ctx, "filtered")

    if (loadFromCache & file.exists(filteredMetaFile) & file.exists(filteredBarcodeFile)) {
        meta <- readSampleData (filteredMetaFile)		#; print(colnames(meta))
        ctx <- meta.setDatasetMeta (ctx, "filtered", meta, store=FALSE)
        barcodeData <- readSampleData (filteredBarcodeFile)
        ctx <- barcode.setDatasetBarcodes (ctx, "filtered", barcodeData, store=FALSE)
        print(paste("Loaded filtered barcodes - Samples:", nrow(barcodeData), "x SNPs:", ncol(barcodeData)))
    } else {
        ctx <- meta.setDatasetMeta (ctx, "filtered", ctx$unfiltered$meta, store=FALSE)
        ctx <- barcode.initializeBarcodes (ctx, "filtered")
        # Trim the metadata to cover the barcodes selected
        sampleNames   <- rownames(ctx$filtered$barcodes)
        filteredMeta  <- filterSampleData (ctx$unfiltered$meta, sampleNames)
        ctx <- meta.setDatasetMeta (ctx, "filtered", filteredMeta)
    }

    # Get the genotypes, distance matrix and execute the pop structure analysis
    ctx <- geno.initialize(ctx, "filtered")
    ctx <- distance.initialize(ctx, "filtered")
    ctx
}
#
###################################################################
# Create a new context, consisting of the same data as the source context, but filtered to a given sample list
###################################################################
#
analysis.trimContext <- function (ctx, sampleNames) {
    trimCtx <- list()
    trimCtx <- analysis.addTrimmedDatasetToContext (ctx, "unfiltered", sampleNames, trimCtx)
    if (!is.null(ctx$filtered)) {
        trimCtx <- analysis.addTrimmedDatasetToContext (ctx, "filtered", sampleNames, trimCtx)
    }
    if (!is.null(ctx$imputed)) {
        trimCtx <- analysis.addTrimmedDatasetToContext (ctx, "imputed", sampleNames, trimCtx)
    }
    trimCtx$config     <- ctx$config
    trimCtx
}
#
analysis.trimContextByTimeInterval <- function (ctx, interval) {		#; print(interval)
    # Get the sample metadata and filter it by time interval
    dataset <- ctx$unfiltered
    sampleMeta <- dataset$meta							#; print(nrow(sampleMeta))
    filterMeta <- meta.filterByDate(sampleMeta, interval$start, interval$end)	#; print(nrow(filterMeta))
    if (is.null(filterMeta)) {
        return(NULL)
    }
    filterSamples <- rownames(filterMeta)					#; print(filterSamples[1:10])
    trimCtx <- analysis.trimContext (ctx, filterSamples)
    trimCtx
}
#
analysis.addTrimmedDatasetToContext <- function (ctx, datasetName, sampleNames, trimCtx) {
    # Use only samples that are shared with the original context
    dataset <- ctx[[datasetName]]
    dsSampleNames <- rownames(dataset$meta)
    sampleNames <- sampleNames[which(sampleNames %in% dsSampleNames)]
    #
    trimCtx[[datasetName]] <- list(name=dataset$name)
    trimCtx <- meta.setDatasetMeta (trimCtx, datasetName, dataset$meta[sampleNames,], store=FALSE)
    if (!is.null(dataset$barcodes)) {
        trimCtx <- barcode.setDatasetBarcodes (trimCtx, datasetName, dataset$barcodes[sampleNames,], store=FALSE)
    }
    if (!is.null(dataset$genos)) {
        trimCtx <- geno.setDatasetGenotypes (trimCtx, datasetName, dataset$genos[sampleNames,], store=FALSE)
    }
    if (!is.null(dataset$distance)) {
        trimCtx <- distance.setDatasetDistance (trimCtx, datasetName, dataset$distance[sampleNames,sampleNames], store=FALSE)
    }
    trimCtx
}
#
###################################################################
# 
###################################################################
#
analysis.selectSampleSet <- function (ctx, sampleSetName, select) {
    sampleMeta <- ctx$unfiltered$meta

    # Select the samples to be analyzed
    sampleMeta <- meta.select(sampleMeta, select)

    # Create a trimmed analysis context containing only data pertaining to the selected samples
    sampleNames <- rownames(sampleMeta)
    trimCtx <- analysis.trimContext (ctx, sampleNames)
    
    sampleSet <- list(
        name=sampleSetName,
        select=select,
        samples=sampleNames,
        ctx=trimCtx,
        clusters=c()
    )
    ctx$sampleSets[[sampleSetName]] <- sampleSet
    
    metaOutFolder <- getOutFolder(ctx, c(sampleSetName, "metadata"))
    metaFilename  <- paste(metaOutFolder, "/meta-", sampleSetName, "-unfiltered.tab", sep="")
    utils::write.table(trimCtx$unfiltered$meta, file=metaFilename, sep="\t", quote=FALSE, row.names=FALSE)
    metaFilename  <- paste(metaOutFolder, "/meta-", sampleSetName, "-filtered.tab", sep="")
    utils::write.table(trimCtx$filtered$meta, file=metaFilename, sep="\t", quote=FALSE, row.names=FALSE)

    unfilteredCount <- length(sampleNames)
    filteredCount <- nrow(trimCtx$filtered$meta)
    print(paste0("Selected ", unfilteredCount, " samples for dataset '", sampleSetName, "', including ", filteredCount, " quality filtered samples"))
    ctx
}
#
###################################################################
# Execute
###################################################################
#
analysis.executeOnSampleSet <- function(ctx, sampleSetName, tasks, params) {
    if (!sampleSetName %in% names(ctx$sampleSets)) {
        stop(paste("Sample set not initialized:", sampleSetName))
    }
    print(paste("Analyzing sample set:", sampleSetName))
    measures    <- analysis.getParamIfDefined ("analysis.measures", params)
    aggregation <- analysis.getParamIfDefined ("aggegation.levels", params)

    # Resolve all the automatic rendering in the plots
    plotList <- NULL	# TODO - plots are not yet implemented, will pass them through the params
    #if (!is.null(plotList)) {
    #    # Get the trimmed analysis context containing only data pertaining to the selected samples
    #    sampleSet <- ctx$sampleSets[[sampleSetName]]
    #    plotList <- resolveAutomaticRenderingInPlots (sampleSet$ctx$meta, plotList)
    #}

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
            pca.execute (ctx, sampleSetName, method, params)

        } else if (task == "tree") {
            tree.execute (ctx, sampleSetName, method, params)

        } else if (task == "graph") {
            clusterGraph.execute (ctx, sampleSetName, plotList, params)

        #} else if (task == "haploNet") {
        #    haploNet.execute (ctx, sampleSetName, plotList, params)
        
        } else if (task == "map") {
            interval <- NULL
            if (method %in% c("drug", "mutation", "diversity", "sampleCount")) {
                intervals <- params$analysis.timeIntervals
                for (idx in 1:length(intervals)) {
                    interval <- intervals[[idx]]
                    map.execute(ctx, sampleSetName, interval, method, aggregation, measures, params)
                }
            } else {
                map.execute(ctx, sampleSetName, interval, method, aggregation, measures, params)
            }

        } else {
            stop(paste("Invalid analysis task:", task))
        }
    }
    print(paste("Analysis", sampleSetName, "completed"))
}
