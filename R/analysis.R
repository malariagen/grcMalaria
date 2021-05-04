###################################################################
# Root Context - Contains all data needed from analysis
###################################################################
#
analysis.createContext <- function (grcData, config) {
    ctx <- list(
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
        ctx <- setDatasetBarcodes (ctx, "filtered", barcodeData, store=FALSE)
        print(paste("Loaded filtered barcodes - Samples:", nrow(barcodeData), "x SNPs:", ncol(barcodeData)))
    } else {
        ctx <- meta.setDatasetMeta (ctx, "filtered", ctx$unfiltered$meta, store=FALSE)
        ctx <- initializeBarcodes (ctx, "filtered")
        # Trim the metadata to cover the barcodes selected
        sampleNames   <- rownames(ctx$filtered$barcodes)
        filteredMeta  <- filterSampleData (ctx$unfiltered$meta, sampleNames)
        ctx <- meta.setDatasetMeta (ctx, "filtered", filteredMeta)
    }

    # Get the genotypes, distance matrix and execute the pop structure analysis
    ctx <- distance.initialize(ctx, "filtered")
    ctx
}
#
###################################################################
# Create a new context, consisting of the same data as the source context, but filtered to a given sample list
###################################################################
#
analysis.selectSamplesInContext <- function (ctx, sampleNames) {
    trimCtx <- list()
    trimCtx <- analysis.addTrimmedDataset (ctx, "unfiltered", sampleNames, trimCtx)
    if (!is.null(ctx$filtered)) {
        trimCtx <- analysis.addTrimmedDataset (ctx, "filtered", sampleNames, trimCtx)
    }
    if (!is.null(ctx$imputed)) {
        trimCtx <- analysis.addTrimmedDataset (ctx, "imputed", sampleNames, trimCtx)
    }
    trimCtx$config     <- ctx$config
    trimCtx
}
#
analysis.addTrimmedDataset <- function (ctx, datasetName, sampleNames, trimCtx) {
    # Use only samples that are shared with the original context
    dataset <- ctx[[datasetName]]
    dsSampleNames <- rownames(dataset$meta)
    sampleNames <- sampleNames[which(sampleNames %in% dsSampleNames)]
    #
    trimCtx[[datasetName]] <- list(name=dataset$name)
    #
    trimMeta <- dataset$meta[sampleNames,]
    trimCtx <- meta.setDatasetMeta (trimCtx, datasetName, trimMeta,store=FALSE)
    #
    if (!is.null(dataset$barcodes)) {
        trimCtx <- barcode.setDatasetBarcodes (trimCtx, datasetName, dataset$barcodes[sampleNames,], store=FALSE)
    }
    if (!is.null(dataset$distance)) {
        trimCtx <- distance.setDatasetDistance (trimCtx, datasetName, ctx$distance[sampleNames,sampleNames], store=FALSE)
    }
    trimCtx
}
#
###################################################################
# 
###################################################################
#
analysis.selectDataset <- function (ctx, name, select) {
    # Create merged metadata fields if needed
    sampleMeta <- ctx$unfiltered$meta
    if (!is.null(dataset$mergeFields)) {
         sampleMeta <- meta.addMergedFields(sampleMeta, dataset$mergeFields)
    }

    # Select the samples to be analyzed
    sampleMeta <- meta.select(sampleMeta, select)

    # Create a trimmed analysis context containing only data pertaining to the selected samples
    sampleNames <- rownames(sampleMeta)
    trimCtx <- analysis.selectSamplesInContext (ctx, sampleNames)
    
    sampleSet <- list(
        name=name,
        select=select,
        samples=sampleNames,
        ctx=trimCtx
    )
    ctx$sampleSets[[name]] <- sampleSet
    ctx
}
#
###################################################################
# Execute
###################################################################
#
analysis.executeOnSampleSet <- function(ctx, sampleSetName, tasks, plotList, aggregation, measures, params) {
    if (sampleSetName %in% names(ctx$sampleSets)) {
        stop(paste("Sample set not initiaized:", sampleSetName))
    }
    sampleSet <- ctx$sampleSets[[sampleSetName]]
    analysisName <- sampleSet$name
    print(paste("Analysis:", analysisName))
    
    # Get the trimmed analysis context containing only data pertaining to the selected samples
    sampleNames <- sampleSet$samples
    trimCtx <- sampleSet$ctx

    # Resolve all the automatic rendering in the plots
    if (!is.null(plotList)) {
        plotList <- resolveAutomaticRenderingInPlots (trimCtx$meta, plotList)
    }

    # Write out the meta
    metaFilename  <- paste(getOutFolder(analysisName), "/meta-", analysisName, "-unfiltered.tab", sep="")
    utils::write.table(trimCtx$unfiltered$meta, file=metaFilename, sep="\t", quote=FALSE, row.names=FALSE)
    metaFilename  <- paste(getOutFolder(analysisName), "/meta-", analysisName, "-filtered.tab", sep="")
    utils::write.table(trimCtx$filtered$meta, file=metaFilename, sep="\t", quote=FALSE, row.names=FALSE)

    # Execute computations and plots with different tasks
    for (tIdx in 1:length(tasks)) {
        tParts <- unlist(strsplit(tasks[tIdx], "/"))
        task <- tParts[1]
        method <- NULL
        if (length(tParts)>1) {
            method <- tParts[2]
        }
        print(paste(analysisName, "-", task, ifelse(is.null(method), "", method)))

        # Execute the task
        analysisCtx <- trimCtx$imputed
        if (task == "njt") {
            njt.execute (analysisCtx, analysisName)
            njt.executePlots (analysisCtx, analysisName, plotList)

        } else if (task == "pca") {
            pca.execute (analysisCtx, analysisName, method)
            pca.executePlots (analysisCtx, analysisName, method, plotList)

        } else if (task == "graph") {
            graph.execute (analysisCtx, analysisName, params)
            graph.executePlots (analysisCtx, analysisName, plotList, params)

        } else if (task == "haploNet") {
            haploNet.execute (analysisCtx, analysisName, plotList, params)

        } else if (task == "map") {
            if (method == "drug") {
                analysisCtx <- trimCtx$unfiltered
            }
            map.execute(analysisCtx, analysisName, method, aggregation, measures, params)

        } else {
            stop(paste("Invalid analysis task:", task))
        }
    }
    print(paste("Analysis", analysisName, "completed"))
}

###################################################################
# Task Parameters retrieval
###################################################################
#
#
#
analysis.getParam <- function (paramName, paramList, defaultValue) {
    value <- defaultValue
    if (!is.null(paramList)) {
        if (!is.null(paramList[[paramName]])) {
            value <- paramList[[paramName]]
        }
    }
    value
}
