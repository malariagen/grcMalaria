###################################################################
# Root Context - Contains all data needed from analysis
###################################################################
#
analysis.createContext <- function (grcData, config) {
    ctx <- list()
    ctx$unfiltered <- analysis.createUnfilteredContext(grcData, config)
    ctx$filtered   <- analysis.createFilteredContext(ctx$unfiltered)
    ctx$imputed    <- impute.createContext(ctx$filtered)
    ctx$config     <- config
    ctx$sampleSets <- list()
    ctx
}
#
###################################################################
# Basic (Unfiltered) Data - complete but raw data - built upon initialization
###################################################################
#
analysis.createUnfilteredContext <- function (grcData, config) {
    print("Initializing Basic Dataset")
    unfilteredCtx <- list(name="unfiltered")
    unfilteredCtx <- meta.setContextSampleMeta(ctx, grcData, store=TRUE)
    print(paste("Loaded metadata - Samples:", nrow(unfilteredCtx$meta)))
    unfilteredCtx
}
#
###################################################################
# Filtered Data, with high-missingness barcodes removed
###################################################################
#
analysis.createFilteredContext <- function (unfilteredCtx, loadFromCache=TRUE) {
    print("Initializing Filtered Dataset")
    filteredCtx <- list(name="filtered")

    filteredMetaFile        <- meta.getMetaDataFile(filteredCtx)
    filteredBarcodeMetaFile <- barcode.getBarcodeMetaFile(filteredCtx)
    filteredBarcodeFile     <- barcode.getBarcodeDataFile(filteredCtx)

    if (loadFromCache & file.exists(filteredMetaFile) & file.exists(filteredBarcodeMetaFile) & file.exists(filteredBarcodeFile)) {
        sampleMeta <- readSampleData (filteredMetaFile)		#; print(colnames(sampleMeta))
        filteredCtx <- meta.setContextSampleMeta (filteredCtx, sampleMeta, store=FALSE)

        barcodeMeta <- read.table(filteredBarcodeMetaFile, as.is=TRUE, header=TRUE, sep="\t")
        barcodeData <- readSampleData (filteredBarcodeFile)
        filteredCtx <- setContextBarcodes (filteredCtx, barcodeData, barcodeMeta, store=FALSE)
        print(paste("Loaded", filteredCtx$name, "barcodes - Samples:", nrow(barcodeData), "x SNPs:", ncol(barcodeData)))
    } else {
        filteredCtx <- meta.setContextSampleMeta (filteredCtx, unfilteredCtx$meta)
        filteredCtx <- initializeBarcodes (filteredCtx)

        # Trim the metadata to cover the barcodes selected
        sampleNames   <- rownames(filteredCtx$barcodes)
        filteredMeta  <- filterSampleData (filteredCtx$meta, sampleNames)
        filteredCtx <- meta.setContextSampleMeta (filteredCtx, filteredMeta)
    }

    # Get the genotypes, distance matrix and execute the pop structure analysis
    filteredCtx <- geno.initialize(filteredCtx)
    filteredCtx <- distance.initialize(filteredCtx)
    filteredCtx
}
#
###################################################################
# Create a new context, consisting of the same data as the source context, but filtered to a given sample list
###################################################################
#
analysis.selectSamplesInContext <- function (ctx, sampleNames) {
    trimCtx <- list()
    trimCtx$unfiltered <- analysis.trimContext (ctx$unfiltered, sampleNames)
    if (!is.null(ctx$filtered)) {
        trimCtx$filtered   <- analysis.trimContext (ctx$filtered, sampleNames)
    }
    if (!is.null(ctx$imputed)) {
        trimCtx$imputed    <- analysis.trimContext (ctx$imputed, sampleNames)
    }
    trimCtx
}
#
analysis.trimContext <- function (ctx, sampleNames) {
    # Use only samples that are shared with the original context
    ctxSampleNames <- rownames(ctx$meta)
    sampleNames <- sampleNames[which(sampleNames %in% ctxSampleNames)]
    #
    trimCtx <- list(name=ctx$name)
    trimCtx <- meta.setContextSampleMeta (trimCtx, ctx$meta[sampleNames,],store=FALSE)
    if (!is.null(ctx$barcodes)) {
        trimCtx <- setContextBarcodes (trimCtx, ctx$barcodes[sampleNames,], ctx$barcodeMeta, store=FALSE)
    }
    if (!is.null(ctx$genos)) {
        trimCtx <- geno.setContextGenotypes (trimCtx, ctx$genos[sampleNames,], ctx$genoMeta, store=FALSE)
    }
    if (!is.null(ctx$distance)) {
        trimCtx <- distance.setContextDistance (trimCtx, ctx$distance[sampleNames,sampleNames], store=FALSE)
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
    write.table(trimCtx$unfiltered$meta, file=metaFilename, sep="\t", quote=FALSE, row.names=FALSE)
    metaFilename  <- paste(getOutFolder(analysisName), "/meta-", analysisName, "-filtered.tab", sep="")
    write.table(trimCtx$filtered$meta, file=metaFilename, sep="\t", quote=FALSE, row.names=FALSE)

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
