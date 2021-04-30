source("graphics.R",     local=FALSE)
source("pca.R",          local=FALSE)
source("njt.R",          local=FALSE)
source("cluster.R",      local=FALSE)
source("graph.R",        local=FALSE)
source("haploNet.R",     local=FALSE)
source("map.common.R",   local=FALSE)
source("map.marker.R",   local=FALSE)
source("map.connect.R",  local=FALSE)
source("map.haplo.R",    local=FALSE)
#
###################################################################
# Root Context - Contains all data needed from analysis
###################################################################
#
analysis.createContext <- function () {
    ctx <- list()
    ctx$unfiltered <- analysis.createUnfilteredContext()
    ctx$filtered   <- analysis.createFilteredContext(ctx$unfiltered)
    ctx$imputed    <- impute.createContext(ctx$filtered)
    ctx
}
#
# Datafile Naming utility- so that datafiles are not overwritten 
#
analysis.getDataFile <- function (context, dataFolder, filename) {
    return(paste(dataFolder, "/", context$name, ".", filename, sep=""))
}
#
###################################################################
# Basic (Unfiltered) Data - complete but raw data - built upon initialization
###################################################################
#
analysis.createUnfilteredContext <- function () {
    print("Initializing Basic Dataset")
    unfilteredCtx <- list(name="unfiltered")
    unfilteredCtx <- meta.initializeMetadata (unfilteredCtx)
    print(paste("Loaded metadata - Samples:", nrow(unfilteredCtx$meta)))
    unfilteredCtx
}
#
###################################################################
# Filtered Data, with high-missingness barcodes removed
###################################################################
#
analysis.createFilteredContext <- function (unfilteredContext, loadFromCache=TRUE) {
    print("Initializing Filtered Dataset")
    filteredContext <- list(name="filtered")
    
    filteredMetaFile <- analysis.getDataFile(filteredContext, folder.data.meta, sampleMetaFname)
    filteredBarcodeMetaFile <- analysis.getDataFile(filteredContext, folder.data.barcode, barcodeMetaFname)
    filteredBarcodeFile <- analysis.getDataFile(filteredContext, folder.data.barcode, barcodeFname)
    
    if (loadFromCache & file.exists(filteredMetaFile) & file.exists(filteredBarcodeMetaFile) & file.exists(filteredBarcodeFile)) {
        sampleMeta <- readSampleData (filteredMetaFile)		#; print(colnames(sampleMeta))
        filteredContext <- meta.setContextSampleMeta (filteredContext, sampleMeta, store=FALSE)
        
        barcodeMeta <- read.table(filteredBarcodeMetaFile, as.is=TRUE, header=TRUE, sep="\t")
        barcodeData <- readSampleData (filteredBarcodeFile)
        filteredContext <- setContextBarcodes (filteredContext, barcodeData, barcodeMeta, store=FALSE)
        print(paste("Loaded", filteredContext$name, "barcodes - Samples:", nrow(barcodeData), "x SNPs:", ncol(barcodeData)))
    } else {
        filteredContext <- meta.setContextSampleMeta (filteredContext, unfilteredContext$meta)
        filteredContext <- initializeBarcodes (filteredContext)
    
        # Trim the metadata to cover the barcodes selected
        sampleNames   <- rownames(filteredContext$barcodes)
        filteredMeta  <- filterSampleData (filteredContext$meta, sampleNames)
        filteredContext <- meta.setContextSampleMeta (filteredContext, filteredMeta)
    }

    # Get the genotypes, distance matrix and execute the pop structure analysis
    filteredContext <- geno.initialize(filteredContext)
    filteredContext <- distance.initialize(filteredContext)
    filteredContext
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
# Execute
###################################################################
#
analysis.execute <- function (ctx, datasetList, tasks, plotList=NULL, aggregation=NULL, measures=NULL, params=NULL) {
    for (aIdx in 1:length(datasetList)) {
        dataset <- datasetList[[aIdx]]
        analysis.executeSingle (ctx, dataset=dataset, tasks=tasks, plotList=plotList, aggregation=aggregation, measures=measures, params=params) 
    }
    print("All datasets analyzed")
}
#
###################################################################
# Run Individual Analysis
###################################################################
#
analysis.executeSingle <- function(ctx, dataset, tasks, plotList, aggregation, measures, params) {
    analysisName <- dataset$name
    print(paste("Analysis:", analysisName))
  
    # Create merged metadata fields if needed
    sampleMeta <- ctx$unfiltered$meta
    if (!is.null(dataset$mergeFields)) {
         sampleMeta <- meta.addMergedFields(sampleMeta, dataset$mergeFields)
    }
  
    # Select the samples to be plotted
    if (!is.null(dataset$select)) {
        sampleMeta <- meta.select(sampleMeta, dataset$select)
    }
    
    # Resolve all the automatic rendering in the plots
    if (!is.null(plotList)) {
        plotList <- resolveAutomaticRenderingInPlots (sampleMeta, plotList)
    }
    
    # Create a trimmed analysis context containing only data pertaining to the selected samples
    sampleNames <- rownames(sampleMeta)
    trimCtx <- analysis.selectSamplesInContext (ctx, sampleNames)
  
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
