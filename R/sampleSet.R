###################################################################
# 
###################################################################
#
sampleSet.selectSampleSet <- function (ctx, sampleSetName, select) {

    # Select the samples to be analyzed
    rootCtx <- ctx$rootCtx
    sampleMeta <- context.getMeta (rootCtx) 
    selMeta <- meta.select(sampleMeta, select)
    selSamples <- rownames(selMeta)

    # Create a trimmed analysis context containing only data pertaining to the selected samples
    trimCtx <- context.trimContext (rootCtx, selSamples)
    trimCtx$sampleSets <- NULL
    
    sampleSet <- new.env()
    sampleSet$name=sampleSetName
    sampleSet$select=select
    sampleSet$samples=selSamples
    sampleSet$ctx=trimCtx
    sampleSet$clusters=new.env()
    sampleSet$baseMaps=new.env()
    
    rootCtx$sampleSets[[sampleSetName]] <- sampleSet
    
    unfilteredCount <- length(trimCtx$unfiltered$samples)
    filteredCount   <- length(trimCtx$filtered$samples)
    imputedCount    <- length(trimCtx$imputed$samples)
    print(paste0("Selected ", unfilteredCount, " samples for sampleset '", sampleSetName, "', including ", 
                              filteredCount, " quality filtered samples and ", imputedCount, " imputed samples"))
}
