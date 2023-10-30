###################################################################
# 
###################################################################
#
sampleSet.selectSampleSet <- function (userCtx, sampleSetName, select) {
    sampleMeta <- userCtx$unfiltered$meta

    # Select the samples to be analyzed
    sampleMeta <- meta.select(sampleMeta, select)

    # Create a trimmed analysis context containing only data pertaining to the selected samples
    sampleNames <- rownames(sampleMeta)
    trimCtx <- context.trimContext (userCtx, sampleNames)
    trimCtx$sampleSets <- NULL
    
    sampleSet <- new.env()
    sampleSet$name=sampleSetName
    sampleSet$select=select
    sampleSet$samples=sampleNames
    sampleSet$ctx=trimCtx
    sampleSet$clusters=new.env()
    
    userCtx$sampleSets[[sampleSetName]] <- sampleSet
    
    metaOutFolder <- getOutFolder(userCtx$config, c(sampleSetName, "metadata"))
    metaFilename  <- paste(metaOutFolder, "/meta-", sampleSetName, "-unfiltered.tab", sep="")
    utils::write.table(trimCtx$unfiltered$meta, file=metaFilename, sep="\t", quote=FALSE, row.names=FALSE)
    metaFilename  <- paste(metaOutFolder, "/meta-", sampleSetName, "-filtered.tab", sep="")
    utils::write.table(trimCtx$filtered$meta, file=metaFilename, sep="\t", quote=FALSE, row.names=FALSE)

    unfilteredCount <- length(sampleNames)
    filteredCount <- nrow(trimCtx$filtered$meta)
    print(paste0("Selected ", unfilteredCount, " samples for dataset '", sampleSetName, "', including ", filteredCount, " quality filtered samples"))
}
