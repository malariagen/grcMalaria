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
            map.execute (userCtx, sampleSetName, method, params)

        } else {
            stop(paste("Invalid analysis task:", task))
        }
    }
    print(paste("Analysis", sampleSetName, "completed"))
}
