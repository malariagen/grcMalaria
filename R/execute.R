###################################################################
# Execute analyses
###################################################################
#
execute.executeOnSampleSet <- function(userCtx, sampleSetName, tasks, params) {
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
