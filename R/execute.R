###################################################################
# Execute analyses
###################################################################
#
execute.executeOnSampleSet <- function(userCtx, sampleSetName, task, params) {
    #
    # Execute computations and plots 
    #
    execute.checkSampleSet (userCtx, sampleSetName)
    print(paste("Analyzing sample set:", sampleSetName))
    #
    # Split the task if there is a method part in it
    #
    print(paste(sampleSetName, "-", task))
    t <- execute.parseTask (task)				#; print(t)
    #
    # Execute the task
    #
    if (t$name == "pca") {
        pca.execute (userCtx, sampleSetName, t$method, params)
    } else if (t$name == "tree") {
        tree.execute (userCtx, sampleSetName, t$method, params)
    } else if (t$name == "graph") {
        clusterGraph.execute (userCtx, sampleSetName, params)
    } else if (t$name == "map") {
        map.execute (userCtx, sampleSetName, t$method, params)
    } else {
        stop(paste("Invalid analysis task:", t$name))
    }
    print(paste("Analysis", sampleSetName, "completed"))
}
#
execute.checkSampleSet <- function(userCtx, sampleSetName) {
    if (!sampleSetName %in% names(userCtx$sampleSets)) {
        stop(paste("Sample set not initialized:", sampleSetName))
    }
}
#
execute.parseTask <- function(task) {
    tParts <- unlist(strsplit(task, "/"))
    taskName <- tParts[1]
    method <- NULL
    if (length(tParts)>1) {
        method <- tParts[2]
    }
    list(name=taskName, method=method)
}
