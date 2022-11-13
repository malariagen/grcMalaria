###############################################################################
# Caching data files
################################################################################
distance.getDistanceDataFile <- function (ctx, datasetName) {
    dataFile <- getContextCacheFile(ctx, datasetName, "distance", "sampleDistance.tab")
    dataFile
}

distance.initialize <- function (ctx, datasetName, loadFromCache=TRUE, store=TRUE) {
    distDataFile <- distance.getDistanceDataFile(ctx, datasetName)
    dataset <- ctx[[datasetName]]
    if (loadFromCache & file.exists(distDataFile)) {
        distData <- readMatrix(distDataFile)
        distance.setDatasetDistance(ctx, datasetName, distData, store=FALSE)
        print(paste("Loaded", dataset$name, "distance matrix - Samples:", nrow(distData)))
    } else {
        genoData <- dataset$genos
        print(paste("Computing pairwise distances for",nrow(genoData),"samples using",ncol(genoData),"SNPs"))
        distData <- computeDistances (as.matrix(genoData))
        distance.setDatasetDistance(ctx, datasetName, distData, store=store)
    }
}

distance.setDatasetDistance <- function (ctx, datasetName, distance, store=TRUE) {
    dataset <- ctx[[datasetName]]
    dataset$distance <- distance
    if (store) {
        distDataFile <- distance.getDistanceDataFile(ctx, datasetName)
        writeMatrix(distance, distDataFile)
    }
}

