###############################################################################
# Caching data files
################################################################################
distance.getCachedDistanceMatrixFilename <- function (ctx, datasetName) {
    distDataFile <- getContextCacheFile(ctx$rootCtx, datasetName, "distance", "sampleDistance")
    distDataFile
}

distance.initializeDistanceMatrix <- function (ctx, datasetName) {		#; print(paste("Initializing distance Matrix: ", datasetName))
    distDataFile <- distance.getCachedDistanceMatrixFilename(ctx, datasetName)
    if (!rdaFileExists(distDataFile)) {
        distData <- distance.createDistanceMatrix(ctx, datasetName)
    }
}

distance.retrieveDistanceMatrix <- function (ctx, datasetName) {
    distDataFile <- distance.getCachedDistanceMatrixFilename(ctx, datasetName)	#; print(distDataFile)
    if (rdaFileExists(distDataFile)) {						#; print("RDA file exists")
        distData <- readRdaMatrix(distDataFile)					#; print(distData[1:10,1:10])
    } else {
        distData <- distance.createDistanceMatrix(ctx, datasetName)
    }
    #
    # Trim the matrix  (eg. if this is a context from a sampleSet)
    #
    dataset <- ctx[[datasetName]]
    sampleNames <- rownames(dataset$meta)
    distData <- distData[sampleNames,sampleNames]
    distData
}

distance.createDistanceMatrix <- function (ctx, datasetName) {		#; print(paste("Creating distance Matrix: ", datasetName))
    dataset <- ctx$rootCtx[[datasetName]]
    genoData <- dataset$genos						#; print(genoData[1:10,30:40])
    print(paste("Computing pairwise distances for",nrow(genoData),"samples using",ncol(genoData),"SNPs"))
    genoMat <- as.matrix(genoData)					#; print(genoMat[1:10,1:10])
    distData <- computeDistances (genoMat)				#; print(head(distData))
    distDataFile <- distance.getCachedDistanceMatrixFilename(ctx, datasetName)
    writeRdaMatrix(distData, distDataFile)
    distData
}

DEFAULT_MOST_SIMILAR_COUNT <- 100
distance.findMostSimilarSamples <- function (ctx, datasetName, mostSimilarCount=DEFAULT_MOST_SIMILAR_COUNT) {
    distData <- distance.retrieveDistanceMatrix (ctx, datasetName)
    sampleNames <- colnames(distData)
    msIdxData  <- data.frame(matrix(nrow=mostSimilarCount,ncol=0))
    msDistData <- data.frame(matrix(nrow=mostSimilarCount,ncol=0))
    for (sIdx in 1:ncol(distData)) {
        sName <- sampleNames[sIdx]
        sDist <- distData[,sIdx]
        sOrder <- order(sDist)
        sOrder <- sOrder[which(sOrder!=sIdx)]
        sOrder <- sOrder[1:mostSimilarCount]
        msIdxData[,sName]  <- sOrder
        msDistData[,sName] <- sDist[sOrder]
    }
    list(indexes=msIdxData, distances=msDistData)
}