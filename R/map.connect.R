###############################################################################
# Measures of connectedness
###############################################################################
#
connectMap.getConnectednessMeasures <- function() {
    allMeasureNames <- c("similarity","meanDistance")
    allMeasureNames
}

###############################################################################
# Map Aggregated Measure Analysis
################################################################################
#
connectMap.executeMap <- function(map) {

    mapMaster   <- map$master
    mapType     <- mapMaster$type
    measureName <- map$measureName
    interval    <- map$interval
    #
    datasetName <- map$datasetName
    sampleSet   <- mapMaster$sampleSet
    userCtx     <- mapMaster$userCtx
    params      <- mapMaster$params
    config      <- userCtx$config
    #
    # Get the context, trimmed by time interval
    #
    ctx        <- map$mapCtx				#; print(str(ctx))
    dataset    <- ctx[[datasetName]]
    sampleMeta <- dataset$meta
    if (nrow(sampleMeta)==0) {
        print(paste("No samples found - skipping interval", interval$name))
        return()
    }
    #
    baseMapInfo <- mapMaster$baseMapInfo
    #
    # Get the output folders
    #
    dataFolder <- mapMaster$dataFolder
    plotFolder <- mapMaster$plotFolder
    #
    # Silly trick to make the package checker happy... :-(
    lon <- lat <- Latitude <- Longitude <- Lon1 <- Lon2 <- Lat1 <- Lat2 <- label <- NULL
    #
    # Now compute the aggregation units, the values to be plotted, and make the map
    # Get the aggregated data for the aggregation units
    # Get the connectedness measure data for the aggregation unit pairs
    #
    aggLevel <- as.integer(map$aggregation)     				#; print(aggLevel)   
    aggUnitData <- map$aggUnitData						#; print(aggUnitData)
    aggUnitPairData <- map$aggUnitPairData					#; print(aggUnitPairData)
    #
    #
    #
    mValues <- aggUnitPairData[,measureName]					#; print(measureName); 
    selAggUnitPairData <- aggUnitPairData					#; print(nrow(selAggUnitPairData))
    #
    # If we're basing connectedness on mean distance, get the threshold below which we must remove the pairs
    # and select the aggregation unit pairs to be plotted
    #
    minValue <- 0
    if (startsWith(measureName, "meanDistance")) {
        minValue <- connectMap.getMeasureLevel (measureName, "meanDistance-ge")	#; print(minValue)
        selAggUnitPairData <- aggUnitPairData[which(mValues > minValue),]	#; print(nrow(selAggUnitPairData))
    }
    #
    # Sort the pairs so that higher values get plotted last
    #
    mValues <- selAggUnitPairData[,measureName]
    selAggUnitPairData <- selAggUnitPairData[order(mValues),]
    #
    # Do the actual plot, starting with the background map
    #
    mapPlot <- baseMapInfo$baseMap
    #
    # Now plot the connections
    #
    mapPlot <- mapPlot +
        ggplot2::geom_curve(ggplot2::aes(x=Lon1, y=Lat1, xend=Lon2, yend=Lat2, 
                                         linewidth=!!rlang::sym(measureName), colour=!!rlang::sym(measureName)),
                            data=selAggUnitPairData, 
                            curvature=0.2, alpha=0.75) +
        ggplot2::scale_linewidth_continuous(guide="none", range=c(params$connectCurveWidthMin, params$connectCurveWidthMax)) +          # scale for edge widths
        ggplot2::scale_colour_gradientn(colours=c("skyblue1","midnightblue"))
    #
    # Now add the markers
    #
    mapPlot <- mapPlot +
        ggplot2::geom_point(ggplot2::aes(x=Longitude, y=Latitude), 
                            data=aggUnitData,
                            size=params$map.markerSize, shape=19, col="red")
    #
    # If we need to show aggregation unit names, we need to compute the label positioning and plot
    #
    showMarkerNames <- param.getParam ("map.markerNames", params)
    if (showMarkerNames) {
        mapPlot <- map.addAggregationUnitNameLayer (mapPlot, aggUnitData, baseMapInfo, params)
    }
    #
    # Return map plot for completion and saving to file
    #
    mapPlot
}

connectMap.resolveMeasureNames <- function(params) {
    measureNames <- param.getParam ("analysis.measures", params)		#; print(measureNames)
    if ("ALL" %in% measureNames) {
        measureNames <- connectMap.getConnectednessMeasures()
    }
    #
    # Create one "expanded" measure for each similarity threshold (e.g. "similarity-ge0.80")
    #
    expanded <- c()
    for (mIdx in 1:length(measureNames)) {
        measureName <- measureNames[mIdx]
        if (measureName == "similarity") {
            levels <- as.numeric(param.getParam("map.connect.identity.min", params))
            prefix <- "similarity-ge"
        } else if (measureName == "meanDistance") {
            levels <- as.numeric(param.getParam("map.connect.meanDistance.min", params))
            prefix <- "meanDistance-ge"
        } else {
            stop(paste("Invalid measure of connectedness:",measureName))
        }
        for (lIdx in 1 : length(levels)) {
            measureStr <- paste(prefix, format(levels[lIdx], digits=2, nsmall=2), sep="")
            expanded <- c(expanded, measureStr)
        }
    }									#; print(expanded)
    expanded
}

connectMap.estimateMeasures <- function (ctx, datasetName, sampleSetName, aggLevel, aggUnitData, mapType, measureNames, params, dataFolder)	{
    dataset <- ctx[[datasetName]]
    sampleMeta   <- dataset$meta
    barcodeData  <- dataset$barcodes
    distData     <- distance.retrieveDistanceMatrix (ctx, datasetName)

    # Create aggregation index for each sample (the id of the aggregation unit where the sample originates)
    sampleGids <- as.character(map.getAggregationUnitIds (aggLevel, sampleMeta))
    sampleIds <- rownames(sampleMeta)
    
    # Get all aggregation units
    aggUnitGids <- rownames(aggUnitData)						#; print(aggUnitGids) ; print(head(aggUnitData))
    
    # Get the data for all aggregation units
    measureData <- NULL
    for (a1Idx in 1:(length(aggUnitGids)-1)) {
        for (a2Idx in (a1Idx+1):length(aggUnitGids)) {
            # Get the sample data
            gid1 <- aggUnitGids[a1Idx]							#; print(gid1)
            name1 <- aggUnitData$AdmDivName[a1Idx]					#; print(name1)
            lat1 <- aggUnitData$Latitude[a1Idx]
            lon1 <- aggUnitData$Longitude[a1Idx]
            samples1 <- sampleIds[which(sampleGids == gid1)]				#; print(length(samples1))
            #
            gid2 <- aggUnitGids[a2Idx]							#; print(gid2)
            name2 <- aggUnitData$AdmDivName[a2Idx]					#; print(name2)
            lat2 <- aggUnitData$Latitude[a2Idx]
            lon2 <- aggUnitData$Longitude[a2Idx]
            samples2 <- sampleIds[which(sampleGids == gid2)]				#; print(length(samples2))
            #
            pairDist <- distData[samples1,samples2]					#; print(dim(pairDist))
            #
            # Get the admin division values from the first sample of this unit (assuming the values are the same for all)
            mValues <- connectMap.estimateDistMeasures (pairDist, measureNames)		#; print (mValues)
            cValues <- c(gid1, gid2, name1, name2, lat1, lon1, lat2, lon2, mValues)
            measureData <- rbind(measureData, cValues)
        }
    }
    aggUnitPairData <- data.frame(measureData)						#; print (dim(aggUnitPairData))
    nc <- ncol(aggUnitPairData)
    aggUnitPairData[,5:nc] <- sapply(aggUnitPairData[,5:nc], as.numeric)		#; print(ncol(aggUnitPairData))
    colnames(aggUnitPairData) <- c("Unit1", "Unit2", "UnitName1", "UnitName2", "Lat1", "Lon1", "Lat2", "Lon2", measureNames)
    
    # Write out the aggregation unit data to file
    aggDataFilename  <- paste(dataFolder, "/AggregatedPairData-", sampleSetName, "-", aggLevel, ".tab", sep="")
    utils::write.table(aggUnitPairData, file=aggDataFilename, sep="\t", quote=FALSE, row.names=FALSE)

    aggUnitPairData
}
#
connectMap.estimateDistMeasures <- function (pairDistData, measureNames) {
    result <- c()
    dist <- unlist(pairDistData)
    for (mIdx in 1:length(measureNames)) {
        measureName <- measureNames[mIdx]
        if (startsWith(measureName, "similarity")) {
            level <- connectMap.getMeasureLevel (measureName, "similarity-ge")	#; print(level)
            maxDist = 1.0 - level						#; print(maxDist)
            value <- length(which(dist <= maxDist)) / length(dist)
        } else if (startsWith(measureName, "meanDistance")) {
            value <- 1 - mean(dist)
        } else {
            stop(paste("Invalid distance connectedness measure:", measureName))
        }
        result <- c(result, value)
    }
    result
}
#
connectMap.getMeasureLevel <- function(measureName, prefix) {
    if (!startsWith(measureName, prefix)) {
        return (NA)				#; print("Incorrect prefix")
    }
    level <- substring(measureName,nchar(prefix)+1)	#; print(level)
    as.numeric(level)
}

