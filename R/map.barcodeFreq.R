###############################################################################
# Map Haplotype Frequency Analysis
################################################################################
#
barcodeFreqMap.executeMap <- function (map) {

    mapMaster   <- map$master
    mapType     <- mapMaster$type
    measureName <- map$measureName
    interval    <- map$interval
    #
    #datasetName <- map$datasetName
    datasetName <- "imputed"		# Cannot be another dataset, barcodes must be imputed.
    sampleSet   <- mapMaster$sampleSet
    userCtx     <- mapMaster$userCtx
    params      <- mapMaster$params
    config      <- context.getConfig(userCtx)
    #
    # Get the context and make sure we have samples to map
    #
    ctx        <- map$mapCtx				#; print(str(ctx))
    dataset    <- ctx[[datasetName]]
    sampleMeta <- context.getMeta (ctx, datasetName)
    if (nrow(sampleMeta)==0) {
        print(paste("No samples found - skipping interval", interval$name))
        return()
    }
    #
    plotTitle <- sampleSet$name
    baseMapInfo <- mapMaster$baseMapInfo
    #
    # Get the output folders
    #
    dataFolder <- mapMaster$dataFolder
    plotFolder <- mapMaster$plotFolder
    #
    # Silly trick to make the package checker happy... :-(
    lon <- lat <- label <- x <- y <- alpha <- NULL
    #
    # Parse the measure to get the visualization type and the minIdentity level
    #
    visType <- barcodeFreqMap.getVisualizationTypeFromMeasure (measureName)	#; print(visType)
    #
    # Now compute the aggregation units, the values to be plotted, and make the map
    # Get the aggregated data for the aggregation units
    #
    aggLevel <- as.integer(map$aggregation)     				#; print(aggLevel)   
    aggUnitData <- map$aggUnitData						#; print(aggUnitData)
    #
    # Work out a "standard" marker size (1/25 of the shortest side) and apply a scaling factor
    # Compute the number of samples that correspond to this standard size
    #
    stdMarkerSize <- clusterShareMap.getStandardMarkerSize(baseMapInfo, params)
    stdMarkerCount <- clusterShareMap.getStandardMarkerSampleCount(aggUnitData, params)
    #
    # Get the sample counts data
    #
    clusterCountData <- barcodeFreqMap.buildCountData (ctx, aggLevel, aggUnitData, params)	#; print(head(clusterCountData)
    #
    # Do the actual plot, starting with the background map
    #
    mapPlot <- baseMapInfo$baseMap
    #
    # If we need to show aggregation unit names, we need to compute the label positioning and plot
    #
    showMarkerNames <- param.getParam ("map.markerNames", params)
    if (showMarkerNames) {
        relSizes <- sqrt(aggUnitData$SampleCount)
        mapPlot <- map.addAggregationUnitNameLayer (mapPlot, aggUnitData, baseMapInfo, params, markerSize=relSizes)
    }
    #
    # Now add the markers and title
    #
    mapPlot <- barcodeFreqMap.addFreqMarkers (mapPlot, visType, clusterCountData, aggUnitData, stdMarkerCount, stdMarkerSize, params)
    mapPlot <- mapPlot + ggplot2::labs(title=plotTitle, subtitle="")
    #
    # Return map plot for completion and saving to file
    #
    mapPlot
}
#
################################################################################
#
# We use visualization types to specify different graphic renditions from the same analyses: "pie" or "bar" markers
#
barcodeFreqMap.resolveMeasureNames <- function(clusterSets, params) {		#; print(names(params))
    visTypes <- param.getParam ("map.cluster.visualizations", params)		#; print(visTypes) # "pie" or "bar"
    measureLabels <- paste("barcodeFrequency", visTypes, sep="/")		#; print(measureLabels)
    measureLabels
}
#
################################################################################
#
barcodeFreqMap.getVisualizationTypeFromMeasure <- function(measureName) {		#; print(measure)
    sParts <- unlist(strsplit(measureName, "/"))
    visType=sParts[2]
    visType
}
#
###############################################################################
# Haplotype Frequency Markers plotting 
################################################################################
#
barcodeFreqMap.addFreqMarkers <- function (mapPlot, visType, countData, aggUnitData, stdMarkerCount, stdMarkerSize, params) {
    if (visType=="bar") {
        mapPlot <- barcodeFreqMap.addFreqBars (mapPlot, countData, aggUnitData, stdMarkerCount, stdMarkerSize, params)
    } else if (visType=="pie") {
        mapPlot <- barcodeFreqMap.addFreqPies (mapPlot, countData, aggUnitData, stdMarkerCount, stdMarkerSize, params)
    }
    mapPlot
}

barcodeFreqMap.addFreqPies <- function (mapPlot, countData, aggUnitData, stdMarkerCount, stdMarkerSize, params) {
    # Silly trick to make the package checker happy... :-(
    Longitude <- Latitude <- Haplo <- ClusterCount <- SampleCount <-NULL
    # Nasty trick we have to do because there is no linewidth aesthetic at present
    ggplot2::update_geom_defaults(ggforce::GeomArcBar, ggplot2::aes(linewidth=!!params$pieLineWidth))

    # Now add the pie chart markers
    if (stdMarkerCount==0) {
        mapPlot <- mapPlot +
                   ggforce::geom_arc_bar(ggplot2::aes(x0=Longitude, y0=Latitude, r0=0, r=stdMarkerSize, 
                                         fill=Haplo, amount=ClusterCount),
                             data=countData, stat="pie", inherit.aes=FALSE,
                             #colour="gray25", stroke=0.5, fill="white", show.legend=FALSE)
                             colour="gray25", fill="white", show.legend=FALSE)
    } else {
        mapPlot <- mapPlot +
                   ggforce::geom_arc_bar(ggplot2::aes(x0=Longitude, y0=Latitude, r0=0, r=stdMarkerSize*sqrt(SampleCount/stdMarkerCount), 
                                         fill=Haplo, amount=ClusterCount),
                             data=countData, stat="pie", inherit.aes=FALSE,
                             #colour="gray25", stroke=0.5, fill="white", show.legend=FALSE)
                             colour="gray25", fill="white", show.legend=FALSE)
    }
    mapPlot
}

barcodeFreqMap.addFreqBars <- function (mapPlot, countData, aggUnitData, stdMarkerCount, stdMarkerSize, params) {

    # Get all aggregation unit ids, in descending order of sample count
    aggUnitData <- aggUnitData[order(-aggUnitData$SampleCount),]
    aggUnitGids <- rownames(aggUnitData)							#; print(aggUnitGids)
    freqBarData <- NULL
    for (aIdx in 1:length(aggUnitGids)) {
        aggUnitGid <- aggUnitGids[aIdx]
        sampleCount <- aggUnitData$SampleCount[aIdx]
        
        # Determine the desired size of the marker
        if (stdMarkerCount==0) {
            mHeight <- 2 * stdMarkerSize
        } else {
            mHeight <- 2 * stdMarkerSize*sqrt(sampleCount/stdMarkerCount)
        }
        mWidth <- mHeight / 2 
        
        # Determine the height representing one sample
        sHeight <- mHeight / sampleCount
        
        # Determine the starting coords (lower left corner)
        y0 <- aggUnitData$Latitude[aIdx] - (mHeight/2)
        x1 <- aggUnitData$Longitude[aIdx] - (mWidth/2)
        x2 <- aggUnitData$Longitude[aIdx] + (mWidth/2)
        
        # Get all the clusters for this unit, ordered by sample count
        unitClusterData <- countData[which(countData$UnitId==aggUnitGid),]
        unitClusterData <- unitClusterData[order(unitClusterData$ClusterCount),]
        unitRowCount <- nrow(unitClusterData)
        
        # For each cluster, work out the y boundaries
        y2 <- vector(mode="numeric", length=unitRowCount)
        y <- y0
        for (i in 1:unitRowCount) {
           y <- y + (unitClusterData$ClusterCount[i] * sHeight)
           y2[i] <- y
        }
        y1 <- c(y0, y2[1:(unitRowCount-1)])
        x1 <- rep(x1, unitRowCount)
        x2 <- rep(x2, unitRowCount)
        unitBarData <- cbind(unitClusterData, x1=x1, x2=x2, y1=y1, y2=y2)
        freqBarData <- rbind(freqBarData, unitBarData)
    }				
    #print(freqBarData)
    mapPlot <- mapPlot +
               ggplot2::geom_rect(ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), 
                                  data=freqBarData, inherit.aes=FALSE,
                                  colour="gray25", linewidth=params$pieLineWidth, fill="white", show.legend=FALSE)
    mapPlot
}

#
# For each aggregation unit, we get a count of each unique cluster, ordered in descending count
#
barcodeFreqMap.buildCountData <- function(ctx, aggLevel, aggUnitData, params) {

    # Get all aggregation unit ids
    aggUnitGids <- rownames(aggUnitData)						#; print(aggUnitGids)
    
    # Create aggregation index for each sample (the id of the aggregation unit where the sample originates)
    sampleMeta <- context.getMeta (ctx, "imputed")
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta)
    clusters <- ctx$rootCtx$imputed$barcodes

    # Get the data for all aggregation units
    countData <- NULL
    for (aIdx in 1:length(aggUnitGids)) {
        # Get the sample data to be aggregated for this unit
        aggUnitGid <- aggUnitGids[aIdx]
        unitCusters <- clusters[which(aggIndex == aggUnitGid)]				#; print(nrow(aggHaplos))
        unitClusterCounts <- as.integer(table(unitCusters))
        unitClusterCounts <- unitClusterCounts[order(-unitClusterCounts)]
        unitHaploNum <- length (unitClusterCounts)

        unitData <- aggUnitData[aIdx,]
        df <- data.frame(matrix(unitData,ncol=ncol(aggUnitData),nrow=unitHaploNum,byrow=TRUE))
        colnames(df) <- colnames(aggUnitData)
        df$Haplo <- paste("Haplo",c(1:unitHaploNum),sep="_")
        df$ClusterCount <- unitClusterCounts
        countData <- rbind(countData, df)
    }
    countData$Haplo <- factor(countData$Haplo)
    countData$Latitude <- as.numeric(countData$Latitude)
    countData$Longitude <- as.numeric(countData$Longitude)
    countData$SampleCount <- as.numeric(countData$SampleCount)
    countData
}
