###############################################################################
# Map Cluster Sharing Analysis
################################################################################
#
CLUSTER_SHARE_MAP_PREFIX <- "clusterSharing-ge"
#
################################################################################
#
clusterShareMap.executeMap <- function (map) {

    mapMaster   <- map$master
    mapType     <- mapMaster$type
    measure     <- map$measure
    interval    <- map$interval
    #
    pp          <- mapMaster$plotParams
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
    mParts <- clusterShareMap.parseMeasure (measure)			#; print(measure)
    minIdentity <- mParts$minIdentity					#; print(minIdentity)
    visType <- mParts$visType						#; print(visType)		
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
    # Get the cluster data
    #
    clusterSetLabel  <- getMinIdentityLabel (minIdentity)
    clusterSet       <- mapMaster$clusterSets[[clusterSetLabel]]
    clusterSetName   <- clusterSet$clusterSetName
    clusterData      <- clusterSet$clusters
    #
    clusterShareData <- clusterMap.getUnitClusterCountData (clusterData, sampleMeta, aggLevel, aggUnitData, params)    #TODO - rename to aggClusterCountData
											#; print(clusterShareData)
    clusterPalette <- mapMaster$clusterSetPalettes[[clusterSetLabel]]
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
        mapPlot <- map.addAggregationUnitNameLayer (mapPlot, aggUnitData, baseMapInfo, pp, markerSize=relSizes)
    }
    #
    # Now add the markers and title
    #
    mapPlot <- clusterShareMap.addShareMarkers (mapPlot, visType, clusterShareData, clusterPalette, aggUnitData, stdMarkerCount, stdMarkerSize, pp)
    mapPlot <- mapPlot + ggplot2::labs(title=plotTitle, subtitle="")
    #
    # Return map plot for completion and saving to file
    #
    mapPlot
}
#
################################################################################
#
# Work out a "standard" marker size (1/25 of the shortest side) and apply a scaling factor
#
clusterShareMap.getStandardMarkerSize <- function(baseMapInfo, params) {
    #
    # Work out a "standard" marker size (1/25 of the shortest side) and apply a scaling factor
    #
    scalingFactor <- param.getParam ("map.cluster.markerScale", params)			#; print(scalingFactor)
    bbox <- baseMapInfo$anBB;
    minSide <- min(abs(bbox$yMax-bbox$yMin),abs(bbox$xMax-bbox$xMin))
    stdMarkerSize <- scalingFactor * minSide / 25					#; print(stdMarkerSize)
    stdMarkerSize
}    
#
################################################################################
#
# Compute the number of samples that correspond to the standard marker size
#
clusterShareMap.getStandardMarkerSampleCount <- function(aggUnitData, params) {
    #
    # Compute the number of samples that correspond to this standard size
    #
    stdMarkerCountParam <- param.getParam ("map.cluster.markerSampleCount", params)	#; print(stdMarkerCountParam)
    stdMarkerCount <- 0
    if (is.numeric(stdMarkerCountParam)) {
        stdMarkerCount <- stdMarkerCountParam
    } else if (stdMarkerCountParam == "mean") {
        stdMarkerCount <- mean(as.numeric(aggUnitData$SampleCount))			#; print(aggUnitData$SampleCount)
    } else if (stdMarkerCountParam == "none") {
        stdMarkerCount <- 0
    } else {
        stop(paste("Invalid map.stdMarkerCount value:", stdMarkerCountParam))        
    }											#; print(stdPieCount)
    stdMarkerCount <- as.numeric(stdMarkerCount)					#; print(stdMarkerCount)
    stdMarkerCount
}    
#
################################################################################
#
# We use visualization types to specify different graphic renditions from the same analyses: "pie" or "bar" markers
#
clusterShareMap.resolveMeasures <- function(clusterSets, params) {
    expanded <- c()
    visTypes <- param.getParam ("map.cluster.visualizations", params)			#; print(visTypes) # "pie" or "bar"
    minIdentityLabels <- names(clusterSets)
    setCount <- length(minIdentityLabels)
    for (idIdx in 1:setCount) {
        minIdentityLabel <- minIdentityLabels[idIdx]
        minIdentity <- getMinIdentityFromLabel (minIdentityLabel)
        measureLabels <- paste(getMinIdentityLabel(minIdentity, prefix=CLUSTER_SHARE_MAP_PREFIX), 
                               visTypes, sep="/")				#; print(measureLabels)
        expanded <- c(expanded, measureLabels)
    }										#; print(expanded)
    expanded
}
#
################################################################################
#
clusterShareMap.parseMeasure <- function(measure) {
    prefix <- CLUSTER_SHARE_MAP_PREFIX
    if (!startsWith(measure, prefix)) {
        return (NULL)				#; print("Incorrect prefix")
    }
    suffix <- substring(measure,nchar(prefix)+1)
    sParts <- unlist(strsplit(suffix, "/"))
    list(minIdentity=as.numeric(sParts[1]), visType=sParts[2])
}
#
###############################################################################
# Haplotype Sharing Markers plotting 
################################################################################
#
clusterShareMap.addShareMarkers <- function (mapPlot, visType, clusterShareData, clusterPalette, aggUnitData, stdMarkerCount, stdMarkerSize, pp) {
    # Do th actual marker plotting for the twp types of markers
    if (visType=="bar") {
        mapPlot <- clusterShareMap.addShareBars (mapPlot, clusterShareData, clusterPalette, aggUnitData, stdMarkerCount, stdMarkerSize, pp)
    } else if (visType=="pie") {
        mapPlot <- clusterShareMap.addSharePies (mapPlot, clusterShareData, clusterPalette, aggUnitData, stdMarkerCount, stdMarkerSize, pp)
    }
    mapPlot
}
#
#
#
clusterShareMap.addSharePies <- function (mapPlot, clusterShareData, clusterPalette, aggUnitData, stdMarkerCount, stdMarkerSize, pp) {
    # Silly trick to make the package checker happy... :-(
    Longitude <- Latitude <- Cluster <- ClusterCount <- SampleCount <-NULL
    # Nasty trick we have to do because there is no linewidth aesthetic at present
    ggplot2::update_geom_defaults(ggforce::GeomArcBar, ggplot2::aes(linewidth=!!pp$pieLineWidth))

    # Now add the pie chart markers
    if (stdMarkerCount==0) {
        mapPlot <- mapPlot +
                   ggforce::geom_arc_bar(ggplot2::aes(x0=Longitude, y0=Latitude, r0=0, r=stdMarkerSize, 
                                         fill=Cluster, amount=ClusterCount),
                                data=clusterShareData, stat="pie", inherit.aes=FALSE,
                                colour="gray25", show.legend=FALSE) +
                   ggplot2::scale_fill_manual(values=clusterPalette)
    } else {
        mapPlot <- mapPlot +
                   ggforce::geom_arc_bar(ggplot2::aes(x0=Longitude, y0=Latitude, r0=0, r=stdMarkerSize*sqrt(SampleCount/stdMarkerCount), 
                                         fill=Cluster, amount=ClusterCount),
                                data=clusterShareData, stat="pie", inherit.aes=FALSE,
                                colour="gray25", show.legend=FALSE) +
                   ggplot2::scale_fill_manual(values=clusterPalette)
    }
    mapPlot
}
#
#
#
clusterShareMap.addShareBars <- function (mapPlot, clusterShareData, clusterPalette, aggUnitData, stdMarkerCount, stdMarkerSize, pp) {

    # Get all aggregation unit ids, in descending order of sample count
    aggUnitData <- aggUnitData[order(-aggUnitData$SampleCount),]
    
    aggUnitGids <- rownames(aggUnitData)							#; print(aggUnitGids)
    freqBarData <- NULL
    for (aIdx in 1:length(aggUnitGids)) {
        aggUnitGid <- aggUnitGids[aIdx]								#; print(aggUnitGid) 
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
        
        # Get all the cluster for this unit, ordered by name
        unitClusterData <- clusterShareData[which(clusterShareData$UnitId==aggUnitGid),]	#; print(unitClusterData)
        unitClusters <- as.character(unitClusterData$Cluster)					#; print(unitClusters)
        unitClusterData <- unitClusterData[order(unitClusters, decreasing=TRUE),]		#; print(unitClusterData)
        unitRowCount <- nrow(unitClusterData)							#; print(unitRowCount)
        
        # For each cluster, work out the y boundaries
        if (unitRowCount==1) {
            y1 <- y0
            y2 <- y0 + (unitClusterData$ClusterCount[1] * sHeight)
	} else {
            y2 <- vector(mode="numeric", length=unitRowCount)
            y <- y0
            for (i in 1:unitRowCount) {
                y <- y + (unitClusterData$ClusterCount[i] * sHeight)
                y2[i] <- y
            }
            y1 <- c(y0, y2[1:(unitRowCount-1)])
            x1 <- rep(x1, unitRowCount)
            x2 <- rep(x2, unitRowCount)
	}											#; print(x1); print(x2); print(y1); print(y2)
        unitBarData <- data.frame(unitClusterData, x1=x1, x2=x2, y1=y1, y2=y2)			#; print(nrow(unitBarData))
        freqBarData <- rbind(freqBarData, unitBarData)
    }												#; print(freqBarData)
    
    # Silly trick to make the package checker happy... :-(
    Cluster <- NULL

    # Got the dataframe for drawing all the rectangles, now do the drawing and colouring
    mapPlot <- mapPlot +
               ggplot2::geom_rect(ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=Cluster), 
                                  data=freqBarData, inherit.aes=FALSE,
                                  colour="gray25", linewidth=pp$pieLineWidth, show.legend=FALSE) +
               ggplot2::scale_fill_manual(values=clusterPalette)
    mapPlot
}
