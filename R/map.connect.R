###############################################################################
# Measures of connectedness
################################################################################
#
connectMap.getConnectednessMeasures <- function() {
    allMeasures <- c("similarity","meanDistance")
    allMeasures
}

###############################################################################
# Map Aggregated Measure Analysis
################################################################################
#
connectMap.execute <- function(ctx, analysisName, mapType, visType, aggregation, measures, params) {
    # Get the output folders
    dataOutFolder <- getOutFolder(ctx, analysisName, c(paste("map", mapType, sep="-"), "data"))
    # Get the sample metadata
    sampleMeta <- ctx$meta
    # Build the map necessary to display these samples
    # Construct a base plot for completing subsequent maps
    baseMapInfo <- map.buildBaseMap (ctx, analysisName, sampleMeta, dataOutFolder, params)
    
    # If "similarity" is one of the measures, remove it and add one measure for each similarity threshold (e.g. "similarity-ge0.80")
    measures <- connectMap.expandMeasures(measures, params)
    
    # Now compute the aggregation units, the values to be plotted, and make the map
    for (aggIdx in 1:length(aggregation)) {
        aggLevel <- as.integer(aggregation[aggIdx])        					#; print(aggLevel)
        aggLevelIdx <- aggLevel + 1

        # Get the aggregated data for the aggregation units
        aggUnitData <- map.getAggregationUnitData (aggLevel, ctx, analysisName, mapType, params, dataOutFolder)	#; print(aggUnitData)
        aggUnitPairData <- connectMap.estimateMeasures (aggLevel, aggUnitData, ctx, analysisName, mapType, measures, params, dataOutFolder)	#; print(aggUnitPairData)

        for (mIdx in 1:length(measures)) {
            measure <- measures[mIdx]						#; print(measure)
            if (connectMap.filterThreshold (measure)) {
                minValues <- as.numeric(analysis.getParam(paste("map.connect", measure, "min", sep="."), params, 0))
            } else {
                minValues <- 0
            }									#; print(minValues)
            for (mvIdx in 1:length(minValues)) {
                minValue <- minValues[mvIdx]
            
                # Select the aggregation unit pairs to be plotted
                # In this case, those thatwith value above the threshold
                mValues <- aggUnitPairData[,measure]
                selAggUnitPairData <- aggUnitPairData[which(mValues > minValue),]
                
                # Sort the pairs so that higher values get plotted last
                mValues <- selAggUnitPairData[,measure]
                selAggUnitPairData <- selAggUnitPairData[order(mValues),]
                
                # Do the actual plot, starting with the background map
                mapPlot <- baseMapInfo$baseMap
                
                # This function replaces aes_strng() allowing the use of column names with dashes
                fn_aesString <- get("aes_string", asNamespace("ggplot2"))
                aes_string2 <- function(...){
                    args <- lapply(list(...), function(x) sprintf("`%s`", x))
                    #do.call(aes_string, args)
                    do.call(fn_aesString, args)
                }
                
                # Now plot the connections
                mapPlot <- mapPlot +
                    ggplot2::geom_curve(aes_string2(x="Lon1", y="Lat1", xend="Lon2", yend="Lat2", size=measure, colour=measure),
                                        data=selAggUnitPairData, curvature=0.2, alpha=0.75) +
                    ggplot2::scale_size_continuous(guide=FALSE, range=c(0.25, 4)) +          					# scale for edge widths
                    ggplot2::scale_colour_gradientn(colours=c("skyblue1","midnightblue"))
                
                # Now add the markers
                mapPlot <- mapPlot +
                    ggplot2::geom_point(aes_string2(x="Longitude", y="Latitude"), data=aggUnitData, size=4, shape=19, col="red")
    	    
                # If we need to show aggregation unit names, we need to compute the label positioning and plot
                showMarkerNames <- analysis.getParam ("map.markerNames",   params, default.map.markerNames)
                if (showMarkerNames) {
                aggColName <- c("Country","AdmDiv1","AdmDiv2")[aggLevelIdx]
                    lp <- map.computeLabelParams (aggUnitData, aggColName, baseMapInfo)
                    mapPlot <- mapPlot +
                        ggplot2::geom_label_repel(ggplot2::aes(x=lon, y=lat, label=label), data=lp, size=4.5, fontface="bold", color="darkgray",
                                                  hjust=lp$just, vjust=0.5, nudge_x=lp$x, nudge_y=lp$y, label.padding=grid::unit(0.2, "lines"))
                }
                	    
                # Now add the decorative elements
                mapPlot <- mapPlot +
                           ggplot2::theme(
                               plot.title = element_text(face = "bold",size = ggplot2::rel(1.2), hjust = 0.5),
                               panel.background = ggplot2::element_rect(colour = NA),
                               plot.background = ggplot2::element_rect(colour = NA),
                               axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
                               axis.title.y = ggplot2::element_text(angle = 90,vjust = 2),
                               axis.title.x = ggplot2::element_text(vjust = -0.2))
    	    
    	        # Save to file. the size in inches is given in the config.
    	        mapSize  <- analysis.getParam ("map.size", params, default.map.size)
    	        plotFolder <- getOutFolder(analysisName, c(paste("map", mapType, sep="-"), "plots"))
                graphicFilenameRoot  <- paste(plotFolder, paste("map", analysisName, aggLevel, measure, sep="-"), sep="/")
                if (connectMap.filterThreshold (measure)) {
    	            graphicFilenameRoot  <- paste(graphicFilenameRoot, paste("ge", format(minValue, digits=2, nsmall=2), sep=""), sep="-")
                }
                ggplot2::ggsave(plot=mapPlot, filename=paste(graphicFilenameRoot,"png",sep="."), device="png", 
                                width=mapSize[1], height=mapSize[2], units="in", dpi=300)
            }
        }
    }
}

connectMap.filterThreshold  <- function(measure) {
    measure %in% c("meanDistance")
}

connectMap.expandMeasures <- function(measures, params) {
    result <- c()
    for (mIdx in 1:length(measures)) {
        measure <- measures[mIdx]
        if (measure == "similarity") {
            levels <- as.numeric(analysis.getParam("map.connect.similarity.min", params, 1.0))
            for (lIdx in 1 : length(levels)) {
                measureStr <- paste("similarity-ge", format(levels[lIdx], digits=2, nsmall=2), sep="")
                result <- c(result, measureStr)
            }
        } else {
            result <- c(result, measure)
        }
    }
    result
}

connectMap.getMeasureLevel <- function(measure, prefix) {
    if (!startsWith(measure, prefix)) {
        return (NA)				#; print("Incorrect prefix")
    }
    level <- substring(measure,nchar(prefix)+1)	#; print(level)
    as.numeric(level)
}

connectMap.estimateMeasures <- function (aggLevel, aggUnitData, ctx, analysisName, mapType, measures, params, dataFolder)	{

    sampleMeta   <- ctx$meta
    barcodeData  <- ctx$barcodes
    distData     <- ctx$distance

    # Create aggregation index for each sample (the id of the aggregation unit where the sample originates)
    aggUnitIds <- as.character(map.getAggregationUnitIds (aggLevel, sampleMeta, params))
    sampleIds <- rownames(sampleMeta)
    
    # Get all aggregation units
    aggUnits <- rownames(aggUnitData)							#; print(aggUnits) ; print(head(aggUnitData))
    
    # Get the data for all aggregation units
    measureData <- NULL
    for (a1Idx in 1:(length(aggUnits)-1)) {
        for (a2Idx in (a1Idx+1):length(aggUnits)) {
            # Get the sample data
            aggUnit1 <- aggUnits[a1Idx]							#; print(aggUnit1)
            lat1 <- aggUnitData$Latitude[a1Idx]
            lon1 <- aggUnitData$Longitude[a1Idx]
            aggSamples1 <- sampleIds[which(aggUnitIds == aggUnit1)]			#; print(length(aggSamples1))
            #aggBarcodes1 <- barcodeData[aggSamples1,]
            #
            aggUnit2 <- aggUnits[a2Idx]							#; print(aggUnit2)
            lat2 <- aggUnitData$Latitude[a2Idx]
            lon2 <- aggUnitData$Longitude[a2Idx]
            aggSamples2 <- sampleIds[which(aggUnitIds == aggUnit2)]			#; print(length(aggSamples2))
            #aggBarcodes2 <- barcodeData[aggSamples2,]
            #
            pairDist <- distData[aggSamples1,aggSamples2]				#; print(dim(pairDist))
        
            # Get the admin division values from the first sample of this unit (assuming the values are the same for all)
            mValues <- connectMap.estimateDistMeasures (pairDist, measures)			#; print (mValues)
            cValues <- c(aggUnit1, aggUnit2, lat1, lon1, lat2, lon2, mValues)
            measureData <- rbind(measureData, cValues)
        }
    }
    aggUnitPairData <- data.frame(measureData)						#; print (dim(aggUnitPairData))
    nc <- ncol(aggUnitPairData)
    aggUnitPairData[,3:nc] <- sapply(aggUnitPairData[,3:nc], as.numeric)		#; print(ncol(aggUnitPairData))
    colnames(aggUnitPairData) <- c("Unit1", "Unit2", "Lat1", "Lon1", "Lat2", "Lon2", measures)
    
    # Write out the aggregation unit data to file
    aggDataFilename  <- paste(dataFolder, "/AggregatedPairData-", analysisName, "-", aggLevel, ".tab", sep="")
    write.table(aggUnitPairData, file=aggDataFilename, sep="\t", quote=FALSE, row.names=FALSE)

    aggUnitPairData
}

connectMap.estimateDistMeasures <- function (distData, measures) {
    result <- c()
    dist <- unlist(distData)
    for (mIdx in 1:length(measures)) {
        measure <- measures[mIdx]
        if (startsWith(measure, "similarity")) {
            level <- connectMap.getMeasureLevel (measure, "similarity-ge")	#; print(level)
            maxDist = 1.0 - level						#; print(maxDist)
            value <- length(which(dist <= maxDist)) / length(dist)
        } else if (measure == "meanDistance") {
            value <- 1 - mean(dist) 
        } else {
            stop(paste("Invalid distance connectedness measure:", measure))
        }
        result <- c(result, value)
    }
    result
}

