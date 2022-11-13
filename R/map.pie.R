###############################################################################
# Map Aggregated Measure Analysis
###############################################################################
#
pieMap.execute <- function(userCtx, datasetName, sampleSetName, interval, mapType, baseMapInfo, aggregation, measures, params) {

    # Get the context and trim it by time interval
    sampleSet <- userCtx$sampleSets[[sampleSetName]]
    plotCtx <- analysis.trimContextByTimeInterval (sampleSet$ctx, interval)
    if (is.null(plotCtx)) {
        print(paste("No samples found- skipping interval", interval$name))
        return()
    }											#; print(str(plotCtx))
    dataset <- plotCtx[[datasetName]]
    sampleMeta <- dataset$meta
    if (nrow(sampleMeta)==0) {
        print(paste("No samples found - skipping interval", interval$name))
        return()
    }
    config <- userCtx$config
    
    # Get the output folders
    dataFolder <- getOutFolder(config, sampleSetName, c(paste("map", mapType, sep="-"), "data"))
    
    # Silly trick to make the package checker happy... :-(
    lon <- lat <- label <- NULL
    #
    # Check the measures specified are valid.
    # If they are, ensure we have colour palettes for these measures, creating them if necessary.
    # For a given measure, all plots for this sample set use the same palette, otherwise the viewer will be confused when looking at multiple maps.
    #
    measures <- pieMap.checkAllelePropMeasures (userCtx, sampleSetName, measures)
    
    #
    # Now compute the aggregation units, the values to be plotted, and make the map
    #
    for (aggIdx in 1:length(aggregation)) {
        aggLevel <- as.integer(aggregation[aggIdx])        						#; print(aggLevel)
        aggLevelIdx <- aggLevel + 1

        # Get the aggregated data for the aggregation units
        aggUnitData <- map.getAggregationUnitData (plotCtx, datasetName, aggLevel, sampleSetName, mapType, params, dataFolder)	#; print(aggUnitData)

        for (mIdx in 1:length(measures)) {
            measure <- measures[mIdx]		           									#; print(measure)
            pieMapData <- pieMap.buildCountData (userCtx, datasetName, sampleSetName, aggLevel, aggUnitData, measure)		#; print(pieMapData)
	    
            # Select the aggregation units to be plotted
            # In this case, those that have allele count data for the measure being plotted
            selAggUnitIds <- unique(pieMapData$UnitId)
            selAggUnitData <- aggUnitData[which(aggUnitData$UnitId %in% selAggUnitIds),]		#; print(nrow(selAggUnitData))
            
            # Compute marker sizes. 
            pieSizes <- pieMap.getAggUnitPieSizes (selAggUnitData, params)				#; print(pieSizes); print(pieMapData$UnitId)
            pieMapData$PieSize <- as.integer(pieSizes[pieMapData$UnitId])				#; print(pieMapData)

            # Do the actual plot, starting with the background map
            mapPlot <- baseMapInfo$baseMap
            
            # If we need to show aggregation unit names, we need to compute the label positioning and plot before the markers
            showMarkerNames <- analysis.getParam ("map.markerNames", params)
            if (showMarkerNames) {
                lp <- map.computeLabelParams (selAggUnitData, baseMapInfo)
                mapPlot <- mapPlot + 
                    ggrepel::geom_label_repel(data=lp, ggplot2::aes(x=lon, y=lat, label=label), 
                                              size=4.5, fontface="bold", color="darkgray",
                                              hjust=lp$just, vjust=0.5, nudge_x=lp$x, nudge_y=lp$y, label.padding=grid::unit(0.2, "lines"))
            }
            
            # This function replaces aes_strng() allowing the use of column names with dashes
            fn_aesString <- get("aes_string", asNamespace("ggplot2"))
            aes_string2 <- function(...){
                args <- lapply(list(...), function(x) sprintf("`%s`", x))
                #do.call(aes_string, args)
                do.call(fn_aesString, args)
            }

            # Now add the markers, coloured according to the palette 
            valuePalette <- pieMap.getMeasurePalette (userCtx, sampleSetName, measure)

            # Now add the pie chart markers
            mapPlot <- mapPlot +
                    ggforce::geom_arc_bar(data=pieMapData, stat="pie", inherit.aes=FALSE, 
                                          ggplot2::aes(x0=Longitude, y0=Latitude, r0=0, r=(PieSize/100),
                                          fill=Allele, amount=AlleleCount),
                                          colour="gray25", show.legend=TRUE) +
                    ggplot2::scale_fill_manual(values=valuePalette)

            # Now add the decorative elements
if (FALSE) {
            if (mapType=="location") {
                valueLabels <- rep("", length(mValues))
            } else {
                valueLabels <- round(mValues, digits=2)
            }
            if (mapType=="sampleCount") {
                mapPlot <- mapPlot +
                           ggplot2::geom_text(data=selAggUnitData, ggplot2::aes_string(x="Longitude", y="Latitude", colour=colourAdmDivCol),
                                              label=valueLabels, hjust=0.5, vjust=0.5, size=4.5, fontface="bold", show.legend=FALSE) +
                           ggplot2::scale_colour_manual(values=admDivTextPalette)
            } else {
                mapPlot <- mapPlot +
                           ggplot2::geom_text(data=selAggUnitData, ggplot2::aes_string(x="Longitude", y="Latitude"), 
                                              label=valueLabels, hjust=0.5, vjust=0.5, size=4.5, fontface="bold")
            }
}            
            mapPlot <- mapPlot +
                       ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",size = ggplot2::rel(1.2), hjust = 0.5),
                       panel.background = ggplot2::element_rect(colour = NA),
                       plot.background = ggplot2::element_rect(colour = NA),
                       axis.title = ggplot2::element_text(face = "bold",size = ggplot2::rel(1)),
                       axis.title.y = ggplot2::element_text(angle=90,vjust =2),
                       axis.title.x = ggplot2::element_text(vjust = -0.2))
	    
            # Save to file. the size in inches is given in the config.
            mapSize  <- analysis.getParam ("plot.size", params)
            plotFolder <- getOutFolder(ctx$config, sampleSetName, c(paste("map", mapType, sep="-"), "plots"))
            aggLabel <- map.getAggregationLabels(aggLevel)
            if (mapType=="sampleCount") {
                aggLabel <- paste(aggLabel, datasetName, sep="-")
            }												#;print(aggLabel)
            graphicFilenameRoot  <- paste(plotFolder, paste("map", sampleSetName, aggLabel, measure, sep="-"), sep="/")
            if (!is.null(interval$name)) {
                graphicFilenameRoot  <- paste(graphicFilenameRoot, interval$name, sep="-")		#;print(graphicFilenameRoot)
            }
            ggplot2::ggsave(plot=mapPlot, filename=paste(graphicFilenameRoot,"png",sep="."), device="png", width=mapSize$width, height=mapSize$height, units="in", dpi=300)
        }
    }
}

#
###############################################################################
# Graphical rendition of site markers
################################################################################
#
pieMap.getAggUnitPieSizes <- function(aggUnitData, params) {
    # Compute pie chart sizes. 
    # If only one size was given in the config, then the markers will be constant size/
    # If there are two sizes, then merkers will be sized proportional to the number of samples, with the smaller
    # size representing 1 sample, and the larger size representing the numer of samples in the largest aggregation
    mSizeParam <- analysis.getParam ("map.markerSize", params)		#; print(mSizeParam)
    if (length(mSizeParam) > 1) {
        minSize  <- mSizeParam[1]
        maxSize  <- mSizeParam[2]
        sizeRange  <- maxSize - minSize
        counts <- aggUnitData$SampleCount
        minCount <- 1
        maxCount <- max(counts)
        countRange <- maxCount - minCount
        markerSizes <- minSize + (sizeRange * ((counts - minCount) / countRange))
    } else {
        markerSizes <- rep(mSizeParam, nrow(aggUnitData))
    }
    markerSizes <- as.integer(markerSizes)
    names(markerSizes) <- aggUnitData$UnitId
    markerSizes								#; print(markerSizes)
}
#
###############################################################################
# Estimation of marker measures (drug resistance and diversity)
################################################################################
#
# Check the allele proportion measures specified are valid.
# If they are, ensure we have colour palettes for these measures, creating them if necessary.
# For a given measure, all plots for this sample set use the same palette, otherwise the viewer will be confused when looking at multiple maps.
#
pieMap.checkAllelePropMeasures <- function(userCtx, sampleSetName, measures) {
    config <- userCtx$config 
    allMeasures <- c(config$cluster.stats.alleleCounts, config$amplificationColumns, config$drugResistancePositions)		#; print(allMeasures)
    if ("ALL" %in% measures) {
        measures <- allMeasures
    } else {
        #
        # Check the measure name is valid and ensure the relevant colour palette is loaded measure name is valid
        #
        sampleSet <- userCtx$sampleSets[[sampleSetName]]
        palettesList <- sampleSet$valuePalettes
        for (i in 1 : length(measures)) {
            measure <- measures[i]							#; print(measure)
            if (!(measure %in% allMeasures)) {
                stop(paste("Invalid allele for allele proportions:", measure))
            }
            p <- palettesList[[measure]]
            if (is.null(p)) {
                p <- pieMap.loadMeasurePalette (userCtx, sampleSetName, measure)	#; print(paste("Loading palette",measure))
                palettesList[[measure]] <- p						#; print(palettesList[[measure]])
            }
        }
    }
    measures
}
#
pieMap.getMeasurePalette  <- function(userCtx, sampleSetName, measure) {
    sampleSet <- userCtx$sampleSets[[sampleSetName]]		#; print(userCtx$sampleSets[[sampleSetName]]$valuePalettes[[measure]])
    p <- sampleSet$valuePalettes[[measure]]			#; print(p)
    p
}
#
pieMap.getMeasureAllValues <- function(userCtx, sampleSetName, measure) {
    sampleSet <- userCtx$sampleSets[[sampleSetName]]
    dataset <- sampleSet$ctx$unfiltered
    sampleMeta <- dataset$meta				#; print(nrow(sampleMeta))
    valueCounts <- meta.getColumnValueCounts (sampleMeta, measure, excludeMultiValues=TRUE)
    values <- names(valueCounts)
    values
}
#
pieMap.loadMeasurePalette <- function(userCtx, sampleSetName, measure) {
    values <- pieMap.getMeasureAllValues (userCtx, sampleSetName, measure)	#; print(values)
    #
    # Put "WT" as the first value for pie charts
    #
    if ("WT" %in% values) {
        wtIdx <- which(values == "WT")
        wtValue <- values[wtIdx]
        values <- c(wtValue, values[-wtIdx])
    }
    userPalette <- graphics.getColourPalette (userCtx)
    measurePalette <- rep_len(userPalette, length(values))	#; print(measurePalette)
    names(measurePalette) <- values				#; print(measurePalette)
    measurePalette
}
#
# For each aggregation unit, we get a count of each unique allele, ordered according to the palette order
#
pieMap.buildCountData <- function(userCtx, datasetName, sampleSetName, aggLevel, aggUnitData, measure) {

    sampleSet <- userCtx$sampleSets[[sampleSetName]]
    dataset <- sampleSet$ctx[[datasetName]]
    sampleMeta <- dataset$meta							#; print(str(sampleMeta))

    # Get all aggregation unit ids
    aggUnitGids <- rownames(aggUnitData)					#; print(aggUnitGids)
    
    # Create aggregation index for each sample (the id of the aggregation unit where the sample originates)
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta)		#; print(aggIndex)
    
    # Get the palette for this measure, which will give us both the colours and the order of the alleles
    valuePalette <- pieMap.getMeasurePalette (userCtx, sampleSetName, measure)	#; print(valuePalette)
    valueOrder <- 1:length(valuePalette)					#; print(valueOrder)
    names(valueOrder) <- names(valuePalette)
    
    # Get the data for all aggregation units
    countData <- NULL
    for (aIdx in 1:length(aggUnitGids)) {
    
        # Get the sample data to be aggregated for this unit
        aggUnitGid <- aggUnitGids[aIdx]						#; print(aggUnitGid)
        aggMeta <- sampleMeta[which(aggIndex == aggUnitGid),]			#; print(nrow(aggMeta))
        
        # Get the allele counts
        vCounts <- meta.getColumnValueCounts (aggMeta, measure, excludeMultiValues=TRUE)		#; print(vCounts)
        if (length(vCounts) == 0) {
            next
        }
	vNames <- names(vCounts)				#; print(vNames)
	vOrder <- as.integer(valueOrder[vNames])		#; print(vOrder)
	vColour <- valuePalette[vNames]				#; print(vColour)
	
	# Create a dataframe with the necessary data for this unit, copying the unit data in each row
        unitData <- aggUnitData[aIdx,]				#; print(unitData)
        df <- data.frame(matrix(unitData,ncol=ncol(aggUnitData),nrow=length(vCounts),byrow=TRUE))
        colnames(df) <- colnames(aggUnitData)			#; print(df)
        # Now add the allele counts, sort and merge with previous units
        df$Allele      <- vNames
        df$AlleleCount <- vCounts
        df$Order       <- vOrder
        df$Colour      <- vColour				#; print(df)
        
        df <- df[order(df$Order),]				#; print(df)
        if (is.null(countData)) {
            countData <- df
        } else {						#; print(str(df))
            countData <- rbind(countData, df)
        }
    }
    countData$UnitId <- as.character(countData$UnitId)
    countData$Allele <- factor(countData$Allele)
    countData$AlleleCount <- as.integer(countData$AlleleCount)
    countData$Order <- as.integer(countData$Order)
    countData$Latitude <- as.numeric(countData$Latitude)
    countData$Longitude <- as.numeric(countData$Longitude)
    countData$SampleCount <- as.numeric(countData$SampleCount)
    
    countData
}

