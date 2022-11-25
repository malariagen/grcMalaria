###############################################################################
# Map Aggregated Measure Analysis
###############################################################################
#
pieMap.execute <- function(userCtx, datasetName, sampleSetName, interval, mapType, baseMapInfo, aggregation, measures, params) {

    # Get the context and trim it by time interval
    sampleSet <- userCtx$sampleSets[[sampleSetName]]
    plotCtx <- analysis.trimContextByTimeInterval (sampleSet$ctx, interval)		#; print(interval$name)
    if (is.null(plotCtx)) {
        print(paste("No samples found- skipping interval", interval$name))
        return()
    }											#; print(str(plotCtx))
    dataset <- plotCtx[[datasetName]]
    sampleMeta <- dataset$meta								#; print(nrow(sampleMeta))
    if (nrow(sampleMeta)==0) {
        print(paste("No samples found - skipping interval", interval$name))
        return()
    }
    config <- userCtx$config								#;print(sampleMeta[,"Country"])
    
    # Get the output folders
    dataFolder <- getOutFolder(config, sampleSetName, c(paste("map", mapType, sep="-"), "data"))
    
    # Silly trick to make the package checker happy... :-(
    lon <- lat <- label <- NULL
    Longitude <- Latitude <- PieSize <- Allele <- AlleleCount <- NULL
    
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
            pieMapData <- pieMap.buildCountData (plotCtx, datasetName, sampleSetName, aggLevel, aggUnitData, measure)		#; print(pieMapData)
	    
            # Select the aggregation units to be plotted
            # In this case, those that have allele count data for the measure being plotted
            selAggUnitIds <- unique(pieMapData$UnitId)
            selAggUnitData <- aggUnitData[which(aggUnitData$UnitId %in% selAggUnitIds),]		#; print(nrow(selAggUnitData))
            
            # Compute pie chart sizes. 
	    pieSizes <- pieMap.getAggUnitPieSizes (selAggUnitData, params)		#; print(pieSizes); print(pieMapData$UnitId)
	    pieMapData$PieSize <- as.numeric(pieSizes[pieMapData$UnitId])		#; print(pieMapData)
            
            # Compute label coordinates for the pie chart segments
            pieMapData <- pieMap.computeLabelCoordinates (pieMapData, 1.2)		#; print(pieMapData)

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
            valuePalette <- pieMap.getMeasurePalette (plotCtx, sampleSetName, measure)

            # Now add the pie chart markers
            mapPlot <- mapPlot +
                ggforce::geom_arc_bar(data=pieMapData, stat="pie", inherit.aes=FALSE, 
                                          ggplot2::aes(x0=Longitude, y0=Latitude, r0=0, r=PieSize,
                                          fill=Allele, amount=AlleleCount),
                                          colour="gray25", show.legend=TRUE) +
                ggplot2::scale_fill_manual(values=valuePalette)
                    
            # Now add the pie slice labels
            showAlleleCounts <- FALSE			# TBD
            if (showAlleleCounts) {
                mapPlot <- mapPlot +
	            ggplot2::geom_text(data=pieMapData, ggplot2::aes(x=LabelX, y=LabelY, label=AlleleCount),
	                               show.legend=FALSE)
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
            plotFolder <- getOutFolder(userCtx$config, sampleSetName, c(paste("map", mapType, sep="-"), "plots"))
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
pieMap.computeLabelCoordinates <- function(pieMapData, labelRadius=0.5) {		#; print(pieMapData)
    resultDf <- NULL
    aggUnitIds <- as.character(unique(pieMapData$UnitId))
    for (uIdx in 1:length(aggUnitIds)) {
        unitId <- aggUnitIds[uIdx]
        unitData <- pieMapData[which(pieMapData$UnitId == unitId),]
        unitData <- unitData[order(unitData$Order),]				#; print(unitData)
        
        counts <- unitData$AlleleCount
        ctotal <- sum(counts)
        counts <- c(0, cumsum(counts))
        c1 <- counts[1:(length(counts)-1)]
        c2 <- counts[2:length(counts)]
        cmid <- (c1 + c2) / 2
        angle <- 2 * pi * cmid / ctotal
        lRadius <- unitData$PieSize * labelRadius

        unitData$LabelX <- unitData$Longitude + (sin(angle) * lRadius)
        unitData$LabelY <- unitData$Latitude  + (cos(angle) * lRadius)

        if (is.null(resultDf)) {
            resultDf <- unitData
        } else {
            resultDf <- rbind(resultDf, unitData)
        }
    }
    resultDf									#; print(resultDf)
}
            

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
    #markerSizes <- as.numeric(markerSizes)
    markerSizes <- as.numeric(markerSizes / 100)		# The /100 is an arbitrary scaling factor, will fix later 
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
pieMap.checkAllelePropMeasures <- function(ctx, sampleSetName, measures) {
    userCtx <- ctx$rootCtx
    config <- userCtx$config 
    allMeasures <- c(config$countColumns, config$amplificationColumns, config$drugResistancePositions)		#; print(allMeasures)
    if ("ALL" %in% measures) {
        measures <- allMeasures
    }
    #
    # Check the measure name is valid and ensure the relevant colour palette is loaded measure name is valid
    #
    sampleSet <- userCtx$sampleSets[[sampleSetName]]
    palettesList <- sampleSet$valuePalettes
    #for (i in 1 : 2) {
    for (i in 1 : length(measures)) {
        measure <- measures[i]							#; print(measure)
        if (!(measure %in% allMeasures)) {
            stop(paste("Invalid allele for allele proportions:", measure))
        }
        p <- palettesList[[measure]]						#; print(p)
        if (is.null(p)) {
            p <- pieMap.loadMeasurePalette (userCtx, sampleSetName, measure)	#; print(paste("Loading palette",measure))
            palettesList[[measure]] <- p					#; print(palettesList[[measure]])
        }
    }
    measures
}
#
pieMap.getMeasurePalette  <- function(ctx, sampleSetName, measure) {		#; print(names(ctx$rootCtx$sampleSets[[sampleSetName]]))
    sampleSet <- ctx$rootCtx$sampleSets[[sampleSetName]]	#; print(ctx$rootCtx$sampleSets[[sampleSetName]]$valuePalettes[[measure]])
    p <- sampleSet$valuePalettes[[measure]]			#; print(p)
    p
}
#
pieMap.getMeasureAttributes <- function(config, measure) {
    if (pieMap.isAminoPosition(config, measure)) {
        posData <- meta.getPositionData (config, measure)
        columnName <- posData$columnName
        wtValue <- posData$ref
    } else {
        columnName <- measure
        wtValue <- "WT"
    }
    list(wtValue=wtValue, columnName=columnName)
}
#
pieMap.isAminoPosition <- function (config, measure) {
    measure %in% config$drugResistancePositions
}
#
pieMap.getMeasureAllValues <- function(ctx, sampleSetName, measure, excludeMultiValues=TRUE, wtFirst=TRUE) {
    userCtx <- ctx$rootCtx
    sampleSet <- userCtx$sampleSets[[sampleSetName]]
    dataset <- sampleSet$ctx$unfiltered
    sampleMeta <- dataset$meta					#; print(nrow(sampleMeta))
    #	
    attr <- pieMap.getMeasureAttributes (userCtx$config, measure)	#; print(attr)
    columnName <- attr$columnName
    wtValue <- attr$wtValue
    #
    valueCounts <- meta.getColumnValueCounts (sampleMeta, columnName, excludeMultiValues)
    values <- names(valueCounts)
    #
    # Put WT allele as the first value for pie charts
    #
    if (wtFirst) {
        if (wtValue %in% values) {
            wtIdx <- which(values == wtValue)
            wtValue <- values[wtIdx]
            values <- c(wtValue, values[-wtIdx])
        }
    }
    values
}
#
pieMap.loadMeasurePalette <- function(ctx, sampleSetName, measure) {
    userCtx <- ctx$rootCtx
    values <- pieMap.getMeasureAllValues (userCtx, sampleSetName, measure)	#; print(values)
    userPalette <- graphics.getColourPalette (userCtx)
    measurePalette <- rep_len(userPalette, length(values))	#; print(measurePalette)
    names(measurePalette) <- values				#; print(measurePalette)
    measurePalette
}
#
# For each aggregation unit, we get a count of each unique allele, ordered according to the palette order
#
pieMap.buildCountData <- function(ctx, datasetName, sampleSetName, aggLevel, aggUnitData, measure) {

    dataset <- ctx[[datasetName]]
    sampleMeta <- dataset$meta							#; print(str(sampleMeta))

    # Get all aggregation unit ids
    aggUnitGids <- rownames(aggUnitData)					#; print(aggUnitGids)
    
    # Create aggregation index for each sample (the id of the aggregation unit where the sample originates)
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta)		#; print(aggIndex)
    
    # Get the palette for this measure, which will give us both the colours and the order of the alleles
    valuePalette <- pieMap.getMeasurePalette (ctx, sampleSetName, measure)	#; print(valuePalette)
    valueOrder <- 1:length(valuePalette)					#; print(valueOrder)
    names(valueOrder) <- names(valuePalette)
    
    # Get the data for all aggregation units
    countData <- NULL
    for (aIdx in 1:length(aggUnitGids)) {
    
        # Get the sample data to be aggregated for this unit
        aggUnitGid <- aggUnitGids[aIdx]						#; print(aggUnitGid)
        aggMeta <- sampleMeta[which(aggIndex == aggUnitGid),]			#; print(nrow(aggMeta))
        
        # Get the allele counts
        attr <- pieMap.getMeasureAttributes (ctx$config, measure)					#; print(attr)
        vCounts <- meta.getColumnValueCounts (aggMeta, attr$columnName, excludeMultiValues=TRUE)	#; print(vCounts)
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

