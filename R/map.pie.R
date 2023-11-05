###############################################################################
# Map Aggregated Measure Analysis
###############################################################################
#
pieMap.executeMap <- function(map) {
    mapMaster   <- map$master
    mapType     <- mapMaster$type
    measure     <- map$measure
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
    ctx        <- map$mapCtx							#; print(str(ctx))
    dataset    <- ctx[[datasetName]]
    sampleMeta <- dataset$meta
    if (nrow(sampleMeta)==0) {
        return(NULL)
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
    lon <- lat <- label <- Longitude <- Latitude <- PieSize <- Allele <- AlleleCount <- LabelX <- LabelY <- NULL
    #
    # Now compute the aggregation units, the values to be plotted, and make the map
    # Get the aggregated data for the aggregation units
    #
    aggLevel <- as.integer(map$aggregation)     				#; print(aggLevel)   
    aggUnitData <- map$aggUnitData						#; print(aggUnitData)
    #
    # Get a table of aggregation units and the sample count for each allele
    #
    pieMapData <- pieMap.buildCountData (ctx, datasetName, sampleSet$name, aggLevel, aggUnitData, measure, params)	#; print(pieMapData)
    selAggUnitIds <- pieMapData$UnitId
    selAggUnitData <- aggUnitData[which(aggUnitData$UnitId %in% selAggUnitIds),]				#; print(selAggUnitData)
    if (nrow(selAggUnitData) == 0) {
        return(NULL)
    }
    #
    # Compute pie chart sizes.
    #
    pieSizes <- pieMap.getAggUnitPieSizes (pieMapData, params)			#; print(pieSizes)   #; print(pieMapData$UnitId)
    pieMapData$PieSize <- as.numeric(pieSizes[pieMapData$UnitId])		#; print(pieMapData)
    #
    # Do the actual plot, starting with the background map
    #
    mapPlot <- baseMapInfo$baseMap
    #
    # If we need to show aggregation unit names, we need to compute the label positioning and plot
    #
    showMarkerNames <- param.getParam ("map.markerNames", params)
    if (showMarkerNames) {
        objectSizes <- (pieSizes * 80) / params$plot.scaleFactor		#; print(objectSizes)
        mapPlot <- map.addAggregationUnitNameLayer (mapPlot, selAggUnitData, baseMapInfo, params, markerSize=objectSizes)
    }
    #
    # Get the categorical palette for this measure, which will give us both the colours and the order of the alleles
    #
    catPalette <- pieMap.getCategoryPalette (ctx, sampleSet$name, measure, params)	#; print(catPalette)
    alleleNames <- names(catPalette)
    #
    # Plot the pies
    #
    # Nasty trick we have to do because there is no linewidth aesthetic at present
    ggplot2::update_geom_defaults(ggforce::GeomArcBar, ggplot2::aes(linewidth=!!params$pieLineWidth))
    mapPlot <- mapPlot +
        scatterpie::geom_scatterpie (ggplot2::aes(x=Longitude, y=Latitude, r=PieSize),
                                     data=pieMapData, cols=names(catPalette), sorted_by_radius=TRUE) +
        ggplot2::scale_fill_manual(values=catPalette, breaks=alleleNames)
    mapPlot
}

###############################################################################
# Allele Counts for Aggregation Units
###############################################################################
#
# For each aggregation unit, we get a count of each unique allele, ordered according to the palette order
#
pieMap.buildCountData <- function(ctx, datasetName, sampleSetName, aggLevel, aggUnitData, measure, params) {
    dataset <- ctx[[datasetName]]
    sampleMeta <- dataset$meta							#; print(str(sampleMeta))
    #
    # Get all aggregation unit ids
    #
    aggUnitGids <- rownames(aggUnitData)					#; print(aggUnitGids)
    aggUnitCount <- length(aggUnitGids)						#; print(aggUnitCount)
    #
    # Create aggregation index for each sample (the id of the aggregation unit where the sample originates)
    #
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta)		#; print(aggIndex)
    #
    # Get the categorical palette for this measure, which will give us both the colours and the order of the alleles
    #
    catPalette <- pieMap.getCategoryPalette (ctx, sampleSetName, measure, params)	#; print(catPalette)
    catNames <- names(catPalette)						#; print(catNames)
    catCount <- length(catPalette)						#; print(catCount)
    catOrder <- 1:catCount
    names(catOrder) <- catNames							#; print(catOrder)
    #
    #  Create a matrix for the allele counts 
    #
    m <- matrix(0, nrow=aggUnitCount, ncol=catCount, dimnames=list(aggUnitGids, catNames))
    otherCounts <- vector(mode="integer", length=aggUnitCount)
    totals      <- vector(mode="integer", length=aggUnitCount)
    #
    # Get the allele counts and other data for each aggregation units
    #
    for (aIdx in 1:aggUnitCount) {
        #
        # Get the sample data to be aggregated for this unit
        #
        aggUnitGid <- aggUnitGids[aIdx]						#; print(aggUnitGid)
        aggMeta <- sampleMeta[which(aggIndex == aggUnitGid),]			#; print(nrow(aggMeta))
        #
        # Get the allele counts and add them to the matrix
        #
        attr <- pieMap.getMeasureAttributes (ctx$config, measure)					#; print(attr)
        vCounts <- meta.getColumnValueCounts (aggMeta, attr$columnName, excludeMultiValues=TRUE)	#; print(vCounts)
        if (length(vCounts) == 0) {
            next
        }
        vNames <- names(vCounts)						#; print(vNames)
        total <- 0
        for (i in 1:length(vCounts)) {						#; print(pieMapData[i,])
            vName <- vNames[i]; vCount <- vCounts[i]
            if (vName %in% catNames) {
                m[aIdx,vName] <- vCount
            } else {
                otherCounts[aIdx] <- otherCounts[aIdx] + vCount
            }
            total <- total + vCount
        }
        totals[aIdx] <- total
    }
    #
    # Add the necessary mapping data to the allele counts
    #
    countData <- as.data.frame(m)
    if ("Other" %in% catNames) {
        countData$Other <- otherCounts
    } else {
        countData <- cbind(countData, Other=otherCounts)
    }
    countData <- cbind(countData, Latitude=aggUnitData$Latitude, Longitude=aggUnitData$Longitude, 
                       SampleCount=aggUnitData$SampleCount, UnitId=aggUnitData$UnitId, AdmDivName=aggUnitData$AdmDivName)
    rownames(countData) <- countData$UnitId
    countData
}

###############################################################################
# Sizes of pie plots
###############################################################################
#
pieMap.getAggUnitPieSizes <- function(pieMapData, params) {			#; print(pieMapData)
    if (nrow(pieMapData) == 0) {
        return (numeric(0))
    }
    # Compute pie chart sizes. 
    # If only one size was given in the config, then the markers will be constant size/
    # If there are two sizes, then merkers will be sized proportional to the number of samples, with the smaller
    # size representing 1 sample, and the larger size representing the numer of samples in the largest aggregation
    mSizeParam <- param.getParam ("map.markerSize", params)			#; print(mSizeParam)
    if (length(mSizeParam) > 1) {
        minSize  <- mSizeParam[1]
        maxSize  <- mSizeParam[2]
        sizeRange  <- maxSize - minSize
        counts <- pieMapData$SampleCount					#; print(counts)
        minCount <- 1
        maxCount <- max(counts)
        countRange <- maxCount - minCount
        markerSizes <- minSize + (sizeRange * ((counts - minCount) / countRange))
    } else {
        markerSizes <- rep(mSizeParam, nrow(pieMapData))
    }
    markerSizes <- as.numeric(markerSizes) / 80		# The /80 is an arbitrary scaling factor, determined by visual comparison with other plots
    names(markerSizes) <- pieMapData$UnitId
    markerSizes									#; print(markerSizes)
}

###############################################################################
# Categorical Colour Palettes
###############################################################################
#
# Get the palette to be used for this measure in this sampleset
# (the same palette may be used for multiple time slices, so we keep it in the context for reuse)
#
pieMap.getCategoryPalette <- function (ctx, sampleSetName, measure, params) {
    userCtx <- ctx$rootCtx
    sampleSet <- userCtx$sampleSets[[sampleSetName]]
    p <- sampleSet$categoryPalettes[[measure]]
    if (is.null(p)) {
        p <- params$map.alleleColours
        if (is.null(p)) {
            values <- pieMap.getMeasureAllValues (userCtx, sampleSetName, measure)	#; print(values)
            userPalette <- graphics.getColourPalette (userCtx)
            p <- rep_len(userPalette, length(values))					#; print(p)
            names(p) <- values
        }
        if (!("Other" %in% names(p))) {
            p <- c(p, Other="white")
        }
        sampleSet$categoryPalettes[[measure]] <- p					#; print(p)
    }
    p
}

###############################################################################
# Estimation of marker measures
################################################################################
#
# Check the allele proportion measures specified are valid.
# If they are, ensure we have colour palettes for these measures, creating them if necessary.
# For a given measure, all plots for this sample set use the same palette, otherwise the viewer will be confused when looking at multiple maps.
#
pieMap.resolveMeasures <- function(ctx, sampleSetName, params) {
    userCtx <- ctx$rootCtx
    config <- userCtx$config 
    measures <- param.getParam ("analysis.measures", params)
    #
    #
    #
    allMeasures <- c(config$countColumns, config$amplificationColumns, config$drugResistancePositions)		#; print(allMeasures)
    if ("ALL" %in% measures) {
        measures <- allMeasures
    }
    #
    # Check the measure name is valid
    #
    #for (i in 1 : 2) {
    for (i in 1 : length(measures)) {
        measure <- measures[i]							#; print(measure)
        if (!(measure %in% allMeasures)) {
            stop(paste("Invalid allele for allele proportions:", measure))
        }
    }
    measures
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
    sampleMeta <- dataset$meta						#; print(nrow(sampleMeta))
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
