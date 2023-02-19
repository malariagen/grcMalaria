###############################################################################
# Map Aggregated Measure Analysis
###############################################################################
#
markerMap.getDiversityMeasures <- function() {
    c("maxHaploFreq",
      "haploHet",
      "meanSnpHet",
      "medianDistance")
}
#
markerMap.executeMap <- function(map) {
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
    lon <- lat <- label <- NULL
    #
    # Now compute the aggregation units, the values to be plotted, and make the map
    # Get the aggregated data for the aggregation units
    #
    aggLevel <- as.integer(map$aggregation)     				#; print(aggLevel)   
    aggUnitData <- map$aggUnitData						#; print(aggUnitData)
    #
    # For sample count markers, the colour may be based on a different admin division level from the aggregation
    #
    if (mapType %in% c("sampleCount","location")) {
        colourAdmDivLevel <- param.getParam ("map.markerColourAggLevel", params)	#; print(colourAdmDivLevel) # should be 0 or 1
        colourAdmDivTitle <- ADM_DIV_LABELS[colourAdmDivLevel+1]			#; print(colourAdmDivCol)
        colourAdmDivCol   <- GID_COLUMNS[colourAdmDivLevel+1]				#; print(colourAdmDivCol)
        colourAdmDivs <- aggUnitData[,colourAdmDivCol]					#; print(colourAdmDivs)
        colourGids    <- unique(colourAdmDivs)						#; print(colourGids)
        colPalette    <- graphics.getColourPalette (userCtx)
        admDivPalette <- rep_len(colPalette, length.out=length(colourGids))
        admDivTextPalette <- graphics.makeTextPalette (admDivPalette)		
        names(admDivPalette) <- names(admDivTextPalette) <- colourGids			#; print(admDivPalette); print(admDivTextPalette)
        admDivPaletteLabels  <- map.getAdmDivNames (colourGids)				#; print(admDivPaletteLabels)
    }
    #
    # Select the aggregation units to be plotted (all, in the case of location markers)
    # In this case, those that have an NA for the measure being plotted (should be only for drug resistance)
    #
    selAggUnitData <- aggUnitData
    if (mapType != "location") {
        vals <- aggUnitData[,measure]					#; print(vals)
        selAggUnitData <- aggUnitData[which(!is.na(vals)),]		#; print(nrow(selAggUnitData))
    }
    #
    # Compute marker sizes. 
    #
    pointSizes <- markerMap.getAggUnitMarkerSizes (selAggUnitData, params)
    #
    # Do the actual plot, starting with the background map
    #
    mapPlot <- baseMapInfo$baseMap
    #
    # If we need to show aggregation unit names, we need to compute the label positioning 
    # and plot before the markers
    #
    showMarkerNames <- param.getParam ("map.markerNames", params)
    if (showMarkerNames) {
        lp <- map.computeLabelParams (selAggUnitData, baseMapInfo)
        mapPlot <- mapPlot + 
            ggrepel::geom_label_repel(data=lp, ggplot2::aes(x=lon, y=lat, label=label), 
                                      size=4.5, fontface="bold", color="darkgray",
                                      hjust=lp$just, vjust=0.5, nudge_x=lp$x, nudge_y=lp$y, 
                                      label.padding=grid::unit(0.2, "lines"))
    }
    #
    # Get the values to be displyed in the markers (except for the location markers)
    #
    if (mapType == "location") {
        valueLabels <- rep("", nrow(selAggUnitData))
    } else {
        mValues <- selAggUnitData[,measure]
        valueLabels <- round(mValues, digits=2)
    }
    #
    # Now add the markers, coloured according to the appropriate scale, depending on the type of map
    #
    if (mapType=="diversity") {
        markerColours <- param.getParam ("map.diversity.markerColours", params)
        # Two marker colours can be specified to create a gradient. If a single marker colour is 
        # specified, then create a gradient from white to that colour.
        if (length(markerColours) == 1) {
                    markerColours <- c("white", markerColours)
        }
        mMax <- max(mValues); scaleMax <- round(mMax,digits=1); scaleMax <- ifelse(scaleMax<mMax, scaleMax+0.1, scaleMax)
        mMin <- min(mValues); scaleMin <- round(mMin,digits=1); scaleMin <- ifelse(scaleMin>mMin, scaleMin-0.1, scaleMin)  #; print(paste(scaleMin, scaleMax))
        mapPlot <- mapPlot +
        ggplot2::geom_point(data=selAggUnitData, aes_string2(x="Longitude", y="Latitude", fill=measure), size=pointSizes, shape=21, stroke=2) +
        ggplot2::scale_fill_gradientn(limits=c(scaleMin,scaleMax), colours=markerColours, values=c(0,1))
    } else if (mapType=="sampleCount") {
        mapPlot <- mapPlot +
            ggplot2::geom_point(data=selAggUnitData, aes_string2(x="Longitude", y="Latitude", fill=colourAdmDivCol),
	                        size=pointSizes, shape=21, stroke=2) +
            ggplot2::scale_fill_manual(values=admDivPalette, labels=admDivPaletteLabels, name=colourAdmDivTitle,
                                       guide=ggplot2::guide_legend(override.aes=list(size=3,stroke=0.5)))
    }  else if (mapType=="location") {
        mapPlot <- mapPlot +
            ggplot2::geom_point(data=selAggUnitData, aes_string2(x="Longitude", y="Latitude", fill=colourAdmDivCol),
	                        size=pointSizes, shape=21, stroke=2)
    } else if (mapType=="drug") {
        mapPlot <- mapPlot +
            ggplot2::geom_point(aes_string2(x="Longitude", y="Latitude", fill=measure), 
	                        data=selAggUnitData, size=pointSizes, shape=21, stroke=2) +
            ggplot2::scale_fill_gradientn(limits=c(0,1), colours=c("green3","orange2","red3","red3"), values=c(0, 0.2, 0.75, 1))
    } else if (mapType=="mutation") {
        markerColours <- param.getParam ("map.prevalence.markerColours", params)
        #
        # Two marker colours can be specified to create a gradient. If a single marker colour is 
        # specified, then create a gradient from white to that colour.
        #
        if (length(markerColours) == 1) {
            markerColours <- c("white", markerColours)
        }
        scaleMin <- 0; scaleMax <- 1
        mapPlot <- mapPlot +
            ggplot2::geom_point(aes_string2(x="Longitude", y="Latitude", fill=measure), 
	                        data=selAggUnitData, size=pointSizes, shape=21, stroke=2) +
            ggplot2::scale_fill_gradientn(limits=c(scaleMin,scaleMax), colours=markerColours, values=c(0,1))
    }
    #	    
    # Now add the decorative elements
    #
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
    mapPlot <- mapPlot +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",size = ggplot2::rel(1.2),hjust = 0.5),
                           panel.background = ggplot2::element_rect(colour = NA),
                           plot.background = ggplot2::element_rect(colour = NA),
                           axis.title = ggplot2::element_text(face = "bold",size = ggplot2::rel(1)),
                           axis.title.y = ggplot2::element_text(angle=90,vjust =2),
                           axis.title.x = ggplot2::element_text(vjust = -0.2))
    #
    # Save to file. the size in inches is given in the config
    #
    mapSize  <- param.getParam ("plot.size", params)
    ggplot2::ggsave(plot=mapPlot, filename=map$plotFile, device="png", 
                    width=mapSize$width, height=mapSize$height, units="in", dpi=300)
}
#
###############################################################################
# Graphical rendition of site markers
################################################################################
#
markerMap.getAggUnitMarkerSizes <- function(aggUnitData, params) {
    # Compute marker sizes. 
    # If only one size was given in the config, then the markers will be constant size/
    # If there are two sizes, then merkers will be sized proportional to the number of samples, with the smaller
    # size representing 1 sample, and the larger size representing the numer of samples in the largest aggregation
    mSizeParam <- param.getParam ("map.markerSize", params)
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
        markerSizes <- mSizeParam
    }
    markerSizes
}
#
###############################################################################
# Estimation of marker measures (drug resistance and diversity)
################################################################################
#
markerMap.resolveMeasures <- function(mapType, measures, config) {
    # Get the admin division values from the first sample of this unit (assuming the values are the same for all)
    if (mapType=="diversity") {
        if ("ALL" %in% measures) {
            measures <- markerMap.getDiversityMeasures()
        }
    } else if (mapType=="drug") {
        if ("ALL" %in% measures) {
            measures <- config$drugs
        }
    } else if (mapType=="mutation") {
        if ("ALL" %in% measures) {
            measures <- config$drugResistanceMutations
        }
    }
    measures
}
#
markerMap.estimateMeasures <- function(ctx, datasetName, aggLevel, aggUnitData, sampleSetName, mapType, measures, params, dataFolder) {	#;print(measures)

    dataset <- ctx[[datasetName]]
    sampleMeta   <- dataset$meta
    barcodeData  <- dataset$barcodes
    distData     <- dataset$distance

    # Create aggregation index for each sample (the id of the aggregation unit where the sample originates)
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta)
    
    # Get all aggregation units
    aggUnitGids <- rownames(aggUnitData)				#; print(aggUnitGids)
    
    # Get the data for all aggregation units
    measureData <- matrix(nrow=0, ncol=length(measures), dimnames=list(c(),measures))
    for (aIdx in 1:length(aggUnitGids)) {

        # Get the sample data to be aggregated for this unit
        aggUnitGid <- aggUnitGids[aIdx]					#; print(aggUnitGid)
        aggSamplesMeta <- sampleMeta[which(aggIndex == aggUnitGid),]	#; print(nrow(aggSamplesMeta))
        aggSamples <- rownames(aggSamplesMeta)
        
        # Get the admin division values from the first sample of this unit (assuming the values are the same for all)
        if (mapType=="diversity") {
            aggBarcodes <- barcodeData[aggSamples,]
            aggDist <- distData[aggSamples,aggSamples]
            cValues <- markerMap.estimateDiversityMeasures (ctx, aggBarcodes, aggDist, measures)
        } else if (mapType=="drug") {
            cValues <- meta.getResistancePrevalence (ctx, aggSamplesMeta, measures, params)
        } else if (mapType=="mutation") {
            cValues <- meta.getMutationPrevalence (ctx, aggSamplesMeta, measures, params)
        }                                                           	#; print(cValues)
        measureData <- rbind(measureData, cValues)                 	#; print(measureData)
    }
    measureData <- as.data.frame(measureData)				#;  print(measureData)
    aggUnitData <- cbind(aggUnitData, measureData)			#;  print(dim(aggUnitData))
    for (mIdx in 1:length(measures)) {
         measure <- measures[mIdx]
         aggUnitData[,measure] <- as.numeric(aggUnitData[,measure])
    }

    # Write out the aggregation unit data to file
    aggDataFilename  <- paste(dataFolder, "/AggregatedData-", sampleSetName, "-", aggLevel, ".tab", sep="")
    utils::write.table(aggUnitData, file=aggDataFilename, sep="\t", quote=FALSE, row.names=FALSE)

    aggUnitData
}
#
# Diversity measure estimates from genetic barcodes.
# Note that this is executed over imputed barcodes, and therefore there is not missingness or het genotypes in the barcodes.
#
markerMap.estimateDiversityMeasures <- function (ctx, barcodeData, distData, measures) {
    result <- c()
    for (mIdx in 1:length(measures)) {
        measure <- measures[mIdx]
        if (measure == "maxHaploFreq") {
            haplos <- apply(barcodeData,1,paste,collapse="")
            value <- max(table(haplos)) / length(haplos)
        } else if (measure == "haploHet") {
            haplos <- apply(barcodeData,1,paste,collapse="")
            value <- pegas::heterozygosity(haplos)
        } else if (measure == "meanSnpHet") {
            hets <- apply(barcodeData, 2, pegas::heterozygosity)
            value <- mean(hets)
        } else if (measure == "medianDistance") {
            mat <- as.matrix(distData)
            mat[lower.tri(mat,diag=TRUE)] <- NA
	    value <- stats::median(mat, na.rm=TRUE)
        } else {
            stop(paste("Invalid diversity measure:", measure))
        }
        result <- c(result, value)
    }
    result
}
