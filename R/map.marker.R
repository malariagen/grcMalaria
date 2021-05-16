###############################################################################
# Map Aggregated Measure Analysis
################################################################################
#
markerMap.getDiversityMeasures <- function() {
    c("maxHaploFreq",
      "haploHet1",
      "haploHet2",
      "haploHet3",
      "meanSnpHet",
      "medianDistance")
}

markerMap.execute <- function(ctx, datasetName, analysisName, mapType, aggregation, measures, params) {
    config <- ctx$config
    dataset <- ctx[[datasetName]]

    # Get the output folders
    dataFolder <- getOutFolder(ctx, analysisName, c(paste("map", mapType, sep="-"), "data"))
    # Get the sample metadata
    sampleMeta <- dataset$meta
    # Build the map necessary to display these samples
    # Construct a base plot for completing subsequent maps
    baseMapInfo <- map.buildBaseMap (ctx, datasetName, analysisName, sampleMeta, dataFolder, params)

    # Now compute the aggregation units, the values to be plotted, and make the map
    for (aggIdx in 1:length(aggregation)) {
        aggLevel <- as.integer(aggregation[aggIdx])        					#; print(aggLevel)
        aggLevelIdx <- aggLevel + 1

        # Get the aggregated data for the aggregation units
        aggUnitData <- map.getAggregationUnitData (ctx, datasetName, aggLevel, analysisName, mapType, params, dataFolder)	#; print(aggUnitData)
        measures <- markerMap.checkMeasures (ctx, mapType, measures)
        aggUnitData <- markerMap.estimateMeasures (ctx, datasetName, aggLevel, aggUnitData, analysisName, mapType, measures, params, dataFolder)	#; print(aggUnitData)

        for (mIdx in 1:length(measures)) {
            measure <- measures[mIdx]		           		#; print(measure)
            
            # Select the aggregation units to be plotted
            # In this case, those that have an NA for the measure being plotted (should be only for drug resistance)
            mValues <- aggUnitData[,measure]				#; print(mValues)
            selAggUnitData <- aggUnitData[which(!is.na(mValues)),]	#; print(nrow(selAggUnitData))
            
            # Compute marker sizes. 
            pointSizes <- markerMap.getAggUnitMarkerSizes (selAggUnitData, params)

            # Do the actual plot, starting with the background map
            mapPlot <- baseMapInfo$baseMap
            
            # If we need to show aggregation unit names, we need to compute the label positioning and plot before the markers
            showMarkerNames <- analysis.getParam ("map.markerNames", params)
            if (showMarkerNames) {
                aggColName <- c("Country","AdmDiv1","AdmDiv2")[aggLevelIdx]
                lp <- map.computeLabelParams (selAggUnitData, aggColName, baseMapInfo)
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

            # Now add the markers, coloured according to the appropriate scale, depending on the type of map
            mValues <- selAggUnitData[,measure]
            if (mapType=="diversity") {
                markerColours <- analysis.getParam ("map.diversity.markerColours", params)
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
            } else if (mapType=="drug") {
                mapPlot <- mapPlot +
	              ggplot2::geom_point(data=selAggUnitData, aes_string2(x="Longitude", y="Latitude", fill=measure), size=pointSizes, shape=21, stroke=2) +
	              ggplot2::scale_fill_gradientn(limits=c(0,1), colours=c("green3","orange2","red3","red3"), values=c(0, 0.2, 0.75, 1))
            } else if (mapType=="mutation") {
                markerColours <- analysis.getParam ("map.prevalence.markerColours", params)
                # Two marker colours can be specified to create a gradient. If a single marker colour is 
                # specified, then create a gradient from white to that colour.
                if (length(markerColours) == 1) {
                    markerColours <- c("white", markerColours)
                }
                scaleMin <- 0; scaleMax <- 1
                mapPlot <- mapPlot +
	              ggplot2::geom_point(data=selAggUnitData, aes_string2(x="Longitude", y="Latitude", fill=measure), size=pointSizes, shape=21, stroke=2) +
	              ggplot2::scale_fill_gradientn(limits=c(scaleMin,scaleMax), colours=markerColours, values=c(0,1))
	          } 
	    
            # Now add the decorative elements
            valueLabels <- round(mValues, digits=2)
            mapPlot <- mapPlot +
                          ggplot2::geom_text(data=selAggUnitData, ggplot2::aes_string(x="Longitude", y="Latitude"), label=valueLabels, hjust=0.5, vjust=0.5, size=4.5, fontface="bold") +
                          ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",size = ggplot2::rel(1.2), hjust = 0.5),
                          panel.background = ggplot2::element_rect(colour = NA),
                          plot.background = ggplot2::element_rect(colour = NA),
                          axis.title = ggplot2::element_text(face = "bold",size = ggplot2::rel(1)),
                          axis.title.y = ggplot2::element_text(angle=90,vjust =2),
                          axis.title.x = ggplot2::element_text(vjust = -0.2))
	    
            # Save to file. the size in inches is given in the config.
            mapSize  <- analysis.getParam ("map.size", params)
            plotFolder <- getOutFolder(ctx, analysisName, c(paste("map", mapType, sep="-"), "plots"))
            aggLabel <- map.getAggregationLabels(aggLevel)
            graphicFilenameRoot  <- paste(plotFolder, paste("map", analysisName, aggLabel, measure, sep="-"), sep="/")
            ggplot2::ggsave(plot=mapPlot, filename=paste(graphicFilenameRoot,"png",sep="."), device="png", width=mapSize[1], height=mapSize[2], units="in", dpi=300)
        }
    }
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
    mSizeParam <- analysis.getParam ("map.markerSize", params)
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
markerMap.checkMeasures <- function(ctx, mapType, measures) {
    # Get the admin division values from the first sample of this unit (assuming the values are the same for all)
    if (mapType=="diversity") {
        if ("ALL" %in% measures) {
            measures <- markerMap.getDiversityMeasures()
        }
    } else if (mapType=="drug") {
        if ("ALL" %in% measures) {
            measures <- ctx$config$drugs
        }
    } else if (mapType=="mutation") {
        if ("ALL" %in% measures) {
            measures <- ctx$config$drugResistanceMutations
        }
    }
    measures
}

markerMap.estimateMeasures <- function(ctx, datasetName, aggLevel, aggUnitData, analysisName, mapType, measures, params, dataFolder) {

    dataset <- ctx[[datasetName]]
    sampleMeta   <- dataset$meta
    barcodeData  <- dataset$barcodes
    distData     <- dataset$distance

    # Create aggregation index for each sample (the id of the aggregation unit where the sample originates)
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta, params)
    
    # Get all aggregation units
    aggUnits <- rownames(aggUnitData)	#; print(aggUnits)
    
    # Get the data for all aggregation units
    measureData <- NULL
    for (aIdx in 1:length(aggUnits)) {

        # Get the sample data to be aggregated for this unit
        aggUnit <- aggUnits[aIdx]
        aggSamplesMeta <- sampleMeta[which(aggIndex == aggUnit),]		#; print(nrow(aggSamplesMeta))
        aggSamples <- rownames(aggSamplesMeta)
        aggBarcodes <- barcodeData[aggSamples,]
        aggDist <- distData[aggSamples,aggSamples]
        
        # Get the admin division values from the first sample of this unit (assuming the values are the same for all)
        if (mapType=="diversity") {
            cValues <- markerMap.estimateDiversityMeasures (ctx, aggBarcodes, aggDist, measures)
        } else if (mapType=="drug") {
            cValues <- meta.getResistancePrevalence (ctx, aggSamplesMeta, measures, params)
        } else if (mapType=="mutation") {
            cValues <- markerMap.estimateMutationPrevalence (ctx, aggSamplesMeta, measures, params)
        }                                                           #; print(cValues)
        measureData <- rbind(measureData, cValues)
    }
    measureData <- data.frame(measureData)
    measureData <- sapply(measureData, as.numeric)
    cNames <- c(colnames(aggUnitData), measures)
    aggUnitData <- cbind(aggUnitData, measureData)
    colnames(aggUnitData) <- cNames
    
    # Write out the aggregation unit data to file
    aggDataFilename  <- paste(dataFolder, "/AggregatedData-", analysisName, "-", aggLevel, ".tab", sep="")
    utils::write.table(aggUnitData, file=aggDataFilename, sep="\t", quote=FALSE, row.names=FALSE)

    aggUnitData
}

markerMap.estimateDiversityMeasures <- function (ctx, barcodeData, distData, measures) {
    result <- c()
    for (mIdx in 1:length(measures)) {
        measure <- measures[mIdx]
        if (measure == "maxHaploFreq") {
            haplos <- apply(barcodeData,1,paste,collapse="")
            value <- max(table(haplos)) / length(haplos)
        } else if (measure == "haploHet1") {
            haplos <- apply(barcodeData,1,paste,collapse="")
            value <- pegas::heterozygosity(haplos)
        } else if (measure == "haploHet2") {
            haplos <- apply(barcodeData,1,paste,collapse="")
            hets <- sapply(1:1000, function(i) pegas::heterozygosity(haplos[sample(1:length(haplos), 10)]))
            value <- mean(hets)
        } else if (measure == "haploHet3") {
            sCount <- ncol(barcodeData)
            hets <- sapply(1:100, function(i) pegas::heterozygosity(apply(barcodeData[,sample(1:sCount,10)], 1, paste, collapse="")))
            #hets <- sapply(1:100, computeHaploHet3, barcodeData)
            value <- mean(hets)
        } else if (measure == "meanSnpHet") {	# Mean SNP Het == Mean Distance
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

markerMap.estimateMutationPrevalence <- function (ctx, sampleMeta, mutationNames, params) {
    positionNames <- c()
    alleles <- c()
    for (mIdx in 1:length(mutationNames)) {
        mut <- mutationNames[mIdx]
        pos <- substring(mut, 1, nchar(mut)-1)
        all <- substring(mut, nchar(mut), nchar(mut))
        positionNames <- c(positionNames, pos)
        alleles <- c(alleles, all)
    }										#; print(positionNames); print(alleles)
    prevData <- meta.getAllelePrevalence (ctx, sampleMeta, positionNames, alleles)	#; print(prevData)
    prevData
}
