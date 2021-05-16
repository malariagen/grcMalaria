###############################################################################
# Map Haplotype Frequency Analysis
################################################################################
#
# We use visualization types to specify different graphic renditions from the same analyses:
# haploFreq has visualization types "pie" and "bar" (pie and bar markers)
# haploShare has visualization types "pie" and "bar" (pie and bar markers) or "group" to show the prevalence of the different haplotype groups
#
haploMap.execute <- function(ctx, datasetName, analysisName, mapType, aggregation, measures, params) {
    visTypes <- analysis.getParam ("map.haplo.visualizations", params)				#; print(visTypes)
    for (vIdx in 1:length(visTypes)) {
        haploMap.executeVisualization (ctx, datasetName, analysisName, mapType, visTypes[vIdx], aggregation, measures, params)
    }
}
#
haploMap.executeVisualization <- function(ctx, datasetName, analysisName, mapType, visType, aggregation, measures, params) {
    #print(visType)
    dataset <- ctx[[datasetName]]
    
    # Create the information object, used to create the plot
    # Start by getting the output folders
    dataFolder <- getOutFolder(ctx, analysisName, c(paste("map", mapType, sep="-"), "data"))
    info <- list(dataFolder=dataFolder)

    info$analysisName <- analysisName;    info$mapType <- mapType;    info$visType <- visType
    info$plotTitle <- analysisName;

    # Build the map necessary to display these samples
    # Construct a base plot for completing subsequent maps
    baseMapInfo <- info$baseMapInfo <- map.buildBaseMap (ctx, datasetName, analysisName, dataset$meta, dataFolder, params)

    # Now compute the aggregation units, the values to be plotted, and make the map
    for (aggIdx in 1:length(aggregation)) {
        aggLevel <- info$aggLevel <- as.integer(aggregation[aggIdx])    			#; print(aggLevel)
        aggLevelIdx <- aggLevel + 1

        # Get the aggregated data for the aggregation units
        aggUnitData <- info$aggUnitData <- map.getAggregationUnitData (ctx, datasetName, aggLevel, analysisName, mapType, params, dataFolder)	#; print(aggUnitData)

        # Work out a "standard" haplotype marker size (1/25 of the shortest side) and apply a scaling factor
        scalingFactor <- analysis.getParam ("map.haplo.markerScale", params)			#; print(scalingFactor)
        bbox <- baseMapInfo$gadmBB;
        minSide <- min(abs(bbox$yMax-bbox$yMin),abs(bbox$xMax-bbox$xMin))
        stdMarkerSize <- info$stdMarkerSize <- scalingFactor * minSide / 25			#; print(info$stdMarkerSize)
        
        # Compute the number of samples that correspond to this standard size
        stdMarkerCountParam <- analysis.getParam ("map.haplo.markerSampleCount", params)	#; print(stdMarkerCountParam)
        stdMarkerCount <- 0
        if (is.numeric(stdMarkerCountParam)) {
            stdMarkerCount <- stdMarkerCountParam
        } else if (stdMarkerCountParam == "mean") {
            stdMarkerCount <- mean(as.numeric(aggUnitData$SampleCount))				#; print(aggUnitData$SampleCount)
        } else if (stdMarkerCountParam == "none") {
            stdMarkerCount <- 0
        } else {
            stop(paste("Invalid map.stdMarkerCount value:", stdMarkerCountParam))        
        }											#; print(stdPieCount)
        stdMarkerCount <- info$stdMarkerCount <- as.numeric(stdMarkerCount)			#; print(info$stdMarkerCount)
        
        # Do the actual plotting by calling the function we've just defined.
        # The haplotype sharing task produces multiple maps, for different thresholds of haplotype identity 
        if (mapType=="haploFreq") {
            # Get the haplotype counts data
	    info$haploCountData <- haploMap.buildCountData (aggLevel, aggUnitData, dataset, params)	#; print(head(info$haploCountData)
            haploMap.plotMap (ctx, info, params)

        } else if (mapType=="haploShare") {
	    identityLevels <- cluster.getIdentityLevels (params)
            for (mIdx in 1:length(identityLevels)) {
                identityLevel <- info$identityLevel <- identityLevels[mIdx]			#; print (identityLevel)

                # Palette for pie charts- determines the number of distinguishable colurs for the chart
	        haploGroupPalette <- ctx$config$defaultPalette
	        maxGroups <- length(haploGroupPalette)
	    
	        # Get the haplotype sharing data
	        groupData <- cluster.findbyIdentity (ctx, datasetName, analysisName, identityLevel, params)
	        haploShareData <- haploMap.buildSharedCountData (aggLevel, aggUnitData, 
	                                     groupData, dataset$meta, maxGroups, params)	#; print(head(haploShareData)
	        info$groupData <- groupData
	        info$haploShareData <- haploShareData
	    
	        # Get the name of all the haplogroups, order them, trim the palette if necessary 
	        # and add white as the last colour for class "Other" (samples not belionging to a haplo group)
	        haploGroupIds <- unique(as.character(haploShareData$HaploGroup))
	        haploGroupIds <- haploGroupIds[order(haploGroupIds)]				#; print(haploGroupIds)
	        haploGroupPalette <- haploGroupPalette[1:(length(haploGroupIds)-1)]
	        haploGroupPalette <- c(haploGroupPalette,"white")
	        names(haploGroupPalette) <- haploGroupIds
	        info$palette <- haploGroupPalette

                if (visType == "group") {
                    # Produce a set of per-haplotype maps to show where the haplotype goup circulates
                    for (gIdx in 1:(length(haploGroupIds)-1)) {				#; print(gIdx)	        
                        haploGroup <- haploGroupIds[gIdx]				#; print(haploGroup)
                        info$haploGroup <- haploGroup
                        info$plotTitle <- paste(analysisName, "-", haploGroup);
                        hgHaploShareData <- info$haploShareData <- haploShareData[which(haploShareData$HaploGroup == haploGroup),]	#; print(nrow(hgHaploShareData))
                        hgAggUnits <- as.character(hgHaploShareData$UnitId)		#; print(hgAggUnits)
                        hgAggUnitData <- info$aggUnitData <- aggUnitData[hgAggUnits,]	#; print(nrow(hgAggUnitData))
                        haploMap.plotMap (ctx, info, params)
                    }
                } else {
                    # Create the plot with the haplosharing markers (pies or bars)
                    haploMap.plotMap (ctx, info, params)
                }
            }
        }
    }
}

# This function cotains the maincode to produce the map plot.
# It will plot the same background map, with different types of markers accordong to map and visualization options
# It has been separated so it can be called multiple times in the case of haplo sharing
haploMap.plotMap <- function (ctx, info, params) {

    # Start with the background map
    baseMapInfo <- info$baseMapInfo
    mapPlot <- baseMapInfo$baseMap

    # If we need to show aggregation unit names, we need to compute the label positioning and plot before the markers
    aggUnitData <- info$aggUnitData
    aggLevel <- info$aggLevel
    aggLevelIdx <- aggLevel+1
    showMarkerNames <- analysis.getParam ("map.markerNames", params)
    if (showMarkerNames) {
        aggColName <- map.getAggregationColumns(aggLevelIdx)				#; print(aggColName)
        lp <- map.computeLabelParams (aggUnitData, aggColName, baseMapInfo)		#; print(lp)
        mapPlot <- mapPlot + 
                   ggplot2::geom_point(ggplot2::aes(x=lon, y=lat), data=lp, colour="red") +
                   ggrepel::geom_label_repel(ggplot2::aes(x=lon, y=lat, label=label), data=lp,
                                    fill="black", size=4.5, fontface="bold", color="darkgray", show.legend=FALSE,
                                    hjust=lp$just, vjust=0.5, nudge_x=lp$x, nudge_y=lp$y, label.padding=grid::unit(0.2, "lines"))
    }

    # Now add the markers
    mapType <- info$mapType
    visType <- info$visType
    stdMarkerCount <- info$stdMarkerCount
    stdMarkerSize <- info$stdMarkerSize
    if (info$mapType == "haploShare") {
        haploShareData <- info$haploShareData
        palette <- info$palette
    }
    
    if (info$mapType=="haploFreq") {
        haploCountData <- info$haploCountData
        mapPlot <- haploMap.addFreqMarkers (mapPlot, visType, haploCountData, aggUnitData, stdMarkerCount, stdMarkerSize)
    } else if (info$mapType=="haploShare") {
        if (info$visType=="group") {
            haploGroup <- info$haploGroup		#; print(haploGroup)
            mapPlot <- haploMap.addConnections (mapPlot, visType, haploGroup, haploShareData, palette, stdMarkerCount, stdMarkerSize)

            groupData <- info$groupData
            groupInfoText <- cluster.getClusterStatsText (ctx, haploGroup, groupData)	#; print(groupInfoText)
            bb <- baseMapInfo$gadmBB
            #mapPlot <- mapPlot +
            #           ggplot2::annotate("label", x=bb$xMin, y=bb$yMax, hjust="left", vjust="top", label=groupInfoText)
            #
            # Dirty trick so we can show an annotation with group info above the legends.
            # We "plot" a couple of points, and create a ficticious alpha legend containing the group info text.
            dummydf <- data.frame(x=c(bb$xMin,bb$xMin), y=c(bb$yMax,bb$yMax), alpha=c(0.1, 0.11))
            mapPlot <- mapPlot +
                       ggplot2::geom_point(ggplot2::aes(x=x, y=y, alpha=as.numeric(alpha), size=0.01), data=dummydf) +
                       ggplot2::scale_alpha_continuous(haploGroup, breaks=c(0.1, 0.11), labels=c(groupInfoText,""), 
                                                        guide=ggplot2::guide_legend(order=1,keywidth=0, keyheight=0.01, nrow=1,
                                                        override.aes=list(shape=NA,fill=NA,size=0.01)))
        } else {
            mapPlot <- haploMap.addShareMarkers (mapPlot, visType, haploShareData, palette, aggUnitData, stdMarkerCount, stdMarkerSize)
        }
    }
    # Now add the decorative elements
    mapPlot <- mapPlot +
    	       ggplot2::labs(title=info$plotTitle, subtitle="")+
               ggplot2::theme(plot.title=ggplot2::element_text(face="bold", size=ggplot2::rel(1.2), hjust=0.5),
                     panel.background=ggplot2::element_rect(colour=NA),
                     plot.background=ggplot2::element_rect(colour=NA),
                     axis.title=ggplot2::element_text(face="bold",size=ggplot2::rel(1)),
                     axis.title.y=ggplot2::element_text(angle=90,vjust=2),
                     axis.title.x=ggplot2::element_text(vjust=-0.2))
     
    # Save to file. the size in inches is given in the params.
    mapSize  <- analysis.getParam ("map.size", params)

    plotFolder <- getOutFolder(ctx, info$analysisName, c(paste("map", mapType, sep="-"), "plots"))
    aggLabel <- map.getAggregationLabels(aggLevel)
    graphicFilename  <- paste("map", info$analysisName, aggLabel, mapType, visType, sep="-")
    if (info$mapType == "haploShare") {
        levelLabel <- cluster.getIdentityLevelLabel(info$identityLevel)
        graphicFilename  <- paste(graphicFilename, levelLabel, sep="-")
        if (info$visType == "group") {
            plotFolder <- getOutFolder(ctx, info$analysisName, c("map-haploGroups", "plots", levelLabel))
            graphicFilename  <- paste(graphicFilename, haploGroup, sep="-")
        }
    }
    graphicFilename  <- paste(plotFolder, paste(graphicFilename,"png",sep="."), sep="/")
    ggplot2::ggsave(plot=mapPlot, filename=graphicFilename, device="png", width=mapSize[1], height=mapSize[2], units="in", dpi=300)
}

###############################################################################
# Haplotype Frequency Markers plotting 
################################################################################
#
haploMap.addFreqMarkers <- function (mapPlot, visType, countData, aggUnitData, stdMarkerCount, stdMarkerSize) {
    if (visType=="bar") {
        mapPlot <- haploMap.addFreqBars (mapPlot, countData, aggUnitData, stdMarkerCount, stdMarkerSize)
    } else if (visType=="pie") {
        mapPlot <- haploMap.addFreqPies (mapPlot, countData, aggUnitData, stdMarkerCount, stdMarkerSize)
    }
    mapPlot
}

haploMap.addFreqPies <- function (mapPlot, countData, aggUnitData, stdMarkerCount, stdMarkerSize) {
    # Now add the pie chart markers
    if (stdMarkerCount==0) {
        mapPlot <- mapPlot +
                   ggforce::geom_arc_bar(ggplot2::aes(x0=Longitude, y0=Latitude, r0=0, r=stdMarkerSize, 
                                         fill=Haplo, amount=HaploCount),
                             data=countData, stat="pie", inherit.aes=FALSE,
                             colour="gray25", stroke=0.5, fill="white", show.legend=FALSE)
    } else {
        mapPlot <- mapPlot +
                   ggforce::geom_arc_bar(ggplot2::aes(x0=Longitude, y0=Latitude, r0=0, r=stdMarkerSize*sqrt(SampleCount/stdMarkerCount), 
                                         fill=Haplo, amount=HaploCount),
                             data=countData, stat="pie", inherit.aes=FALSE,
                             colour="gray25", stroke=0.5, fill="white", show.legend=FALSE)
    }
    mapPlot
}

haploMap.addFreqBars <- function (mapPlot, countData, aggUnitData, stdMarkerCount, stdMarkerSize) {

    # Get all aggregation unit ids, in descending order of sample count
    aggUnitData <- aggUnitData[order(-aggUnitData$SampleCount),]
    aggUnits <- rownames(aggUnitData)							#; print(aggUnits)
    freqBarData <- NULL
    for (aIdx in 1:length(aggUnits)) {
        aggUnit <- aggUnits[aIdx]
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
        
        # Get all the haplos for this unit, ordered by sample count
        unitHaploData <- countData[which(countData$UnitId==aggUnit),]
        unitHaploData <- unitHaploData[order(unitHaploData$HaploCount),]
        unitRowCount <- nrow(unitHaploData)
        
        # For each haplo, worh out the y boundaries
        y2 <- vector(mode="numeric", length=unitRowCount)
        y <- y0
        for (i in 1:unitRowCount) {
           y <- y + (unitHaploData$HaploCount[i] * sHeight)
           y2[i] <- y
        }
        y1 <- c(y0, y2[1:(unitRowCount-1)])
        x1 <- rep(x1, unitRowCount)
        x2 <- rep(x2, unitRowCount)
        unitBarData <- cbind(unitHaploData, x1=x1, x2=x2, y1=y1, y2=y2)
        freqBarData <- rbind(freqBarData, unitBarData)
    }				
    #print(freqBarData)
    mapPlot <- mapPlot +
               ggplot2::geom_rect(ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), data=freqBarData, inherit.aes=FALSE,
                                  colour="gray25", size=0.5, fill="white", show.legend=FALSE)
    mapPlot
}

#
# For each aggregation unit, we get a count of each unique haplotype, ordered in descending count
#
haploMap.buildCountData <- function(aggLevel, aggUnitData, dataset, params) {
    sampleMeta   <- dataset$meta
    barcodeData  <- dataset$barcodes

    # Get all aggregation unit ids
    aggUnits <- rownames(aggUnitData)	#; print(aggUnits)
    
    # Create aggregation index for each sample (the id of the aggregation unit where the sample originates)
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta, params)
    haplos <- apply(barcodeData,1,paste,collapse="")

    # Get the data for all aggregation units
    countData <- NULL
    for (aIdx in 1:length(aggUnits)) {
        # Get the sample data to be aggregated for this unit
        aggUnit <- aggUnits[aIdx]
        unitHaplos <- haplos[which(aggIndex == aggUnit)]		#; print(nrow(aggHaplos))
        unitHaploCounts <- as.integer(table(unitHaplos))
        unitHaploCounts <- unitHaploCounts[order(-unitHaploCounts)]
        unitHaploNum <- length (unitHaploCounts)

        unitData <- aggUnitData[aIdx,]
        df <- data.frame(matrix(unitData,ncol=ncol(aggUnitData),nrow=unitHaploNum,byrow=TRUE))
        colnames(df) <- colnames(aggUnitData)
        df$Haplo <- paste("Haplo",c(1:unitHaploNum),sep="_")
        df$HaploCount <- unitHaploCounts
        countData <- rbind(countData, df)
    }
    countData$Haplo <- factor(countData$Haplo)
    countData$Latitude <- as.numeric(countData$Latitude)
    countData$Longitude <- as.numeric(countData$Longitude)
    countData$SampleCount <- as.numeric(countData$SampleCount)
    countData
}

###############################################################################
# Haplotype Sharing Markers plotting 
################################################################################
#
haploMap.addShareMarkers <- function (mapPlot, visType, haploShareData, haploGroupPalette, aggUnitData, stdMarkerCount, stdMarkerSize) {
    # Do th actual marker plotting for the twp types of markers
    if (visType=="bar") {
        mapPlot <- haploMap.addShareBars (mapPlot, haploShareData, haploGroupPalette, aggUnitData, stdMarkerCount, stdMarkerSize)
    } else if (visType=="pie") {
        mapPlot <- haploMap.addSharePies (mapPlot, haploShareData, haploGroupPalette, aggUnitData, stdMarkerCount, stdMarkerSize)
    }
    mapPlot
}

haploMap.addSharePies <- function (mapPlot, haploShareData, haploGroupPalette, aggUnitData, stdMarkerCount, stdMarkerSize) {

    # Now add the pie chart markers
    if (stdMarkerCount==0) {
        mapPlot <- mapPlot +
                   ggforce::geom_arc_bar(ggplot2::aes(x0=Longitude, y0=Latitude, r0=0, r=stdMarkerSize, 
                                         fill=HaploGroup, amount=HaploCount),
                                data=haploShareData, stat="pie", inherit.aes=FALSE,
                                colour="gray25", stroke=0.5, show.legend=FALSE) +
                   ggplot2::scale_fill_manual(values=haploGroupPalette)
    } else {
        mapPlot <- mapPlot +
                   ggforce::geom_arc_bar(ggplot2::aes(x0=Longitude, y0=Latitude, r0=0, r=stdMarkerSize*sqrt(SampleCount/stdMarkerCount), 
                                         fill=HaploGroup, amount=HaploCount),
                                data=haploShareData, stat="pie", inherit.aes=FALSE,
                                colour="gray25", stroke=0.5, show.legend=FALSE) +
                   ggplot2::scale_fill_manual(values=haploGroupPalette)
    }
    mapPlot
}


haploMap.addShareBars <- function (mapPlot, haploShareData, haploGroupPalette, aggUnitData, stdMarkerCount, stdMarkerSize) {

    # Get all aggregation unit ids, in descending order of sample count
    aggUnitData <- aggUnitData[order(-aggUnitData$SampleCount),]
    
    aggUnits <- rownames(aggUnitData)							#; print(aggUnits)
    freqBarData <- NULL
    for (aIdx in 1:length(aggUnits)) {
        aggUnit <- aggUnits[aIdx]
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
        
        # Get all the haplos for this unit, ordered by Haplo Group name
        unitHaploData <- haploShareData[which(haploShareData$UnitId==aggUnit),]		#; print(unitHaploData)
        unitHaploGroups <- as.character(unitHaploData$HaploGroup)			#; print(unitHaploGroups)
        unitHaploData <- unitHaploData[order(unitHaploGroups, decreasing=TRUE),]	#; print(unitHaploData)
        unitRowCount <- nrow(unitHaploData)
        
        # For each haplo, worh out the y boundaries
        y2 <- vector(mode="numeric", length=unitRowCount)
        y <- y0
        for (i in 1:unitRowCount) {
           y <- y + (unitHaploData$HaploCount[i] * sHeight)
           y2[i] <- y
        }
        y1 <- c(y0, y2[1:(unitRowCount-1)])
        x1 <- rep(x1, unitRowCount)
        x2 <- rep(x2, unitRowCount)
        unitBarData <- cbind(unitHaploData, x1=x1, x2=x2, y1=y1, y2=y2)
        freqBarData <- rbind(freqBarData, unitBarData)
    }											#; print(freqBarData)
    
    # Got the dataframe for drawing all the rectangles, now do the drawing and colouring
    mapPlot <- mapPlot +
               ggplot2::geom_rect(ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=HaploGroup), 
                                  data=freqBarData, inherit.aes=FALSE,
                                  colour="gray25", size=0.5, show.legend=FALSE) +
               ggplot2::scale_fill_manual(values=haploGroupPalette)
    mapPlot
}

haploMap.buildSharedCountData <- function(aggLevel, aggUnitData, groupData, sampleMeta, maxGroups, params) {
    # Get the subgraph membership file for this identity threshold, which is created by the clusteing code
    # (shared with the "graph" analysis task)
    groupMemberData <- cluster.getMemberData (groupData)

    groupData <- groupData[order(groupData$Cluster),]
    if (nrow(groupData) > maxGroups) {
        groupData <- groupData[1:maxGroups,]
    }
    groupIds <- as.character(groupData$Cluster)						#; print(groupIds)
    groupMemberData <- groupMemberData[which(groupMemberData$Cluster %in% groupIds),]
    groupMembers <- groupMemberData$Sample						#; print(groupMembers)
    memberGroups <- groupMemberData$Cluster						#; print(memberGroups)
   
    # Assign group Ids to the sample metadata
    sampleNames <- rownames(sampleMeta)
    hGroup <- rep("-", length(sampleNames))
    names(hGroup) <- sampleNames
    hGroup[groupMembers] <- memberGroups
    sampleMeta$HaploGroup <- hGroup
    
    # Get all aggregation unit ids
    aggUnits <- rownames(aggUnitData)	#; print(aggUnits)
    
    # Create aggregation index for each sample (the id of the aggregation unit where the sample originates)
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta, params)

    # Get the data for all aggregation units
    countData <- NULL
    for (aIdx in 1:length(aggUnits)) {
        # Get the sample data to be aggregated for this unit
        aggUnit <- aggUnits[aIdx]
        unitMeta <- sampleMeta[which(aggIndex == aggUnit),]		#; print(nrow(unitMeta))
        gpUnitMeta <- unitMeta[which(unitMeta$HaploGroup != "-"),]
        unitGroupCounts <- table(gpUnitMeta$HaploGroup)
        unitGroupCounts <- unitGroupCounts[order(-unitGroupCounts)]
        unitGroupNum <- length(unitGroupCounts)+1

        unitData <- aggUnitData[aIdx,]
        df <- data.frame(matrix(unitData,ncol=ncol(aggUnitData),nrow=unitGroupNum,byrow=TRUE))
        colnames(df) <- colnames(aggUnitData)
        df$HaploGroup <- c(names(unitGroupCounts), "Other")
        df$HaploCount <- c(unitGroupCounts, nrow(unitMeta)-nrow(gpUnitMeta))
        df$HaploProp <- df$HaploCount / unitData$SampleCount
        countData <- rbind(countData, df)
    }
    countData$HaploGroup <- factor(countData$HaploGroup)
    countData$Latitude <- as.numeric(countData$Latitude)
    countData$Longitude <- as.numeric(countData$Longitude)
    countData$SampleCount <- as.numeric(countData$SampleCount)
    countData
}

###############################################################################
# Maps of connections between Haplotype Sharing sites plotting 
################################################################################
#
haploMap.addConnections <- function (mapPlot, visType, haploGroup, haploShareData, haploGroupPalette, stdMarkerCount, stdMarkerSize) {
    # Work out the values to display for every marker (the fractional prevalence)
    mValues <- as.numeric(haploShareData$HaploProp)			#; print(haploGroup); print(haploShareData)
    valueLabels <- round(mValues, digits=2)
    
    # Work out size, colour and colour gradient for the markers
    pointSizes <- 16
    haploColour <- haploGroupPalette[haploGroup]
    mMax <- max(mValues); scaleMax <- round(mMax,digits=1); scaleMax <- ifelse(scaleMax<mMax, scaleMax+0.1, scaleMax)
    mMin <- min(mValues); scaleMin <- round(mMin,digits=1); scaleMin <- ifelse(scaleMin>mMin, scaleMin-0.1, scaleMin)  #; print(paste(scaleMin, scaleMax))

    # For ech pair, work out a connectedness measure, in this case the fraction of sample pairs in which both samples carry this haplo     
    aggUnitCount <- nrow(haploShareData)
    if (aggUnitCount > 1) {
        connectData <- NULL
        for (a1Idx in 1:(aggUnitCount-1)) {
            for (a2Idx in (a1Idx+1):aggUnitCount) {
                # Get the sample data
                lat1 <- haploShareData$Latitude[a1Idx]
                lon1 <- haploShareData$Longitude[a1Idx]
                sampleCount1 <- haploShareData$SampleCount[a1Idx]
                haploCount1 <- haploShareData$HaploCount[a1Idx]
                #
                lat2 <- haploShareData$Latitude[a2Idx]
                lon2 <- haploShareData$Longitude[a2Idx]
                sampleCount2 <- haploShareData$SampleCount[a2Idx]
                haploCount2 <- haploShareData$HaploCount[a2Idx]
                #
                pairCount <- sampleCount1 * sampleCount2
                matchCount <- haploCount1 * haploCount2
                prop <- matchCount / pairCount
            
                cValues <- c(lat1, lon1, lat2, lon2, prop)
                connectData <- rbind(connectData, cValues)
            }
        }
        aggUnitPairData <- data.frame(connectData)					#; print (dim(aggUnitPairData))
        nc <- ncol(aggUnitPairData)
        aggUnitPairData[,1:nc] <- sapply(aggUnitPairData[,1:nc], as.numeric)		#; print(ncol(aggUnitPairData))
        colnames(aggUnitPairData) <- c("Lat1", "Lon1", "Lat2", "Lon2", "Weight")	#; print(aggUnitPairData)
        
        hc <- grDevices::rgb2hsv(grDevices::col2rgb(haploColour))
        weakHaploColour <- grDevices::hsv(h=hc[1,1],s=hc[2,1]*0.2,v=hc[3,1])
	
        # Now plot the connections
        mapPlot <- mapPlot +
                   ggplot2::geom_curve(ggplot2::aes(x=Lon1, y=Lat1, xend=Lon2, yend=Lat2, size=Weight, colour=Weight),	# draw edges as arcs
                               data=aggUnitPairData, curvature=0.2, alpha=0.75) +
                   ggplot2::scale_size_continuous(guide=FALSE, range=c(0.25, 4)) +          			# scale for edge widths
                   ggplot2::scale_colour_gradientn(name="Sharing", colours=c(weakHaploColour,haploColour))
    }
     
    # Plot the Aggregation Unit markers
    mapPlot <- mapPlot +
	       ggplot2::geom_point(ggplot2::aes(x=Longitude, y=Latitude, fill=HaploProp), 
	                           data=haploShareData, size=pointSizes, shape=21, stroke=2) +
	       ggplot2::scale_fill_gradientn(name="Prevalence", limits=c(scaleMin,scaleMax), colours=c("white",haploColour), values=c(0,1)) +
               ggplot2::geom_text(ggplot2::aes(x=Longitude, y=Latitude), 
                                  data=haploShareData, label=valueLabels, hjust=0.5, vjust=0.5, size=4.5, fontface="bold")
    mapPlot
}

