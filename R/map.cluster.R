###############################################################################
# Map Cluster Analysis
################################################################################
#
# We use visualization types to specify different graphic renditions from the same analyses:
# barcodeFrequency has visualization types "pie" and "bar" (pie and bar markers)
# clusterPrevalence has visualization type "cluster" to show the prevalence of the different clusters
#
clusterMap.execute <- function(userCtx, datasetName, sampleSetName, mapType, aggregation, measures, params) {
    visTypes <- analysis.getParam ("map.cluster.visualizations", params)			#; print(visTypes)
    for (vIdx in 1:length(visTypes)) {
        clusterMap.executeVisualization (userCtx, datasetName, sampleSetName, mapType, visTypes[vIdx], aggregation, measures, params)
    }
}
#
clusterMap.executeVisualization <- function(userCtx, datasetName, sampleSetName, mapType, visType, aggregation, measures, params) {
    sampleSet <- userCtx$sampleSets[[sampleSetName]]
    ctx <- sampleSet$ctx
    dataset <- ctx[[datasetName]]
    
    # Create the information object, used to create the plot
    # Start by getting the output folders
    dataFolder <- getOutFolder(ctx, sampleSetName, c(paste("map", mapType, sep="-"), "data"))
    info <- list(dataFolder=dataFolder)

    info$sampleSetName <- sampleSetName;    info$mapType <- mapType;    info$visType <- visType
    info$plotTitle <- sampleSetName;

    # Build the map necessary to display these samples
    # Construct a base plot for completing subsequent maps
    baseMapInfo <- info$baseMapInfo <- map.buildBaseMap (ctx, datasetName, sampleSetName, dataset$meta, dataFolder, params)

    # Now compute the aggregation units, the values to be plotted, and make the map
    for (aggIdx in 1:length(aggregation)) {
        aggLevel <- info$aggLevel <- as.integer(aggregation[aggIdx])    			#; print(aggLevel)
        aggLevelIdx <- aggLevel + 1

        # Get the aggregated data for the aggregation units
        aggUnitData <- info$aggUnitData <- map.getAggregationUnitData (ctx, datasetName, aggLevel, sampleSetName, mapType, params, dataFolder)	#; print(aggUnitData)

        # Work out a "standard" marker size (1/25 of the shortest side) and apply a scaling factor
        scalingFactor <- analysis.getParam ("map.cluster.markerScale", params)			#; print(scalingFactor)
        bbox <- baseMapInfo$gadmBB;
        minSide <- min(abs(bbox$yMax-bbox$yMin),abs(bbox$xMax-bbox$xMin))
        stdMarkerSize <- info$stdMarkerSize <- scalingFactor * minSide / 25			#; print(info$stdMarkerSize)
        
        # Compute the number of samples that correspond to this standard size
        stdMarkerCountParam <- analysis.getParam ("map.cluster.markerSampleCount", params)	#; print(stdMarkerCountParam)
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
        if (mapType=="barcodeFrequency") {
            # Get the sample counts data
	    info$clusterCountData <- clusterMap.buildCountData (aggLevel, aggUnitData, dataset, params)	#; print(head(info$clusterCountData)
            clusterMap.plotMap (ctx, info, params)

        } else if (mapType %in% c("clusterSharing", "clusterPrevalence")) {
            # Get the cluster data
            clusterSetName  <- analysis.getParam ("cluster.clusterSet.name", params)	#; print(clusterSetName)
            clusterSetInfos <- cluster.getClustersSetFromContext (userCtx, sampleSetName, clusterSetName)

	    # Process the clusters for different thresholds of min identity 
            setCount <- length(clusterSetInfos)
            for (idIdx in 1:setCount) {
                clusterSetInfo <- clusterSetInfos[[idIdx]]
                if (is.null(clusterSetInfo$clusters)) {	# No clusters present at this identity level
                    next
                }

                # We have clusters, determine how they are shared between aggregation units
                clusterData <- clusterSetInfo$clusters
                clusterShareData <- clusterMap.buildSharedCountData (aggLevel, aggUnitData,
                                               clusterData, dataset$meta, params)	#; print(head(clusterShareData))
                info$clusterShareData <- clusterShareData
                info$minIdentity <- clusterSetInfo$minIdentity				#; print(info$minIdentity)
            
                # Construct a palette using the default palette, recycling it if there are too many clusters
                # and adding white as the last colour for class "Other" (samples not belionging to a cluster)
                clusterIds <- unique(as.character(clusterShareData$Cluster))		#; print(clusterIds)
                clusterPalette <- cluster.getClusterPalette (ctx, clusterIds)		#; print(clusterPalette)
                info$palette <- clusterPalette

                if (visType == "cluster") {
                    # Produce a set of per-cluster maps to show where the cluster circulates
                    for (clIdx in 1:length(clusterIds)) {				#; print(clIdx)	        
                        cluster <- clusterIds[clIdx]					#; print(cluster)
                        if (cluster == "Other") {
                            next
                        }
                        info$cluster <- cluster
                        info$plotTitle <- paste(sampleSetName, "-", cluster);
                        info$clusterSetInfo <- clusterSetInfo
                        clClusterShareData <- clusterShareData[which(clusterShareData$Cluster == cluster),]
                        info$clusterShareData <- clClusterShareData			#; print(nrow(clClusterShareData))
                        clAggUnits <- as.character(clClusterShareData$UnitId)		#; print(clAggUnits)
                        info$aggUnitData <- aggUnitData[clAggUnits,]			#; print(nrow(info$aggUnitData))
                        clusterMap.plotMap (ctx, info, params)
                    }
                } else {
                    # Create the plot with the cluster sharing markers (pies or bars)
                    clusterMap.plotMap (ctx, info, params)
                }
            }
        }
    }
}
#
# This function cotains the maincode to produce the map plot.
# It will plot the same background map, with different types of markers accordong to map and visualization options
# It has been separated so it can be called multiple times in the case of cluster sharing
clusterMap.plotMap <- function (ctx, info, params) {

    # Start with the background map
    baseMapInfo <- info$baseMapInfo
    mapPlot <- baseMapInfo$baseMap

    # If we need to show aggregation unit names, we need to compute the label positioning and plot before the markers
    aggUnitData <- info$aggUnitData
    aggLevel <- info$aggLevel
    showMarkerNames <- analysis.getParam ("map.markerNames", params)
    if (showMarkerNames) {
        lp <- map.computeLabelParams (aggUnitData, baseMapInfo)		#; print(lp)
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
    if (mapType %in% c("clusterSharing", "clusterPrevalence")) {
        clusterShareData <- info$clusterShareData
        palette <- info$palette
    }
    
    clusterSetInfo <- info$clusterSetInfo
    if (info$mapType=="barcodeFrequency") {
        clusterCountData <- info$clusterCountData
        mapPlot <- clusterMap.addFreqMarkers (mapPlot, visType, clusterCountData, aggUnitData, stdMarkerCount, stdMarkerSize)
    } else if (info$mapType=="clusterPrevalence") {
        cluster <- info$cluster			
        clusterStats <- clusterSetInfo$stats
        clusterInfoText <- cluster.getClusterStatsText(clusterStats, cluster)		#; print(clusterInfoText)
            
        mapPlot <- clusterMap.addConnections (mapPlot, visType, cluster, clusterShareData, palette, stdMarkerCount, stdMarkerSize)
            
        # Dirty trick so we can show an annotation with cluster info above the legends.
        # We "plot" a couple of points, and create a ficticious alpha legend containing the cluster info text.
        bb <- baseMapInfo$gadmBB
        dummydf <- data.frame(x=c(bb$xMin,bb$xMin), y=c(bb$yMax,bb$yMax), alpha=c(0.1, 0.11))
        mapPlot <- mapPlot +
                   ggplot2::geom_point(ggplot2::aes(x=x, y=y, alpha=as.numeric(alpha), size=0.01), data=dummydf) +
                   ggplot2::scale_alpha_continuous(cluster, breaks=c(0.1, 0.11), labels=c(clusterInfoText,""), 
                                                   guide=ggplot2::guide_legend(order=1,keywidth=0, keyheight=0.01, nrow=1,
                                                   override.aes=list(shape=NA,fill=NA,size=0.01)))

    } else if (info$mapType=="clusterSharing") {
        mapPlot <- clusterMap.addShareMarkers (mapPlot, visType, clusterShareData, palette, aggUnitData, stdMarkerCount, stdMarkerSize)
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
    mapSize  <- analysis.getParam ("plot.size", params)

    clusterSetLabel <- clusterSetInfo$clusterSetName
    aggLabel <- map.getAggregationLabels(aggLevel)
    levelLabel <- getMinIdentityLabel(info$minIdentity)
    mapLabel <- mapType
    if (visType != "cluster") {
        mapLabel <- paste(mapType, visType, sep="-")
    }
    plotSubFolders <- c(paste("map", mapType, sep="-"), "plots", clusterSetLabel)
    plotFilename <- paste("map", info$sampleSetName, aggLabel, mapLabel, levelLabel, sep="-")
    if (info$mapType == "clusterPrevalence") {
        plotSubFolders <- c(plotSubFolders, levelLabel, aggLabel)
        plotFilename <- paste(plotFilename, info$cluster, sep="-")
    }
    plotFolder <- getOutFolder(ctx, info$sampleSetName, plotSubFolders)
    graphicFilename  <- paste(plotFolder, paste(plotFilename,"png",sep="."), sep="/")
    ggplot2::ggsave(plot=mapPlot, filename=graphicFilename, device="png", width=mapSize$width, height=mapSize$height, units="in", dpi=300)
}

###############################################################################
# Haplotype Frequency Markers plotting 
################################################################################
#
clusterMap.addFreqMarkers <- function (mapPlot, visType, countData, aggUnitData, stdMarkerCount, stdMarkerSize) {
    if (visType=="bar") {
        mapPlot <- clusterMap.addFreqBars (mapPlot, countData, aggUnitData, stdMarkerCount, stdMarkerSize)
    } else if (visType=="pie") {
        mapPlot <- clusterMap.addFreqPies (mapPlot, countData, aggUnitData, stdMarkerCount, stdMarkerSize)
    }
    mapPlot
}

clusterMap.addFreqPies <- function (mapPlot, countData, aggUnitData, stdMarkerCount, stdMarkerSize) {
    # Now add the pie chart markers
    if (stdMarkerCount==0) {
        mapPlot <- mapPlot +
                   ggforce::geom_arc_bar(ggplot2::aes(x0=Longitude, y0=Latitude, r0=0, r=stdMarkerSize, 
                                         fill=Haplo, amount=ClusterCount),
                             data=countData, stat="pie", inherit.aes=FALSE,
                             colour="gray25", stroke=0.5, fill="white", show.legend=FALSE)
    } else {
        mapPlot <- mapPlot +
                   ggforce::geom_arc_bar(ggplot2::aes(x0=Longitude, y0=Latitude, r0=0, r=stdMarkerSize*sqrt(SampleCount/stdMarkerCount), 
                                         fill=Haplo, amount=ClusterCount),
                             data=countData, stat="pie", inherit.aes=FALSE,
                             colour="gray25", stroke=0.5, fill="white", show.legend=FALSE)
    }
    mapPlot
}

clusterMap.addFreqBars <- function (mapPlot, countData, aggUnitData, stdMarkerCount, stdMarkerSize) {

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
               ggplot2::geom_rect(ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), data=freqBarData, inherit.aes=FALSE,
                                  colour="gray25", size=0.5, fill="white", show.legend=FALSE)
    mapPlot
}

#
# For each aggregation unit, we get a count of each unique cluster, ordered in descending count
#
clusterMap.buildCountData <- function(aggLevel, aggUnitData, dataset, params) {
    sampleMeta   <- dataset$meta
    barcodeData  <- dataset$barcodes

    # Get all aggregation unit ids
    aggUnitGids <- rownames(aggUnitData)						#; print(aggUnitGids)
    
    # Create aggregation index for each sample (the id of the aggregation unit where the sample originates)
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta)
    clusters <- apply(barcodeData,1,paste,collapse="")

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

###############################################################################
# Haplotype Sharing Markers plotting 
################################################################################
#
clusterMap.addShareMarkers <- function (mapPlot, visType, clusterShareData, clusterPalette, aggUnitData, stdMarkerCount, stdMarkerSize) {
    # Do th actual marker plotting for the twp types of markers
    if (visType=="bar") {
        mapPlot <- clusterMap.addShareBars (mapPlot, clusterShareData, clusterPalette, aggUnitData, stdMarkerCount, stdMarkerSize)
    } else if (visType=="pie") {
        mapPlot <- clusterMap.addSharePies (mapPlot, clusterShareData, clusterPalette, aggUnitData, stdMarkerCount, stdMarkerSize)
    }
    mapPlot
}

clusterMap.addSharePies <- function (mapPlot, clusterShareData, clusterPalette, aggUnitData, stdMarkerCount, stdMarkerSize) {
    # Now add the pie chart markers
    if (stdMarkerCount==0) {
        mapPlot <- mapPlot +
                   ggforce::geom_arc_bar(ggplot2::aes(x0=Longitude, y0=Latitude, r0=0, r=stdMarkerSize, 
                                         fill=Cluster, amount=ClusterCount),
                                data=clusterShareData, stat="pie", inherit.aes=FALSE,
                                colour="gray25", stroke=0.5, show.legend=FALSE) +
                   ggplot2::scale_fill_manual(values=clusterPalette)
    } else {
        mapPlot <- mapPlot +
                   ggforce::geom_arc_bar(ggplot2::aes(x0=Longitude, y0=Latitude, r0=0, r=stdMarkerSize*sqrt(SampleCount/stdMarkerCount), 
                                         fill=Cluster, amount=ClusterCount),
                                data=clusterShareData, stat="pie", inherit.aes=FALSE,
                                colour="gray25", stroke=0.5, show.legend=FALSE) +
                   ggplot2::scale_fill_manual(values=clusterPalette)
    }
    mapPlot
}


clusterMap.addShareBars <- function (mapPlot, clusterShareData, clusterPalette, aggUnitData, stdMarkerCount, stdMarkerSize) {

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
        
        # Get all the cluster for this unit, ordered by name
        unitClusterData <- clusterShareData[which(clusterShareData$UnitId==aggUnitGid),]	#; print(unitClusterData)
        unitClusters <- as.character(unitClusterData$Cluster)					#; print(unitClusters)
        unitClusterData <- unitClusterData[order(unitClusters, decreasing=TRUE),]		#; print(unitClusterData)
        unitRowCount <- nrow(unitClusterData)
        
        # For each cluster, worh out the y boundaries
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
    }												#; print(freqBarData)
    
    # Got the dataframe for drawing all the rectangles, now do the drawing and colouring
    mapPlot <- mapPlot +
               ggplot2::geom_rect(ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=Cluster), 
                                  data=freqBarData, inherit.aes=FALSE,
                                  colour="gray25", size=0.5, show.legend=FALSE) +
               ggplot2::scale_fill_manual(values=clusterPalette)
    mapPlot
}

clusterMap.buildSharedCountData <- function(aggLevel, aggUnitData, clusterData, sampleMeta, params) {
#clusterMap.buildSharedCountData <- function(aggLevel, aggUnitData, clusterData, sampleMeta, maxGroups, params) {
    # Get the subgraph membership file for this identity threshold, which is created by the clusteing code
    # (shared with the "graph" analysis task)
    memberData <- cluster.getMemberData (clusterData)

    clusterData <- clusterData[order(clusterData$Cluster),]
    #if (nrow(clusterData) > maxGroups) {
    #    clusterData <- clusterData[1:maxGroups,]
    #}
    clusterIds <- as.character(clusterData$Cluster)						#; print(clusterIds)
    memberData <- memberData[which(memberData$Cluster %in% clusterIds),]
    memberIds      <- memberData$Sample								#; print(memberIds)
    memberClusters <- memberData$Cluster							#; print(memberClusters)
   
    # Assign cluster Ids to the sample metadata
    sampleNames <- rownames(sampleMeta)
    sampleCluster <- rep("-", length(sampleNames))
    names(sampleCluster) <- sampleNames
    sampleCluster[memberIds] <- memberClusters
    sampleMeta$Cluster <- sampleCluster
    
    # Get all aggregation unit ids
    aggUnitGids <- rownames(aggUnitData)							#; print(aggUnitGids)
    
    # Create aggregation index for each sample (the id of the aggregation unit where the sample originates)
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta)

    # Get the data for all aggregation units
    countData <- NULL
    for (aIdx in 1:length(aggUnitGids)) {
        # Get the sample data to be aggregated for this unit
        aggUnitGid <- aggUnitGids[aIdx]
        unitMeta <- sampleMeta[which(aggIndex == aggUnitGid),]					#; print(nrow(unitMeta))
        clUnitMeta <- unitMeta[which(unitMeta$Cluster != "-"),]
        unitGroupCounts <- table(clUnitMeta$Cluster)
        unitGroupCounts <- unitGroupCounts[order(-unitGroupCounts)]
        unitGroupNum <- length(unitGroupCounts)+1

        unitData <- aggUnitData[aIdx,]
        df <- data.frame(matrix(unitData,ncol=ncol(aggUnitData),nrow=unitGroupNum,byrow=TRUE))
        colnames(df) <- colnames(aggUnitData)
        df$Cluster <- c(names(unitGroupCounts), "Other")
        df$ClusterCount <- c(unitGroupCounts, nrow(unitMeta)-nrow(clUnitMeta))
        df$ClusterProp <- df$ClusterCount / unitData$SampleCount
        countData <- rbind(countData, df)
    }
    countData$Cluster <- factor(countData$Cluster)
    countData$Latitude <- as.numeric(countData$Latitude)
    countData$Longitude <- as.numeric(countData$Longitude)
    countData$SampleCount <- as.numeric(countData$SampleCount)
    countData
}

###############################################################################
# Maps of connections between Haplotype Sharing sites plotting 
################################################################################
#
clusterMap.addConnections <- function (mapPlot, visType, cluster, clusterShareData, clusterPalette, stdMarkerCount, stdMarkerSize) {
    # Work out the values to display for every marker (the fractional prevalence)
    mValues <- as.numeric(clusterShareData$ClusterProp)			#; print(cluster); print(clusterShareData)
    valueLabels <- round(mValues, digits=2)
    
    # Work out size, colour and colour gradient for the markers
    pointSizes <- 16
    clusterColour <- clusterPalette[cluster]
    mMax <- max(mValues); scaleMax <- round(mMax,digits=1); scaleMax <- ifelse(scaleMax<mMax, scaleMax+0.1, scaleMax)
    mMin <- min(mValues); scaleMin <- round(mMin,digits=1); scaleMin <- ifelse(scaleMin>mMin, scaleMin-0.1, scaleMin)  #; print(paste(scaleMin, scaleMax))

    # For ech pair, work out a connectedness measure, in this case the fraction of sample pairs in which both samples carry this cluster     
    aggUnitCount <- nrow(clusterShareData)
    if (aggUnitCount > 1) {
        connectData <- NULL
        for (a1Idx in 1:(aggUnitCount-1)) {
            for (a2Idx in (a1Idx+1):aggUnitCount) {
                # Get the sample data
                lat1 <- clusterShareData$Latitude[a1Idx]
                lon1 <- clusterShareData$Longitude[a1Idx]
                sampleCount1 <- clusterShareData$SampleCount[a1Idx]
                clusterCount1 <- clusterShareData$ClusterCount[a1Idx]
                #
                lat2 <- clusterShareData$Latitude[a2Idx]
                lon2 <- clusterShareData$Longitude[a2Idx]
                sampleCount2 <- clusterShareData$SampleCount[a2Idx]
                clusterCount2 <- clusterShareData$ClusterCount[a2Idx]
                #
                pairCount <- sampleCount1 * sampleCount2
                matchCount <- clusterCount1 * clusterCount2
                prop <- matchCount / pairCount
            
                cValues <- c(lat1, lon1, lat2, lon2, prop)
                connectData <- rbind(connectData, cValues)
            }
        }
        aggUnitPairData <- data.frame(connectData)						#; print (dim(aggUnitPairData))
        nc <- ncol(aggUnitPairData)
        aggUnitPairData[,1:nc] <- sapply(aggUnitPairData[,1:nc], as.numeric)			#; print(ncol(aggUnitPairData))
        colnames(aggUnitPairData) <- c("Lat1", "Lon1", "Lat2", "Lon2", "Weight")		#; print(aggUnitPairData)
        
        hc <- grDevices::rgb2hsv(grDevices::col2rgb(clusterColour))
        weakHaploColour <- grDevices::hsv(h=hc[1,1],s=hc[2,1]*0.2,v=hc[3,1])
	
        # Now plot the connections
        mapPlot <- mapPlot +
                   ggplot2::geom_curve(ggplot2::aes(x=Lon1, y=Lat1, xend=Lon2, yend=Lat2, size=Weight, colour=Weight),	# draw edges as arcs
                               data=aggUnitPairData, curvature=0.2, alpha=0.75) +
                   ggplot2::scale_size_continuous(guide=FALSE, range=c(0.25, 4)) +          			# scale for edge widths
                   ggplot2::scale_colour_gradientn(name="Sharing", colours=c(weakHaploColour,clusterColour))
    }
     
    # Plot the Aggregation Unit markers
    mapPlot <- mapPlot +
	       ggplot2::geom_point(ggplot2::aes(x=Longitude, y=Latitude, fill=ClusterProp), 
	                           data=clusterShareData, size=pointSizes, shape=21, stroke=2) +
	       ggplot2::scale_fill_gradientn(name="Prevalence", limits=c(scaleMin,scaleMax), colours=c("white",clusterColour), values=c(0,1)) +
               ggplot2::geom_text(ggplot2::aes(x=Longitude, y=Latitude), 
                                  data=clusterShareData, label=valueLabels, hjust=0.5, vjust=0.5, size=4.5, fontface="bold")
    mapPlot
}

