ADM_DIV_COLUMNS <- c("Country", "AdmDiv1", "AdmDiv2")
ADM_DIV_LABELS  <- c("Country", "Province", "District")
GID_COLUMNS     <- c("Country", "AdmDiv1_GID", "AdmDiv2_GID")

###############################################################################
# Common Routines for Map Generation.
################################################################################
#
map.execute <- function(userCtx, sampleSetName, mapType, params) {		#;print(mapType)	;print(params)
    #
    # For pie chart maps, we will need to create "category palettes" to ensure consistent colour mappings across time periods
    #
    if (mapType=="alleleProp") {
        sampleSet <- userCtx$sampleSets[[sampleSetName]] 
        sampleSet$categoryPalettes <- new.env()
    }
    #
    # Get the set of specs for the maps to be produced
    #
    mapSpecs <- map.createMapSpecs (userCtx, sampleSetName, mapType, params)	#;print(mapSpecs)
    #
    # Now create the map plots one by one
    #
    for (mIdx in 1:length(mapSpecs)) {
        mapSpec <- mapSpecs[[mIdx]]
        #
        # Generate the ggplot2 plot object
        #
        mapPlot <- map.generateMapPlot (mapSpec, params)
        if (is.null(mapPlot)) {
            msg <- paste("Skipping map:", mapSpec$master$type, mapSpec$measure)
            if (!is.null(mapSpec$interval)) {
                 msg <- paste(msg, mapSpec$interval$name)
            }
            print(msg)
            next
        }
        #
        # Add the "look and feel" elements to the map plot
        #
        mapPlot <- mapPlot + mapSpec$master$theme				#; plot(mapPlot)
        #
        # Deal with the legend 
        #
        mapPlot <- map.processLegend (mapPlot, mapSpec, params)			#; plot(mapPlot)
        #
        # Save to file. The size in inches is given in the config.
        #
        ggplot2::ggsave(plot=mapPlot, filename=mapSpec$plotFile, device=params$plot.file.format, 
                        width=params$plot.width, height=params$plot.height, units=params$plot.units, dpi=params$plot.dpi)
    }
}
#
map.processLegend <- function(mapPlot, mapSpec, params) {			#; print(mapSpec$master$type)
    if (mapSpec$master$type %in% c("drug", "mutation", "alleleProp", "diversity", "sampleCount", "location")) {
        legPos <- params$plot.legend.pos					#; print(legPos)
        legDir <- params$plot.legend.dir					#; print(legDir)
        if (legDir == "horizontal") {
            mapPlot <- mapPlot + ggplot2::theme(legend.position="bottom")
        } else {
            mapPlot <- mapPlot + ggplot2::theme(legend.position="right")
        }
        if (legPos == "separate") {
            #
            # Separate the legend from the map, and set the map without the legend to be the plot to be wirtten to file
            #
            legendPlot <- cowplot::get_legend(mapPlot)				#; gtable::gtable_show_layout(legendPlot)
            mapPlot <- mapPlot + ggplot2::theme(legend.position="none")
            #
            # Find out the size of the legend (note, we write it out so we have a device for the size estimation)
            #
            legendFile <- paste0(mapSpec$plotFile, ".legend.png")
            grDevices::png(filename=legendFile, width=params$plot.width, height=params$plot.height, units=params$plot.units, res=params$plot.dpi)
            legW <- grid::convertWidth(grid::widthDetails(legendPlot), params$plot.units, valueOnly=TRUE)	#; print(paste(legW, params$plot.units))
            legH <- grid::convertHeight(grid::heightDetails(legendPlot), params$plot.units, valueOnly=TRUE)	#; print(paste(legH, params$plot.units))
            grDevices::dev.off()
            #
            # Write out the legend to a file of the exact size
            #
            legendPlot <- ggpubr::as_ggplot(legendPlot)
            #legendPlot <- legendPlot + 
            #              ggplot2::theme_minimal() + 
            #              #ggplot2::theme(plot.margin=grid::unit(c(5.5, 20, 5.5, 20), "pt"))
            #              ggplot2::theme(legend.margin=ggplot2::margin(l=20, r=20, unit='pt'))
            ggplot2::ggsave(plot=legendPlot, filename=legendFile, device="png", bg="white",
                            width=legW*1.1, height=legH*1.1, units=params$plot.units, dpi=params$plot.dpi)
        }
    }
    mapPlot
}
#
###############################################################################
# Creation of a list of specifications for generating the maps needed
################################################################################
#
map.createMapSpecs <- function(userCtx, sampleSetName, mapType, params) {		#;print(mapType)
    #
    # Initialize the master spec on which to base the rest
    #
    mapMaster <- map.createMapMaster (userCtx, sampleSetName, mapType, params)
    #
    # Create lists of maps to be produced, adding the master as the first element
    #
    mapSpecs <- list()
    #
    # Create one map object for each time interval/dataset/aggregation/measure combination, and add it to the maps object
    #
    datasetNames <- mapMaster$datasetNames
    for (dnIdx in 1:length(datasetNames)) {					#; print(paste(length(datasetNames),"datasetNames:",datasetNames))
        datasetName <- datasetNames[dnIdx]					#; print(paste("datasetName:",datasetName))
        intervals <- mapMaster$intervals
        for (iIdx in 1:length(intervals)) {					#; print(paste(length(intervals),"intervals:",intervals))
            interval <- intervals[[iIdx]]					#; print(paste("interval:",interval))
            #
            # Trim the dataset according to time intervel, and skip if there are no samples
            #
            mapCtx <- context.trimContextByTimeInterval (mapMaster$sampleSet$ctx, interval)
            if (is.null(mapCtx)) {
                print(paste("No samples found- skipping interval", interval$name))
                next
            }
            aggLevels <- mapMaster$aggregations
            for (aIdx in 1:length(aggLevels)) {
                aggLevel <- as.integer(aggLevels[aIdx])				#; print(paste("aggLevel:",aggLevel))
                aggUnitData <- map.getAggregationUnitData (mapCtx, datasetName, aggLevel, sampleSetName, 
		                   mapType, params, mapMaster$dataFolder)	#; print(aggUnitData)
                measures    <- mapMaster$measures				#; print(paste(length(measures),"measures:",measures))
                #
                aggUnitPairData <- NULL
                if (map.isMarkerMap(mapType)) {
                    #
                    # Get aggregation units and values for their measures, except for for sample counts which are already in the dataframe
                    #
                    if (!(mapType %in% c("sampleCount","location"))) {
		        aggUnitData <- markerMap.estimateMeasures (mapCtx, datasetName, aggLevel, aggUnitData, sampleSetName, 
		                                                   mapType, measures, params, mapMaster$dataFolder)	#; print(aggUnitData)
		    }
                } 
                if (mapType=="connect") {
		    aggUnitPairData <- connectMap.estimateMeasures (mapCtx, datasetName, sampleSetName, aggLevel, aggUnitData, 
		    	                                            mapType, measures, params, mapMaster$dataFolder)	#; print(aggUnitPairData)
		}
                #
                #measures <- measures[1]
                for (mIdx in 1:length(measures)) {
                    measure <- measures[mIdx]			#; print(measure)
                    #
                    mapSpec <- list()
                    mapSpec$master      <- mapMaster
                    mapSpec$mapCtx      <- mapCtx
                    mapSpec$datasetName <- datasetName
                    mapSpec$aggregation <- aggLevel
                    mapSpec$interval    <- interval
                    mapSpec$measure     <- measure
                    mapSpec$aggUnitData <- aggUnitData
                    mapSpec$aggUnitPairData <- aggUnitPairData
                    mapSpec$plotFile    <- map.getMapFilepath(mapSpec, params)
                    #
                    mapSpecs[[length(mapSpecs)+1]] <- mapSpec
                }
            }
        }
    }
    mapSpecs
}
#
###############################################################################
# Master object for generating the maps
################################################################################
#
map.createMapMaster <- function (userCtx, sampleSetName, mapType, params) {		#; print(names(params))
    #
    # Make a template config to be used for each time slice interval
    #
    sampleSet  <- userCtx$sampleSets[[sampleSetName]]
    #
    mapMaster <- new.env()
    mapMaster$type        <- mapType
    mapMaster$sampleSet   <- sampleSet
    mapMaster$userCtx     <- userCtx
    mapMaster$params      <- params
    #
    #
    #
    if (mapType %in% c("drug", "mutation", "alleleProp", "location")) {
        mapMaster$datasetNames <- "unfiltered"
    } else if (mapType  %in% c("diversity", "connect", "barcodeFrequency", "clusterSharing", "clusterPrevalence")) {
        mapMaster$datasetNames <- "imputed"
    } else if (mapType == "sampleCount") {
        mapMaster$datasetNames <- c("unfiltered","filtered")
    } else {
        stop(paste("Invalid map type:", mapType))
    }    
    #
    #
    #
    aggregations <- param.getParam ("aggregation.levels", params)	#;print(aggregation)
    mapMaster$aggregations <- aggregations
    #
    # Get the output folders and file name elements
    #
    mapMaster$dataFolder <- getOutFolder(userCtx$config, sampleSet$name, c(paste("map", mapType, sep="-"), "data"))
    mapMaster$plotFolder <- getOutFolder(userCtx$config, sampleSet$name, c(paste("map", mapType, sep="-"), "plots"))	#; print(mapMaster$plotFolder)
    #
    # Get the output image formatting parameters
    #
    mapMaster$theme <- map.getTheme(params)				#; plot(mapMaster$theme)
    #
    # Handle the "ALL" values for markers
    # Check the measures specified are valid; and if they are, ensure we have colour palettes for these measures, creating them if necessary.
    # For a given measure, all plots for this sample set use the same palette, otherwise the viewer will be confused when looking at multiple maps.
    #
    if (mapType %in% c("drug", "mutation", "diversity")) {
        measures <- markerMap.resolveMeasures (userCtx, mapType, params)

    } else if (mapType=="alleleProp") {
        measures <- pieMap.resolveMeasures (userCtx, sampleSetName, params)
        
    } else if (mapType=="connect") {
        measures <- connectMap.resolveMeasures(params)

    } else if (mapType %in% c("location", "sampleCount")) {
        measures <- param.getParam ("analysis.measures", params)
        
    } else if (mapType %in% c("clusterPrevalence","clusterSharing")) {
        # Add one measure for each similarity threshold (e.g. "clusterPrevalence-ge0.80" or "clusterSharing-ge0.80")
        clusterSetName  <- param.getParam ("cluster.clusterSet.name", params)			#; print(clusterSetName)
        mapMaster$clusterSets <- cluster.getClustersSetFromContext (userCtx, sampleSetName, clusterSetName)
        mapMaster$clusterSetPalettes <- clusterMap.getClusterSetsPalettes (userCtx, mapMaster$clusterSets)
        if (mapType=="clusterPrevalence") {
            measures <- clusterMap.resolveMeasures (mapMaster$clusterSets)
        } else if (mapType=="clusterSharing") {
            measures <- clusterShareMap.resolveMeasures (mapMaster$clusterSets, params)
        }

    } else if (mapType=="barcodeFrequency") {
        measures <- barcodeFreqMap.resolveMeasures (mapMaster$clusterSets, params)		#; print(names(params))

    }
    mapMaster$measures <- measures					#; print(paste(length(measures),"measures:",measures))
    #
    # If no time interval is specified, assign a default one (all samples)
    #
    intervals <- NULL
    if (mapType %in% c("drug", "mutation", "alleleProp", "diversity", "sampleCount")) {
        intervals <- param.getParam ("analysis.timeIntervals", params)
    }
    if (is.null(intervals)) {
        intervals <- list(getDefaultTimeInterval())
    }
    mapMaster$intervals <- intervals
    #
    # Get the base map to use as a background for all plots
    #
    baseMapInfo <- baseMap.buildBaseMap (mapMaster)
    mapMaster$baseMapInfo <- baseMapInfo
    mapMaster
}
#
#
#
map.getMapFilepath <- function (map, params) {
    mapMaster <- map$master
    mapType <- mapMaster$type
    #
    # Compose a plot filename
    #
    mapLabel <- map$measure
    plotFolder <- mapMaster$plotFolder
    if (mapType %in% c("clusterPrevalence","clusterSharing")) {
        if (mapType == "clusterPrevalence") {
            mParts <- clusterMap.parseMeasure (map$measure)				#;print(map$measure)
        } else if (mapType == "clusterSharing") {
            mParts <- clusterShareMap.parseMeasure (map$measure)			#;print(map$measure)
        }
        minIdentityLabel <- getMinIdentityLabel(mParts$minIdentity)			#;print(minIdentityLabel)
        clusterSet <- mapMaster$clusterSets[[minIdentityLabel]]				#;print(names(clusterSet))
        clusterSetName <- clusterSet$clusterSetName					#;print(clusterSetName)
        plotFolder <- paste(plotFolder, clusterSetName, sep="/")			#;print(plotFolder)
        if (mapType %in% c("clusterPrevalence")) {
            plotFolder <- paste(plotFolder, minIdentityLabel, sep="/")			#;print(plotFolder)
        }
    }
    if (mapType %in% c("barcodeFrequency","clusterPrevalence","clusterSharing")) {
        mapLabel <- gsub("/", "-", mapLabel)
    }
    #
    aggLabel <- map.getAggregationLabels(map$aggregation)
    if (mapType=="sampleCount") {
        aggLabel <- paste(aggLabel, map$datasetName, sep="-")
    }
    #
    plotFilename  <- paste("map", mapMaster$sampleSet$name, aggLabel, mapLabel, sep="-")
    if (!is.null(map$interval$name)) {
        plotFilename  <- paste(plotFilename, map$interval$name, sep="-")
    }
    
    fileExt <- param.getParam ("plot.file.format", params)				#; print(fileExt)
    plotFilename <- paste0(plotFilename,".", fileExt)
    plotFile <- paste(plotFolder, plotFilename, sep="/")				#;print(plotFile)
    plotFile
}    
#
#
#
map.getTheme <- function(params) {
    th <- ggplot2::theme_grey(base_size=params$themeBaseSize)
    th <- th +
          ggplot2::theme(panel.background = ggplot2::element_rect(fill=map.colour.sea, colour=NA),
                         panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),
                         plot.background  = ggplot2::element_rect(fill="white", colour=NA),
                         plot.title       = ggplot2::element_text(face="bold",size=ggplot2::rel(params$plotTitleSize), hjust=0.5),
                         axis.title       = ggplot2::element_text(face="bold",size=ggplot2::rel(params$axisTitleSize)),
                         axis.title.y     = ggplot2::element_text(angle=90, vjust=2),
                         axis.title.x     = ggplot2::element_text(vjust=-0.2),
                         legend.key.size  = ggplot2::unit(params$legendKeySize, "lines")
                         )
    th
}
#
###############################################################################
# Creation of a list of specifications for generating the maps needed
################################################################################
#
map.generateMapPlot <- function(mapSpec, params) {	#; print(names(mapSpec))
    mapType <- mapSpec$master$type 
    #
    # Create the map plot
    #
    mapPlot <- NULL
    if (map.isMarkerMap(mapType)) {
        mapPlot <- markerMap.executeMap (mapSpec)
    } else if (mapType=="alleleProp") {
        mapPlot <- pieMap.executeMap (mapSpec)
    } else if (mapType=="connect") {
        mapPlot <- connectMap.executeMap (mapSpec)
    } else if (mapType=="clusterPrevalence") {
        mapPlot <- clusterMap.executeMap (mapSpec)
    } else if (mapType=="clusterSharing") {
        mapPlot <- clusterShareMap.executeMap (mapSpec)
    } else if (mapType=="barcodeFrequency") {
        mapPlot <- barcodeFreqMap.executeMap (mapSpec)
    } else {
        stop(paste("Invalid map type:", mapType))
    }							#; plot (mapPlot)
    mapPlot
}
#
map.isMarkerMap <- function(mapType) {
    mapType %in% c("drug", "mutation", "diversity", "sampleCount", "location")
}
#
###############################################################################
# Aggregation of site data
################################################################################
#
map.getAggregationColumns <- function(aggLevels=c()) {
    if (length(aggLevels) == 0) {
        return (ADM_DIV_COLUMNS)
    }
    colIdx <- aggLevels + 1
    ADM_DIV_COLUMNS[colIdx]
}
#
map.getAggregationLabels <- function(aggLevels=c()) {
    if (length(aggLevels) == 0) {
        return (ADM_DIV_LABELS)
    }
    labelIdx <- aggLevels + 1
    ADM_DIV_LABELS[labelIdx]
}
#
map.getAggregationLevelsFromLabels <- function(aggLabels) {
    result <- c()
    for (i in 1:length(aggLabels)) {
        aggLabel <- aggLabels[i]
        idx <- which(ADM_DIV_LABELS == aggLabel)
        if (length(idx) == 0) {
            stop (paste("Invalid aggregation level specified:",aggLabel))
        }
        result <- c(result, (idx-1))
    }
    result <- as.integer(result)
    result 
}
#
map.getAdmDivNames <- function(gids) {		#; print(gids)
    geo <- map.getGeoTables()
    admDivs <- geo$admDivs			#; print(admDivs)
    admDivs <- admDivs[gids,]
    admDivNames <- admDivs$Name
    admDivNames
}
#
map.getAdmDivNamesFromMeta <- function (admDivLevel, sampleMeta) {
    nameCol <- ADM_DIV_COLUMNS[admDivLevel+1]				#; print(nameCol)
    gidCol  <- GID_COLUMNS[admDivLevel+1]				#; print(gidCol)
    df <- data.frame(name=as.character(sampleMeta[,nameCol]),
                     id=as.character(sampleMeta[,gidCol]))
    df <- df[!duplicated(df),]
    admDivNames <- df$name
    names(admDivNames) <- df$id
    admDivNames
}
#
map.getAggregationUnitData <- function(plotCtx, datasetName, aggLevel, analysisName, mapType, params, dataFolder) {

    # Trim all data to discard samples that have incomplete geographical data
    validSamples <- map.getAggregableSamples (plotCtx, datasetName, aggLevel)
    
    dataset <- plotCtx[[datasetName]]
    sampleMeta   <- dataset$meta[validSamples,]

    adminLevelCols  <- map.getAggregationColumns()

    # Create aggregation unit id; this is the GID, unique for each aggregation unit 
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta)

    # Get all aggregation units, in order, and keep only those that have enough samples
    aggregateCountMin <- param.getParam ("map.aggregateCountMin", params)
    aggUnitCounts <- table(aggIndex)
    aggUnitGids <- names(aggUnitCounts[aggUnitCounts >= aggregateCountMin])	#; print(aggUnitGids)
    aggUnitGids <- aggUnitGids[order(aggUnitGids)]				#; print(aggUnitGids)
    
    # Get the data for all aggregation units
    aggUnitCnames <- c("UnitId", "Country", "AdmDivName", "AdmDiv1_GID", "AdmDiv2_GID", "Latitude", "Longitude", "SampleCount")
    aggUnitData <- matrix(nrow=0, ncol=length(aggUnitCnames))

    geo <- map.getGeoTables()
    admDivs <- geo$admDivs
    admDivNames <- map.getAdmDivNamesFromMeta (aggLevel, sampleMeta)
    for (aIdx in 1:length(aggUnitGids)) {
        #
        # Get the sample data to be aggregated for this unit
        #
        aggUnitGid <- aggUnitGids[aIdx]
        aggUnitName <- admDivNames[aggUnitGid]
        aggSamplesMeta <- sampleMeta[which(aggIndex == aggUnitGid),]		#; print(nrow(aggSamplesMeta))
        aggSamples <- rownames(aggSamplesMeta)
        #
        # Get the admin division values from the first sample of this unit
        #
        #admDiv <- admDivs[which(admDivs$GID==aggUnitGid),]			#; print(admDiv)
        admDiv <- admDivs[aggUnitGid,]						#; print(admDiv)
        if (aggLevel == 1) {
            admDiv1_GID <- aggUnitGid
            admDiv2_GID <- "-"
        } else if (aggLevel == 2) {
            admDiv1_GID <- admDiv$Parent
            admDiv2_GID <- aggUnitGid
        } else {
            stop (paste("Unsupported aggregation level."))
        }									#; print(paste(admDiv1_GID, admDiv2_GID))
        cValues <- c(aggUnitGid, admDiv$Country, aggUnitName, admDiv1_GID, admDiv2_GID, 
                     admDiv$Latitude, admDiv$Longitude, nrow(aggSamplesMeta))	#; print(cValues)
        aggUnitData <- rbind(aggUnitData, cValues)
    }
    aggUnitData <- data.frame(aggUnitData)
    rownames(aggUnitData) <- aggUnitGids
    colnames(aggUnitData) <- aggUnitCnames
    aggUnitData$Latitude    <- as.numeric(as.character(aggUnitData$Latitude))
    aggUnitData$Longitude   <- as.numeric(as.character(aggUnitData$Longitude))
    aggUnitData$SampleCount <- as.integer(as.character(aggUnitData$SampleCount))

    # Write out the aggregation unit data to file
    aggDataFilename  <- paste(dataFolder, "/AggregationUnits-", analysisName, "-", aggLevel, ".tab", sep="")
    utils::write.table(aggUnitData, file=aggDataFilename, sep="\t", quote=FALSE, row.names=FALSE)

    aggUnitData
}
#
map.getAggregableSamples <- function(plotCtx, datasetName, aggLevel) {
    dataset <- plotCtx[[datasetName]]
    sampleMeta <- dataset$meta
    for (aIdx in 0:aggLevel) {
        cName  <- map.getAggregationColumns(aIdx);
        missingIdx <- which (sampleMeta[,cName] == "-")
        if (length(missingIdx) > 0) {
            sampleMeta <- sampleMeta[-missingIdx,]
        }
    }
    rownames(sampleMeta)
}
#
map.getAggregationUnitIds <- function(aggLevel, sampleMeta) {	#; print(sampleMeta)
    gidCol <- GID_COLUMNS[aggLevel+1]
    aggUnitId <- sampleMeta[,gidCol]
    aggUnitId
}
#
###############################################################################
# Geo Data (e.g. GADM names)
################################################################################
#
map.getGeoTables <- function () {
    grcMalariaGeodata::getGeoTables()	# From the .rda file
}
map.iso2ToIso3 <- function (iso2Countries) {
    geo <- map.getGeoTables()
    countryTable <- geo$country.codes
    iso3 <- countryTable[iso2Countries,"IsoCode3"]
    iso3
}
map.getCountryData <- function (iso2Country) {
    geo <- map.getGeoTables()
    countryTable <- geo$country.codes
    countryData <- countryTable[iso2Country,]
    countryData
}
#
###############################################################################
# Labelling of aggregated sites
################################################################################
#
map.addAggregationUnitNameLayer <- function (mapPlot, aggUnitData, baseMapInfo, params, markerSize=NULL) {
    #
    cnt <- nrow(aggUnitData)						#; print(aggUnitData) 
    lon <- as.numeric(aggUnitData$Longitude)				#; print(lon)
    lat <- as.numeric(aggUnitData$Latitude)				#; print(lat) 
    label <- as.character(aggUnitData$AdmDivName)			#; print(label) 
    #
    # The following segment seems to be unnecessary, so I commented it for the time being
    #
    #bbox <- baseMapInfo$anBB						#; print(bbox)
    #xNudgeSize <- (bbox$xMax - bbox$xMin)/15				#; print(xNudgeSize)
    #yNudgeSize <- (bbox$yMax - bbox$yMin)/15				#; print(yNudgeSize)
    #nudge_x <- nudge_y <- vector(mode="integer", length=cnt)
    #just <- rep(0.5,cnt)
    #qlon <- stats::quantile(lon)						#; print(qlon) 
    #nudge_x[which(lon<=qlon[2])] <- -xNudgeSize; just[which(lon<=qlon[2])] <- 1
    #nudge_x[which(lon<=qlon[4])] <-  xNudgeSize; just[which(lon<=qlon[4])] <- 0
    #qlat <- stats::quantile(lat)
    #nudge_y[which(lat<qlat[3])]  <- -yNudgeSize
    #nudge_y[which(lat>=qlat[3])] <-  yNudgeSize
    #y[which(x!=0)] <- 0
    #lp <- data.frame(lat=lat, lon=lon, label=label, nudge_x=c(nudge_x), nudge_y=c(nudge_y), just=just, markerSize=markerSize)	#; print(lp)
    #
    #
    #
    if (is.null(markerSize)) {
        markerSize <- rep (10, cnt)
    }
    lp <- data.frame(lat=lat, lon=lon, label=label, markerSize=markerSize)	#; print(lp)
    mapPlot <- mapPlot + 
        #ggplot2::geom_point(ggplot2::aes(x=lon, y=lat), data=lp, colour="red") +
        ggrepel::geom_label_repel(ggplot2::aes(x=lon, y=lat, label=label, point.size=markerSize), 
                                  data=lp,
                                  size=params$map.admNameLabelFontSize, fontface="bold", color="grey50",
                                  #hjust=lp$just, vjust=0.5, nudge_x=lp$nudge_x, nudge_y=lp$nudge_y, force=4,
                                  show.legend=FALSE,
                                  segment.size=params$admNameLineWidth, label.size=params$admNameLineWidth,
                                  label.padding=grid::unit(params$admNameLabelPadding, "lines"))
    mapPlot
}
