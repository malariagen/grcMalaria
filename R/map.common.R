
ADM_DIV_COLUMNS <- c("Country", "AdmDiv1", "AdmDiv2")
ADM_DIV_LABELS  <- c("Country", "Province", "District")
GID_COLUMNS     <- c("Country", "AdmDiv1_GID", "AdmDiv2_GID")

###############################################################################
# Common Routines for Map Generation.
################################################################################
#
map.execute <- function(userCtx, sampleSetName, mapType, aggregation, measures, params) {

    if (mapType %in% c("drug", "mutation")) {
        markerMap.execute (userCtx, "unfiltered", sampleSetName, mapType, aggregation, measures, params)
        
    } else if (mapType == "diversity") {
        markerMap.execute (userCtx, "imputed",    sampleSetName, mapType, aggregation, measures, params)
        
    } else if (mapType == "sampleCount") {
        markerMap.execute (userCtx, "unfiltered", sampleSetName, mapType, aggregation, measures, params)
        markerMap.execute (userCtx, "filtered",   sampleSetName, mapType, aggregation, measures,   params)
        
    } else if (mapType == "connect") {
        connectMap.execute (userCtx, "imputed",   sampleSetName, mapType, aggregation, measures, params)
        
    } else if (mapType == "barcodeFrequency") {
        clusterMap.execute (userCtx, "imputed",   sampleSetName, mapType, aggregation, measures, params)

    } else if (mapType == "clusterSharing") {
        clusterMap.execute (userCtx, "imputed",   sampleSetName, mapType, aggregation, measures, params)

    } else if (mapType == "clusterPrevalence") {
        clusterMap.execute (userCtx, "imputed",   sampleSetName, mapType, aggregation, measures, params)

    } else {
        stop(paste("Invalid map type:", mapType))
    }
}

#
###############################################################################
# Creation of the physical/political background map
################################################################################
#
map.buildBaseMap <- function(ctx, datasetName, analysisName, sampleMeta, dataFolder, params) {

    # Get relevant column names
    cCountry  <- map.getAggregationColumns(0)

    # Now read the countries and provinces so we can get the contours from GADM
    unitList <- list()
    countryValues <- sampleMeta[,cCountry]
    countries <- unique(countryValues)
    for (cIdx in 1:length(countries)) {
        country <- countries[cIdx]					#; print(country)
        countryData <- map.getCountryData (country)			#; print(countryData)
        cMeta <- sampleMeta[which(countryValues == country),]
        adm1GIDValues <- as.character(cMeta$AdmDiv1_GID)
        adm1GIDs <- unique(adm1GIDValues)				#; print(provIds)
        unitList[[country]] <- list(name=countryData$name, iso2=countryData$iso2, iso3=countryData$iso3, adm1GIDs=adm1GIDs) 
    }
    
    # Read the country borders for the countries involved
    gadmFolder <- getCacheFolder(ctx, c("map", "gadm"))
    gadmFolder <- paste(gadmFolder, "/", sep="")
    cIso3 <- map.iso2ToIso3 (countries)
    gadm0 <- GADMTools::gadm_sp_loadCountries(cIso3, level=0, basefile=gadmFolder)
    
    # Read the province borders for the countries involved
    gadm1Spdf <- NULL
    gadmBB <- list(xMin=1000, xMax=-1000, yMin=1000, yMax=-1000)

    for (cIdx in 1:length(countries)) {
        country <- countries[cIdx]								#; print(country)
        cl <- unitList[[country]]								#; print(cl)
        cGadm1 <- GADMTools::gadm_sp_loadCountries(cl$iso3, level=1, basefile=gadmFolder)

        # Select the provinces we need
        adm1GIDs <- cl$adm1GIDs									#; print(adm1GIDs)
        cGadm1 <- GADMTools::gadm_subset(cGadm1, level=1, regions=adm1GIDs, usevar="GID_1")	#; print(cGadm1)

        # Append the data to that of other countries
        cl$gadmAdm1Data <- cGadm1
        if (is.null(gadm1Spdf)) {
            gadm1Spdf <- cGadm1$spdf
        } else {
            gadm1Spdf <- rbind(gadm1Spdf, cGadm1$spdf)
        }
        # Get the bounding box and merge with the other countries
        bb <- GADMTools::gadm_getBbox (cGadm1)
        gadmBB <- list(xMin=min(gadmBB$xMin,bb[1]), xMax=max(gadmBB$xMax,bb[3]), 
                       yMin=min(gadmBB$yMin,bb[2]), yMax=max(gadmBB$yMax,bb[4]))		#; print(gadmBB)
    }
    
    # Adjust the bounding box to give some margin
    xMar <- (gadmBB$xMax-gadmBB$xMin)/20
    yMar <- (gadmBB$yMax-gadmBB$yMin)/20
    gadmBB <- list(xMin=(gadmBB$xMin-xMar), xMax=(gadmBB$xMax+xMar), yMin=(gadmBB$yMin-yMar), yMax=(gadmBB$yMax+yMar))
    gadmBB$tl <- c(gadmBB$yMax,gadmBB$xMin);    gadmBB$br <- c(gadmBB$yMin,gadmBB$xMax)
    gadmBB$bl <- c(gadmBB$yMin,gadmBB$xMin);    gadmBB$tr <- c(gadmBB$yMax,gadmBB$xMax)
    
    # Crop the country boundaries
    gadm0 <- GADMTools::gadm_crop(gadm0, xmin=gadmBB$xMin, ymin=gadmBB$yMin, xmax=gadmBB$xMax, ymax=gadmBB$yMax)

    # Prepare the GADM polygons for plotting
    gadm0_df <- ggplot2::fortify(gadm0$spdf)    	#; print(colnames(gadm0$spdf)); print(colnames(gadm0_df))
    gadm1_df <- ggplot2::fortify(gadm1Spdf)		#; print(colnames(gadm1Spdf))print(colnames(gadm1_df))

    # Get the background map, and adjust coordinates
    bgMap <- OpenStreetMap::openmap(gadmBB$tl, gadmBB$br, zoom=NULL, type=c("nps"), minNumTiles=4, mergeTiles=TRUE)
    ## OSM CRS :: "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"
    bgMap <- OpenStreetMap::openproj(bgMap, projection="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

    # Construct a base plot for completing subsequent maps
    baseMapPlot <- OpenStreetMap::autoplot.OpenStreetMap(bgMap)  +
            ggplot2::annotate("rect", xmin=gadmBB$xMin, ymin=gadmBB$yMin, xmax=gadmBB$xMax, ymax=gadmBB$yMa, fill="white", alpha=0.3) +
    	    #ggplot2::labs(title=analysisName, subtitle="", x="Longitude", y="Latitude")+
    	    ggplot2::labs(x="Longitude", y="Latitude")+
    	    ggplot2::geom_polygon(data=gadm1_df, ggplot2::aes(x=long, y=lat, group=group), bg=NA, col="black", size=1) +
    	    ggplot2::geom_polygon(data=gadm0_df, ggplot2::aes(x=long, y=lat, group=group), bg=NA, col="black", size=1.5)

    # Return all the elements
    list(baseMap=baseMapPlot, bgMap=bgMap, gadmBB=gadmBB, gadm0_df=gadm0_df, gadm1_df=gadm1_df)
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
map.getAdmDivNames <- function(gids) {
    geo <- map.getGeoTables()
    admDivs <- geo$admDivs
    admDivs <- admDivs[gids,]
    admDivNames <- admDivs$AdmDivName
    admDivNames
}
#
map.getAggregationUnitData <- function(ctx, datasetName, aggLevel, analysisName, mapType, params, dataFolder) {

    # Trim all data to discard samples that have incomplete geographical data
    validSamples <- map.getAggregableSamples (ctx, datasetName, aggLevel)
    
    dataset <- ctx[[datasetName]]
    sampleMeta   <- dataset$meta[validSamples,]
    barcodeData  <- dataset$barcodes[validSamples,]
    distData     <- dataset$distance[validSamples,validSamples]

    aggLevelIdx <- aggLevel + 1
    adminLevelCols  <- map.getAggregationColumns()

    # Create aggregation unit id; this is the GID, unique for each aggregation unit 
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta)

    # Get all aggregation units, in order, and keep only those that have enough samples
    aggregateCountMin <- analysis.getParam ("map.aggregateCountMin", params)
    aggUnitCounts <- table(aggIndex)
    aggUnitGids <- names(aggUnitCounts[aggUnitCounts >= aggregateCountMin])	#; print(aggUnitGids)
    aggUnitGids <- aggUnitGids[order(aggUnitGids)]				#; print(aggUnitGids)
    
    # Get the data for all aggregation units
    aggUnitCnames <- c("UnitId", "Country", "AdmDivName", "AdmDiv1_GID", "AdmDiv2_GID", "Latitude", "Longitude", "SampleCount")
    aggUnitData <- matrix(nrow=0, ncol=length(aggUnitCnames))
    geo <- map.getGeoTables()
    admDivs <- geo$admDivs
    for (aIdx in 1:length(aggUnitGids)) {
        #
        # Get the sample data to be aggregated for this unit
        #
        aggUnitGid <- aggUnitGids[aIdx]
        aggSamplesMeta <- sampleMeta[which(aggIndex == aggUnitGid),]		#; print(nrow(aggSamplesMeta))
        aggSamples <- rownames(aggSamplesMeta)
        aggBarcodes <- barcodeData[aggSamples,]
        aggDist <- distData[aggSamples,aggSamples]
        #
        # Get the admin division values from the first sample of this unit
        #
        admDiv <- admDivs[which(admDivs$GID==aggUnitGid),]			#; print(admDiv)
        cValues <- c(aggUnitGid, admDiv$Country, admDiv$AdmDivName, admDiv$AdmDiv1_GID, admDiv$AdmDiv2_GID, 
                     admDiv$Latitude, admDiv$Longitude, nrow(aggSamplesMeta))	#; print(cValues)
        aggUnitData <- rbind(aggUnitData, cValues)
    }
    aggUnitData <- data.frame(aggUnitData)
    aggUnitData[, 6:8] <- sapply(aggUnitData[, 6:8], as.numeric)
    rownames(aggUnitData) <- aggUnitGids
    colnames(aggUnitData) <- aggUnitCnames

    # Write out the aggregation unit data to file
    aggDataFilename  <- paste(dataFolder, "/AggregationUnits-", analysisName, "-", aggLevel, ".tab", sep="")
    utils::write.table(aggUnitData, file=aggDataFilename, sep="\t", quote=FALSE, row.names=FALSE)

    aggUnitData
}
#
map.getAggregableSamples <- function(ctx, datasetName, aggLevel) {
    dataset <- ctx[[datasetName]]
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
map.getAggregationUnitIds <- function(aggLevel, sampleMeta) {
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
    map.geoTables	# From the .rda file
}
map.iso2ToIso3 <- function (iso2Countries) {
    geo <- map.getGeoTables()
    countryTable <- geo$countries
    iso3 <- countryTable[iso2Countries,"iso3"]
    iso3
}
map.getCountryData <- function (iso2Country) {
    geo <- map.getGeoTables()
    countryTable <- geo$countries
    countryData <- countryTable[iso2Country,]
    countryData
}
#
###############################################################################
# Labelling of aggregated sites
################################################################################
#
map.computeLabelParams <- function (aggUnitData, baseMapInfo) {
    bbox <- baseMapInfo$gadmBB
    xNudge <- (bbox$xMax - bbox$xMin)/15
    yNudge <- (bbox$yMax - bbox$yMin)/15
    
    cnt <- nrow(aggUnitData)
    lon <- aggUnitData$Longitude;
    lat <- aggUnitData$Latitude
    lab <- aggUnitData$AdmDivName
    
    x <- y <- vector(mode="integer", length=cnt)
    just <- rep(0.5,cnt)
    qlon <- stats::quantile(lon)		#; print (qlat)
    qlat <- stats::quantile(lat)		#; print (qlon)
    x[which(lon<=qlon[2])] <- -xNudge;	x[which(lon>=qlon[4])] <- xNudge
    just[which(lon<=qlon[2])] <- 1;	just[which(lon>=qlon[4])] <- 0
    
    y[which(lat<qlat[3])]  <- -yNudge;  y[which(lat>=qlat[3])] <- yNudge
    y[which(x!=0)] <- 0			#; print(x); print(y)

    #data.frame(lat=c(lat,lat), lon=c(lon,lon), label=c(rep("aa",cnt),lab), x=c(rep(0,cnt),x), y=c(rep(0,cnt),y), pad=pad)
    data.frame(lat=lat, lon=lon, label=lab, x=c(x), y=c(y), just=just)
}
