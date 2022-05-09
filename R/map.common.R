
map.colour.border.country <- "black"
map.colour.border.admdiv1 <- "#B4B4B4"
map.colour.land  <- "#F1F1F1"
map.colour.sea   <- "#CDDBDB"
map.colour.river <- "#CDDBDB"

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
    #
    # Read the countries needed in this analysis, so we can get the boundary contours
    #
    unitList <- list()
    cCountry  <- map.getAggregationColumns(0)
    countryValues <- sampleMeta[,cCountry]
    countries <- unique(countryValues)
    for (cIdx in 1:length(countries)) {
        country <- countries[cIdx]					#; print(country)
        countryData <- map.getCountryData (country)			#; print(countryData)
        cMeta <- sampleMeta[which(countryValues == country),]
        #
        # Get a list of provinces (AdmDiv1) in this country (references by GID)
        #
        adm1GIDValues <- as.character(cMeta$AdmDiv1_GID)
        adm1GIDs <- unique(adm1GIDValues)
        unitList[[cIdx]] <- list(name=countryData$CountryName, iso2=countryData$IsoCode2, iso3=countryData$IsoCode3, adm1GIDs=adm1GIDs) 
    }
    #
    # Get the boundaries for all provinces (AdmDiv1) needed for this map, and calculate the bounding box
    #
    library(sp)
    geo <- map.getGeoTables()
    adm1Spdf <- NULL
    xMin <- 1000; xMax <- -1000; yMin <- 1000; yMax <- -1000
    for (cIdx in 1:length(countries)) {
        cl <- unitList[[cIdx]]					#; print(cl);
        #
        # Filter to get the provinces for this country (used for drawing province boundaries)
        #
        cAdm1Lines <- geo$admDiv1.lines[[cl$iso2]]
        #
        # Filter to select the provinces we need in the present analysis (used to define a bounding box)
        #
        anAdm1GIDs <- cl$adm1GIDs						#; print(adm1GIDs)
        anAdm1Lines <- cAdm1Lines[cAdm1Lines@data$GADM_GID_1 %in% anAdm1GIDs,]
        for (idx in 1:nrow(anAdm1Lines)) {
            #
            # Get the bounding box for the province
            # (a 2-column matrix; the first column has the minimum, the second the maximum values; rows represent the spatial dimensions)
            #
            #print(class(anAdm1Lines@polygons))
            #print(class(anAdm1Lines@polygons[[idx]]))
            bbx <- sp::bbox(anAdm1Lines@polygons[[idx]])	#; print(bbx) 
            xMin <- min(xMin,bbx[1,1]); xMax <- max(xMax,bbx[1,2]); yMin <- min(yMin,bbx[2,1]); yMax <- max(yMax,bbx[2,2])	#; print(c(xMin,xMax,yMin,yMax)) 
        }
        #
        # Append the data to those of other countries
        #
        cl$gadmAdm1Data <- anAdm1Lines
        if (is.null(adm1Spdf)) {
            adm1Spdf <- anAdm1Lines
        } else {
            adm1Spdf <- rbind(adm1Spdf, anAdm1Lines)
        }
    }
    adm1_df <- ggplot2::fortify(adm1Spdf)		#; print(colnames(adm1Spdf))
    #
    # Construct the bounding box for this analysis
    # Adjust the bounding box to give some margin
    #
    xMar <- (xMax-xMin)/20;    yMar <- (yMax-yMin)/20
    anBB <- list(xMin=(xMin-xMar), xMax=(xMax+xMar), yMin=(yMin-yMar), yMax=(yMax+yMar))
    anBB$tl <- c(anBB$yMax, anBB$xMin);    anBB$br <- c(anBB$yMin, anBB$xMax)
    anBB$bl <- c(anBB$yMin, anBB$xMin);    anBB$tr <- c(anBB$yMax, anBB$xMax)
    #
    # Get and Crop the country boundaries
    #
    adm0 <- geo$country.lines 
    adm0 <- raster::crop(adm0, raster::extent(anBB$xMin, anBB$xMax, anBB$yMin, anBB$yMax))
    adm0_df <- ggplot2::fortify(adm0)    		#; print(colnames(adm0$spdf))
    #
    rivers <- geo$river.lines 
    rivers <- raster::crop(rivers, raster::extent(anBB$xMin, anBB$xMax, anBB$yMin, anBB$yMax))
    river_df <- NULL
    if (!is.null(rivers)) {
        river_df <- ggplot2::fortify(rivers)		#; print(colnames(rivers))
    }
    #
    lakes <- geo$lake.lines				#; print(lakes)
    lakes <- raster::crop(lakes, raster::extent(anBB$xMin, anBB$xMax, anBB$yMin, anBB$yMax))
    lakes_df <- NULL
    if (!is.null(lakes)) {
        lakes_df <- ggplot2::fortify(lakes)		#; print(colnames(lakes))
    }
    #
    # Construct a base plot for completing subsequent maps
    #
    baseMapPlot <- ggplot2::ggplot(bg=map.colour.sea) +
    	    ggplot2::coord_quickmap(xlim=c(anBB$xMin, anBB$xMax), ylim=c(anBB$yMin, anBB$yMax), expand=FALSE) +
            ggplot2::geom_polygon(data=adm0_df, ggplot2::aes(x=long, y=lat, group=group), 
                                  fill=map.colour.land, col=map.colour.border.country, size=1.5) +
    	    ggplot2::labs(x="Longitude", y="Latitude") +
    	    ggplot2::geom_polygon(data=adm1_df, ggplot2::aes(x=long, y=lat, group=group),
    	                          fill=NA, col=map.colour.border.admdiv1, size=1)
    if (!is.null(river_df)) {	                         
        baseMapPlot <- baseMapPlot +
    	    ggplot2::geom_path(data=river_df, ggplot2::aes(x=long, y=lat, group=group), 
    	                       col=map.colour.river, size=1)
    }
    if (!is.null(lakes_df)) {	                         
        baseMapPlot <- baseMapPlot +
    	    ggplot2::geom_polygon(data=lakes_df, ggplot2::aes(x=long, y=lat, group=group), 
    	                          fill=map.colour.river, col=map.colour.river, size=1)
    }
    baseMapPlot <- baseMapPlot +
    	    ggplot2::theme(panel.background=ggplot2::element_rect(fill=map.colour.sea))
    	                         
    # Return all the elements
    list(baseMap=baseMapPlot, anBB=anBB, adm0_df=adm0_df, adm1_df=adm1_df)
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
map.getAdmDivNames <- function(gids) {		; print(gids)
    geo <- map.getGeoTables()
    admDivs <- geo$admDivs			; print(admDivs)
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
map.computeLabelParams <- function (aggUnitData, baseMapInfo) {
    bbox <- baseMapInfo$anBB
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
