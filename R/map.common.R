default.map.adminLevelColumns <- c("Country","AdmDiv1","AdmDiv2")
default.map.locationColumns <- c("Location","Latitude","Longitude")
default.map.aggregateCountMin <- 5
default.map.size <- c(12,12)

###############################################################################
# Common Routines for Map Generation.
################################################################################
#
map.execute <- function(analysisContext, analysisName, mapType, aggregation, measures, params) {
    parts <- unlist(strsplit(mapType,"\\."))	#; print(mapType)
    mapType <- parts[1]				#; print(mapType)
    visType <- parts[2]				#; print(visType)

    if (mapType %in% c("diversity", "drug", "mutation")) {
        markerMap.execute (analysisContext, analysisName, mapType, visType, aggregation, measures, params)
        
    } else if (mapType == "connect") {
        connectMap.execute (analysisContext, analysisName, mapType, visType, aggregation, measures, params)
        
    } else if (mapType %in% c("haploFreq","haploShare")) {
        haploMap.execute (analysisContext, analysisName, mapType, visType, aggregation, measures, params)

    } else {
        stop(paste("Invalid map type:", mapType))
    }
}

#
###############################################################################
# Creation of the physical/political background map
################################################################################
#
map.buildBaseMap <- function(ctx, analysisName, sampleMeta, dataFolder, params) {

    # Get relevant column names
    adminLevelCols  <- analysis.getParam ("map.adminLevelColumns", params, default.map.adminLevelColumns)
    locationCols    <- analysis.getParam ("map.locationColumns",   params, default.map.locationColumns)
    cCountry  <- adminLevelCols[1]; cProvince <- adminLevelCols[2]; cDistrict <- adminLevelCols[3]
    cLocation <- locationCols[1];   cLat      <- locationCols[2];   cLon      <- locationCols[3]

    # Now read the countries and provinces so we can get the contours from GADM
    unitList <- list()
    countryValues <- sampleMeta[,cCountry]
    countries <- unique(countryValues)
    for (cIdx in 1:length(countries)) {
        country <- countries[cIdx]					#; print(country)
        countryData <- getCountryData (country)				#; print(countryData)
        cMeta <- sampleMeta[which(countryValues == country),]
        provValues <- as.character(cMeta[,cProvince])
        prov <- unique(provValues)					#; print(prov)
        gadmProv <- getGADMNames (country, prov)			#; print(gadmProv)
        unitList[[country]] <- list(name=countryData$name, iso2=countryData$iso2, iso3=countryData$iso3, adm1Names=prov, gadmAdm1Names=gadmProv)
    }
    
    # Read the country borders for the countries involved
    gadmFolder <- getDataFolder(ctx, c("map", "gadm"))
    gadmFolder <- paste(gadmFolder, "/", sep="")
    cIso3 <- iso2ToIso3 (countries)
    gadm0 <- GADMTools::gadm_sp_loadCountries(cIso3, level=0, basefile=gadmFolder)
    
    # Read the province borders for the countries involved
    gadm1Spdf <- NULL
    gadmBB <- list(xMin=1000, xMax=-1000, yMin=1000, yMax=-1000)

    for (cIdx in 1:length(countries)) {
        country <- countries[cIdx]
        cl <- unitList[[country]]							#; print(country)
        cGadm1 <- GADMTools::gadm_sp_loadCountries(cl$iso3, level=1, basefile=gadmFolder)	#; print(1)
        # Select the provinces we need
        cGadm1 <- GADMTools::gadm_subset(cGadm1, level=1, regions=cl$gadmAdm1Names)		#; print(cGadm1)
        # Append the data to that of other countries
	    cl$gadmAdm1Data <- cGadm1							#; print(3)
	    if (is.null(gadm1Spdf)) {
	    gadm1Spdf <- cGadm1$spdf
	    } else {
	    gadm1Spdf <- rbind(gadm1Spdf, cGadm1$spdf)
	    }										#; print(4)
        # Get the bounding box and merge with the other countries
        bb <- GADMTools::gadm_getBbox (cGadm1)							#; print(5)
        gadmBB <- list(xMin=min(gadmBB$xMin,bb[1]), xMax=max(gadmBB$xMax,bb[3]), 
                       yMin=min(gadmBB$yMin,bb[2]), yMax=max(gadmBB$yMax,bb[4]))	#; print(gadmBB)
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
    gadm0_df <- ggplot2::fortify(gadm0$spdf)
    gadm1_df <- ggplot2::fortify(gadm1Spdf)
        
    # Get the background map, and adjust coordinates
    bgMap <- OpenStreetMap::openmap(gadmBB$tl, gadmBB$br, zoom=NULL, type=c("nps"), mergeTiles=TRUE)
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

###############################################################################
# Aggregation of site data
################################################################################
#
map.getAggregationUnitData <- function(aggLevel, analysisContext, analysisName, mapType, params, dataFolder) {

    sampleMeta   <- analysisContext$meta
    barcodeData  <- analysisContext$barcodes
    distData     <- analysisContext$distance

    adminLevelCols  <- analysis.getParam ("map.adminLevelColumns", params, default.map.adminLevelColumns)
    aggLevelIdx <- aggLevel + 1

    # Create aggregation unit id; this is unique for each aggregation unit 
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta, params)
    
    # Get all aggregation units, in order
    aggregateCountMin <- analysis.getParam ("map.aggregateCountMin", params, default.map.aggregateCountMin)
    aggUnitCounts <- table(aggIndex)
    aggUnits <- names(aggUnitCounts[aggUnitCounts >= aggregateCountMin])	#; print(aggUnits)
    aggUnits <- aggUnits[order(aggUnits)]					#; print(aggUnits)
    
    # Get the data for all aggregation units
    aggUnitData <- NULL
    
    for (aIdx in 1:length(aggUnits)) {

        # Get the sample data to be aggregated for this unit
        aggUnit <- aggUnits[aIdx]
        aggSamplesMeta <- sampleMeta[which(aggIndex == aggUnit),]		#; print(nrow(aggSamplesMeta))
        aggSamples <- rownames(aggSamplesMeta)
        aggBarcodes <- barcodeData[aggSamples,]
        aggDist <- distData[aggSamples,aggSamples]
        
        # Get the admin division values from the first sample of this unit (assuming the values are the same for all)
        cValues <- c(aggUnit)
        
        admNames <- adminLevelCols[1:aggLevelIdx]
        cValues <- c(cValues, as.character(aggSamplesMeta[1,admNames]))
        
        cValues <- c(cValues, mean(as.numeric(aggSamplesMeta$Latitude)))
        cValues <- c(cValues, mean(as.numeric(aggSamplesMeta$Longitude)))	#; print(cValues)
        cValues <- c(cValues, nrow(aggSamplesMeta))				#; print(cValues)

        aggUnitData <- rbind(aggUnitData, cValues)
    }
    aggUnitData <- data.frame(aggUnitData)
    numericColIdx <- c((aggLevelIdx+2):ncol(aggUnitData))
    aggUnitData[, numericColIdx] <- sapply(aggUnitData[, numericColIdx], as.numeric)
    rownames(aggUnitData) <- aggUnits
    colnames(aggUnitData) <- c("UnitId", c("Country","AdmDiv1","AdmDiv2")[1:aggLevelIdx],
                               "Latitude","Longitude","SampleCount")		#; print(aggUnitCnames)

    # Write out the aggregation unit data to file
    aggDataFilename  <- paste(dataFolder, "/AggregationUnits-", analysisName, "-", aggLevel, ".tab", sep="")
    write.table(aggUnitData, file=aggDataFilename, sep="\t", quote=FALSE, row.names=FALSE)

    aggUnitData
}

map.getAggregationUnitIds <- function(aggLevel, sampleMeta, params) {
    # Get relevant column names
    adminLevelCols  <- analysis.getParam ("map.adminLevelColumns", params, default.map.adminLevelColumns)
    locationCols    <- analysis.getParam ("map.locationColumns",   params, default.map.locationColumns)
    cCountry  <- adminLevelCols[1]; cProvince <- adminLevelCols[2]; cDistrict <- adminLevelCols[3]
    cLocation <- locationCols[1];   cLat <- locationCols[2];        cLon <- locationCols[3]
    
    # Create aggregation unitI IDs; this is unique for each aggregation unit 
    aggUnitId <- sampleMeta[,cCountry]
    if (aggLevel > 0) {
        aggUnitId <- paste(aggUnitId, sampleMeta[,cProvince], sep="__")
    }
    if (aggLevel > 1) {
        aggUnitId <- paste(aggUnitId, sampleMeta[,cDistrict], sep="__")
    }										#; print(head(aggIndex))
    aggUnitId
}

map.geoTables <- NULL

map.getGeoTables <- function () {
    if (is.null(map.geoTables)) {
        geoFile <- paste(folder.code, "GeoData.xlsx", sep="/")
        
        print("Loading GADM province names")
        gadmProvTable <- data.frame(read_excel(geoFile, sheet="GADM-provinces", col_types="text"))
        
        print("Loading Country Data")
        countryDataTable <- data.frame(read_excel(geoFile, sheet="ISO-3166", col_types="text"))
        rownames(countryDataTable) <- countryDataTable$iso2

        map.geoTables <<- list(gadmProv=gadmProvTable, countries=countryDataTable)
    }
    map.geoTables
}
getGADMNames <- function (country, provinces) {
    geo <- map.getGeoTables()
    gadmProvTable <- geo$gadmProv
    countryTable <- gadmProvTable[which(gadmProvTable$Iso2==country),]
    rownames(countryTable) <- as.character(countryTable$AdmDiv1)
    gadmNames <- countryTable[provinces,"GADMname"]			; print(provinces); print(gadmNames)
    gadmNames
}
iso2ToIso3 <- function (iso2Countries) {
    geo <- map.getGeoTables()
    countryTable <- geo$countries
    iso3 <- countryTable[iso2Countries,"iso3"]
    iso3
}
getCountryData <- function (iso2Country) {
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
map.computeLabelParams <- function (aggUnitData, aggColName, baseMapInfo) {

    bbox <- baseMapInfo$gadmBB
    xNudge <- (bbox$xMax - bbox$xMin)/15
    yNudge <- (bbox$yMax - bbox$yMin)/15
    
    cnt <- nrow(aggUnitData)
    lon <- aggUnitData$Longitude;
    lat <- aggUnitData$Latitude
    lab <- aggUnitData[,aggColName]
    
    x <- y <- vector(mode="integer", length=cnt)
    just <- rep(0.5,cnt)
    qlon <- quantile(lon)		#; print (qlat)
    qlat <- quantile(lat)		#; print (qlon)
    x[which(lon<=qlon[2])] <- -xNudge;	x[which(lon>=qlon[4])] <- xNudge
    just[which(lon<=qlon[2])] <- 1;	just[which(lon>=qlon[4])] <- 0
    
    y[which(lat<qlat[3])]  <- -yNudge;  y[which(lat>=qlat[3])] <- yNudge
    y[which(x!=0)] <- 0			#; print(x); print(y)

    #data.frame(lat=c(lat,lat), lon=c(lon,lon), label=c(rep("aa",cnt),lab), x=c(rep(0,cnt),x), y=c(rep(0,cnt),y), pad=pad)
    data.frame(lat=lat, lon=lon, label=lab, x=c(x), y=c(y), just=just)
}
