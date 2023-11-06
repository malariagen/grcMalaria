
map.colour.border.country <- "black"
map.colour.border.admdiv1 <- "#B4B4B4"
map.colour.land  <- "#F1F1F1"
map.colour.sea   <- "#CDDBDB"
map.colour.river <- "#CDDBDB"
#
###############################################################################
# Creation of the physical/political background map
################################################################################
#
baseMap.getBaseMap <- function(userCtx, sampleSetName, params) {
    #
    sampleSet <- userCtx$rootCtx$sampleSets[[sampleSetName]]
    ctx <- sampleSet$ctx
    sampleMeta  <- ctx$unfiltered$meta
    #
    # Check if the base map already exists for this aspect ratio; if so, skip creation
    #
    aspectRatioLabel <- paste0("",params$plot.aspectRatio)
    if (aspectRatioLabel %in% names(sampleSet$baseMaps)) {
        return (sampleSet$baseMaps[[aspectRatioLabel]])
    }
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
    # Silly trick to make the package checker happy... :-(
    long <- lat <- group <- NULL
    # The following line loads the sf package, or else we get an error later
    dummy <- sf::st_point(1:2)
    #
    # Get the boundaries for all provinces (AdmDiv1) needed for this map, and calculate the bounding box
    #
    geo <- map.getGeoTables()
    adm1Spdf <- NULL
    xMin <- 1000; xMax <- -1000; yMin <- 1000; yMax <- -1000
    for (cIdx in 1:length(countries)) {
        #
        # Get the list of units (provinces) in this country
        #
        cUnitList <- unitList[[cIdx]]					#; print(cUnitList);
        cIso2 <- cUnitList$iso2
        #
        # Get the provinces for this country (used for drawing province boundaries)
        #
        cAdm1Lines <- geo$admDiv1.lines[[cIso2]]			#; print(class(cAdm1Lines))
        #
        # Filter to select the provinces we need in the present analysis (used to define a bounding box)
        #
        anAdm1GIDs <- cUnitList$adm1GIDs				#; print(adm1GIDs)
        anAdm1Lines <- cAdm1Lines[which(cAdm1Lines$GADM_GID_1 %in% anAdm1GIDs),]	#; print(class(cAdm1Lines))
        #
        # TEMPORARY FIX (hopefully)
        # Not all GADM_GID_1 identifiers have an entry in the naturalearth dataset. For this reason, anAdm1Lines may
        # contain fewer rows, or none.This is not something we can fix in a hurry, if at all.
        # It affects two thing: the calculation of the bounding box, and the drawing of the provice boundaries.
        # Here, we fix the bounding box which otherwise causes a fatal error; the drawing of the province boundaries
        # cannot be solved as easily, and for now the province boundaries will not be drawn for the provinces affected.
        #
        if (!is.null(anAdm1Lines) && (nrow(anAdm1Lines) == length(anAdm1GIDs))) {
            # All the rows are there, we can calculate the bounding box normally
            bbx <- sf::st_bbox(anAdm1Lines)
            xMin <- min(xMin,bbx$xmin); yMin <- min(yMin,bbx$ymin) 
            xMax <- max(xMax,bbx$xmax); yMax <- max(yMax,bbx$ymax)	#; print(c(xMin,xMax,yMin,yMax)) 
        } else {
            #
            # Remedial code: we roughly estimate the bounding boxes by adding a certain amount to the "notional"
            # latitude and longitude of the admdiv (estimated at 0.7 degrees for now)
            #
            adj <- 0.7
            divs <- geo$admDivs
            divs <- divs[which(divs$GID %in% anAdm1GIDs),]
            lat <- divs$Latitude; lon <- divs$Longitude
            xMin <- min(xMin,(min(lon)-adj)); yMin <- min(yMin,(min(lat)-adj)) 
            xMax <- max(xMax,(max(lon)+adj)); yMax <- max(yMax,(max(lat)+adj))	#; print(c(xMin,xMax,yMin,yMax)) 
        }
        #
        # Append the data to those of other countries
        #
        cUnitList$gadmAdm1Data <- anAdm1Lines
        if (is.null(adm1Spdf)) {
            adm1Spdf <- cAdm1Lines
        } else {
            adm1Spdf <- rbind(adm1Spdf, cAdm1Lines)
        }
    }
    adm1_df <- suppressMessages(ggplot2::fortify(adm1Spdf))	#; print(colnames(adm1Spdf@data)); print(colnames(adm1_df))
    #
    # Create a bounding box, specifying WGS84 (EPSG:4326) to be the coordinates system 
    #
    CRS.WGS84 <- geo$crs
    anBB <- baseMap.adjustMapBoundingBox(xMin, xMax, yMin, yMax, params$plot.aspectRatio)	#; print(anBB)
    anBBCoords <- list(matrix(c(
                      anBB$xMin, anBB$yMin,
                      anBB$xMin, anBB$yMax,
                      anBB$xMax, anBB$yMax,
                      anBB$xMax, anBB$yMin,
                      anBB$xMin, anBB$yMin), ncol = 2, byrow = TRUE))
    anBBExt <- sf::st_polygon(anBBCoords)
    anBBExt <- sf::st_geometry(anBBExt)
    anBBExt <- sf::st_set_crs(anBBExt, CRS.WGS84)
    #
    # Get and Crop the country boundaries
    #
    adm0 <- sf::st_make_valid(geo$country.lines)						#; print("adm0"); print(colnames(adm0))
    adm0 <- sf::st_transform(adm0, CRS.WGS84)
    adm0 <- suppressWarnings(sf::st_intersection(adm0, anBBExt))
    adm0_df <- ggplot2::fortify(adm0)    							#; print(colnames(adm0_df))
    #
    rivers <- sf::st_make_valid(geo$river.lines)						#; print("rivers"); print(colnames(rivers))
    rivers <- sf::st_transform(rivers, CRS.WGS84)
    rivers <- suppressWarnings(sf::st_intersection(rivers, anBBExt))
    river_df <- NULL
    if (!is.null(rivers) & (nrow(rivers)>0)) {
        river_df <- ggplot2::fortify(rivers)	    						#; print(colnames(river_df))				
    }
    #
    lakes <- sf::st_make_valid(geo$lake.lines)							#; print("lakes"); print(colnames(lakes))
    lakes <- sf::st_transform(lakes, CRS.WGS84)
    lakes <- suppressWarnings(sf::st_intersection(lakes, anBBExt))
    lakes_df <- NULL
    if (!is.null(lakes) & (nrow(lakes)>0)) {
        lakes_df <- ggplot2::fortify(lakes)							#; print(colnames(lakes_df))
    }
    #
    mapLineWidth    <- baseMap.getBordersThickness (anBB, params)	#; print(mapLineWidth)
    adm0BorderWidth <- mapLineWidth * params$adm0BorderWidth	#; print(adm0BorderWidth)
    adm1BorderWidth <- mapLineWidth * params$adm1BorderWidth	#; print(adm1BorderWidth)
    riverWidth      <- mapLineWidth * params$riverWidth		#; print(riverWidth)
    lakeshoreWidth  <- mapLineWidth * params$lakeshoreWidth	#; print(lakeshoreWidth)
    #
    # Construct a base plot for completing subsequent maps
    #
    baseMapPlot <- ggplot2::ggplot(bg=map.colour.sea) +
            ggplot2::geom_sf(data=adm0_df, fill=map.colour.land, col=NA) +
    	    ggplot2::labs(x="Longitude", y="Latitude") +
            ggplot2::geom_sf(data=adm1_df, fill=NA, col=map.colour.border.admdiv1, linewidth=adm1BorderWidth) +
            ggplot2::geom_sf(data=adm0_df, fill=NA, col=map.colour.border.country, linewidth=adm0BorderWidth)
    if (!is.null(river_df)) {	                         
        baseMapPlot <- baseMapPlot +
            ggplot2::geom_sf(data=river_df, col=map.colour.river, linewidth=riverWidth)
    }
    if (!is.null(lakes_df)) {	                         
        baseMapPlot <- baseMapPlot +
            ggplot2::geom_sf(data=lakes_df, fill=map.colour.river, col=map.colour.river, linewidth=lakeshoreWidth)
    }
    baseMapPlot <- baseMapPlot +
            ggplot2::coord_sf(xlim=c(anBB$xMin,anBB$xMax), ylim=c(anBB$yMin,anBB$yMax), expand=FALSE) +
    	    ggplot2::theme(panel.background=ggplot2::element_rect(fill=map.colour.sea))
    # Return all the elements
    baseMapInfo <- list(baseMap=baseMapPlot, anBB=anBB, adm0_df=adm0_df, adm1_df=adm1_df)
    sampleSet$baseMaps[[aspectRatioLabel]] <- baseMapInfo
    baseMapInfo
}
#
# Construct the bounding box
#
baseMap.adjustMapBoundingBox <- function (xMin, xMax, yMin, yMax, aspectRatio) {
    #
    # Adjust the bounding box to give some margin
    #
    xMar <- (xMax-xMin)/20;    yMar <- (yMax-yMin)/20
    xMin <- xMin-xMar; xMax <- xMax+xMar; yMin <- yMin-yMar; yMax <- yMax+yMar
    #
    # Now extend the bounding box so that it has the correct aspect ratio
    #
    mapW <- xMax - xMin; mapH <- yMax - yMin
    mapAR <- mapW / mapH
    if (mapAR < aspectRatio) {
        newW <- mapH * aspectRatio
        extraW <- newW - mapW
        xMin <- xMin - (extraW / 2); xMax <- xMax + (extraW / 2)
    } else if (mapAR > aspectRatio) {
        newH <- mapW / aspectRatio
        extraH <- newH - mapH
        yMin <- yMin - (extraH / 2); yMax <- yMax + (extraH / 2)
    }
    #
    # Create the bounding box
    #
    bb <- list(xMin=xMin, xMax=xMax, yMin=yMin, yMax=yMax, 
               tl=c(yMax, xMin), br=c(yMin, xMax), bl=c(yMin, xMin), tr=c(yMax, xMax))
    bb
}
#
#
#
baseMap.getBordersThickness <- function(bbox, params) {
    plotWidth <- bbox$xMax - bbox$xMin						#; print(plotWidth)
    borderThickness <- params$baseBorderWidth / sqrt(plotWidth / 3.0)		#; print(borderThickness)
    borderThickness
}
