
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
    # Silly trick to make the package checker happy... :-(
    long <- lat <- group <- NULL
    # The following line loads the sf package, or else we get an error later
    dummy <- sf::st_point(1:2)
    #
    # Get all the possible marker positions for this map and work out the boundaries
    #
    geo <- map.getGeoTables()
    admDivInfo <- geo$admDivs
    mkGids <- unique(c(sampleMeta$AdmDiv1_GID, sampleMeta$AdmDiv2_GID))
    mkLat <- admDivInfo[mkGids,"Latitude"]
    mkLon <- admDivInfo[mkGids,"Longitude"]
    #
    # Work out the extent of the markers, and add a margin to determine a bounding box
    #
    xMin <- min(mkLon); xMax <- max(mkLon);
    yMin <- min(mkLat); yMax <- max(mkLat); 			#; print(paste(xMin, xMax, yMin, yMax))
    marginProp <- 0.15
    bb <- baseMap.adjustMapBoundingBox(xMin, xMax, yMin, yMax, params$plot.aspectRatio, marginProp)	#; print(bb)
    bbCoords <- list(matrix(c(bb$xMin, bb$yMin,
                              bb$xMin, bb$yMax,
                              bb$xMax, bb$yMax,
                              bb$xMax, bb$yMin,
                              bb$xMin, bb$yMin), ncol = 2, byrow = TRUE))
    CRS.WGS84 <- geo$crs
    bbPoly <- sf::st_polygon(bbCoords)
    bbExt <- sf::st_geometry(bbPoly)
    bbExt <- sf::st_set_crs(bbExt, CRS.WGS84)
    #
    # Get and Crop the layers
    #
    constructAnnotatedLayer <- function(sfObj, bb, bbExt) {				#;print(nrow(sfObj))
        xmin <- sfObj$xMin; xmax <- sfObj$xMax; ymin <- sfObj$yMin; ymax <- sfObj$yMax	#;print(head(xmin));print(head(xmax));print(head(ymin));print(head(ymax))
        include <- !((xmin>bb$xMax)|(xmax<bb$xMin)|(ymin>bb$yMax)|(ymax<bb$yMin))	#;print(head(include))
        sfObj <- sfObj[include,]							#;print(nrow(sfObj))
        layer <- sf::st_transform(sfObj, CRS.WGS84)
        df <- NULL
        if (!is.null(layer) & (nrow(layer)>0)) {
            df <- ggplot2::fortify(layer)				
        }
        df
    }
    adm0_df  <- constructAnnotatedLayer (geo$ann.country.lines, bb, bbExt)
    adm1_df  <- constructAnnotatedLayer (geo$ann.admDiv1.lines, bb, bbExt)
    river_df <- constructAnnotatedLayer (geo$ann.river.lines,   bb, bbExt)
    lakes_df <- constructAnnotatedLayer (geo$ann.lake.lines,    bb, bbExt)
    #
    mapLineWidth    <- baseMap.getBordersThickness (bb, params)	#; print(mapLineWidth)
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
            ggplot2::coord_sf(xlim=c(bb$xMin,bb$xMax), ylim=c(bb$yMin,bb$yMax), expand=FALSE) +
    	    ggplot2::theme(panel.background=ggplot2::element_rect(fill=map.colour.sea))
    # Return all the elements
    baseMapInfo <- list(baseMap=baseMapPlot, anBB=bb)
    sampleSet$baseMaps[[aspectRatioLabel]] <- baseMapInfo
    baseMapInfo
}
#
# Construct the bounding box
#
baseMap.adjustMapBoundingBox <- function (xMin, xMax, yMin, yMax, aspectRatio, marginProp) {
    #
    # Adjust the bounding box to give some margin
    #
    xMar <- (xMax-xMin)*marginProp;    yMar <- (yMax-yMin)*marginProp
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
    bb <- list(xMin=xMin, xMax=xMax, yMin=yMin, yMax=yMax)
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
