graphics.defaultColour <- "gray"
graphics.defaultLcolour <- "black"
graphics.defaultPch <- 21
graphics.defaultSize <- 1
graphics.defaultLwd <- 1
graphics.defaultOrder <- 0

###############################################################################
# Initialization and configuration
################################################################################
#
graphics.loadAttributes <- function (userCtx, attrId, attrField, attrFile, attrSheet) {		#; print(attrId)
    if (is.null(attrId)) {
        stop("A name was not specified for the attribute list, please provide one.")
    }												#; print(colnames(userCtx$unfiltered$meta)); print(attrField)
    if (!(attrField %in% colnames(userCtx$unfiltered$meta))) {
        stop(paste0("The field name specified ('", attrField, "') was not found in the GRC data."))
    }
    if (!file.exists(attrFile)) {
        stop(paste("Attribute list file", attrFile, "does not exist."))
    }
    if (is.null(attrSheet)) {
        stop("A sheet name was not specified for the attribute list file, please provide one.")
    }
    
    attrData <- readExcelData (attrFile, attrSheet)
    checkField <- function (data, field) {
        if (!(field %in% colnames(data))) {
            stop(paste0("Field '", field, "' not found in sheet '", attrSheet, "'."))
        }
    }
    checkField(attrData, "label")
    checkField(attrData, "values")
    checkField(attrData, "order")
    if (!("legendOrder" %in% colnames(attrData))) {
        attrData$legendOrder <- attrData$order
    }
    
    userAttr <- list(id=attrId, field=attrField, data=attrData)

    # Store the attributes list in the context
    userCtx$userGraphicAttributes[[attrId]] <- userAttr
}
#
###############################################################################
# Attributes Resolution
################################################################################
#
graphics.getGraphicalAttributes <- function (ctx, sampleMeta, userAttrNames) {
    #
    # Retrieve the use attributes
    #
    userAttrList <- list()
    for (uIdx in 1 : length(userAttrNames)) {
        userAttrName <- userAttrNames[uIdx]
        userAttr <- ctx$userGraphicAttributes[[userAttrName]]
        if (is.null(userAttr)) {
            stop(paste0("Could not find Graphic Attributes named '", userAttrName, "'"))
        }
        userAttrList[[uIdx]] <- userAttr
    }
    #
    # Make the Graphic Attributes table for the samples
    #
    sampleCount <- nrow(sampleMeta)
    gaData <- data.frame(plot__colour=rep(graphics.defaultColour, sampleCount), 
                         plot__pch=rep(graphics.defaultPch, sampleCount), 
                         plot__size=rep(graphics.defaultSize, sampleCount), 
                         plot__lwd=rep(graphics.defaultLwd, sampleCount),
                         plot__lcolour=rep(graphics.defaultLcolour, sampleCount),
                         plot__order=rep(graphics.defaultOrder, sampleCount))
    rownames(gaData) <- rownames(sampleMeta)
    allLegendData <- data.frame(matrix(nrow=0, ncol=7))
    #
    # Now apply the properties
    #
    userAttrColNames   <- c("colour", "pch", "size", "lwd", "lcolour", "order")
    sampleAttrColNames <- paste("plot", userAttrColNames, sep="__")
    defaultValues <- c(graphics.defaultColour,
                       graphics.defaultPch,
                       graphics.defaultSize,
		       graphics.defaultLwd,
		       graphics.defaultLcolour,
		       graphics.defaultOrder)
    legendDataColNames <- c("label", userAttrColNames)
    
    for (uIdx in 1 : length(userAttrList)) {
         userAttr <- userAttrList[[uIdx]]				#; print(userAttr)
         userAttrId   <- userAttr$id
         fldName      <- userAttr$field
         userAttrData <- userAttr$data					#; print(userAttrData)
         #
         # Put the default properties as the last rule
         #
         userAttrData$values <- as.character(userAttrData$values)	#; print(userAttrData$values)
         if ("<DEFAULT>" %in% userAttrData$values) {
             isDefault <- (userAttrData$values == "<DEFAULT>")		#; print(isDefault); print(userAttrData[!isDefault,]); print(userAttrData[isDefault,])
             userAttrData <- rbind(userAttrData[!isDefault,], userAttrData[isDefault,])
         }
         #
         processed <- rep(FALSE, sampleCount)
         fldValues <- sampleMeta[,fldName]
         legendData <- data.frame(matrix(nrow=0, ncol=7))
         #
         # Apply the rules one by one
         #
         for (rIdx in 1 : nrow(userAttrData)) {
              legendEntry <- c(userAttrData$label[rIdx])
              ruleValuesStr <- userAttrData$values[rIdx]
              matchIdx <- NULL
              if (ruleValuesStr == "<DEFAULT>") {
                  matchIdx <- which(!processed)
              } else if (ruleValuesStr == "<MULTIVALUE>") {
 	          matchIdx <- which(grepl(",", fldValues) & !processed)
              } else {
                  ruleValues <-  trimws(unlist(strsplit(ruleValuesStr, ",")))
                  matchIdx <- which(fldValues %in% ruleValues)
              }
              if (length(matchIdx) == 0) {
                  next
              }
              for (aIdx in 1 : length(userAttrColNames)) {
                  userAttrColName <- userAttrColNames[aIdx]
                  if (!(userAttrColName %in% colnames(userAttrData))) {
                      legendEntry <- c(legendEntry, defaultValues[aIdx])
                      next
                  }
                  newValue <- userAttrData[rIdx,userAttrColName]
                  if (newValue != "-") {
                      sampleAttrColName <- sampleAttrColNames[aIdx]	#; print(sampleAttrColName)
                      gaData[matchIdx,sampleAttrColName] <- newValue
                      legendEntry <- c(legendEntry, newValue)
                  } else {
                      legendEntry <- c(legendEntry, defaultValues[aIdx])
                  }
              }
              processed[matchIdx] <- TRUE
              if (userAttrData$legendOrder[rIdx] > 0) {
                  legendData <- rbind(legendData, legendEntry)
              }
          }
          colnames(legendData) <- legendDataColNames
          legendOrder <- as.numeric(legendData$order)
          legendData <- legendData[order(legendOrder),]
          allLegendData <- rbind(allLegendData, legendData)
    }									#; print(legendData)
    gaData$plot__order   <- as.numeric(gaData$plot__order)
    gaData$plot__colour  <- as.character(gaData$plot__colour)
    gaData$plot__pch     <- as.integer(gaData$plot__pch)
    gaData$plot__size    <- as.numeric(gaData$plot__size)
    gaData$plot__lwd     <- as.numeric(gaData$plot__lwd)
    gaData$plot__lcolour <- as.character(gaData$plot__lcolour)
    #
    colnames(allLegendData) <- legendDataColNames
    #
    list(sampleAttrData=gaData, legendData=allLegendData)
}

###############################################################################
# Plot Configuration
###############################################################################
getGraphicsFilename <- function(graphicFilenameRoot) {
    fileExt <- ".png"    # Default
    return (paste(graphicFilenameRoot, fileExt, sep=""))
}

initializeGraphics <- function(graphicFilename, params) {
    width  <- params$plot.width
    height <- params$plot.height
    units  <- params$plot.units
    dpi   <- params$plot.dpi
    #format <- params$plot.file.format
    defaultFontSize <- 14

    grDevices::png(filename=graphicFilename, width=width, height=height, units=units, res=dpi, 
                   pointsize=defaultFontSize, bg="white")
}

###############################################################################
# Colour Palettes
###############################################################################
graphics.getColourPalette <- function (ctx) {
    userCtx <- ctx$rootCtx
    palette <- userCtx$userColourPalette			#; print(palette)
    if (is.null(palette)) {
        palette <- userCtx$config$defaultPalette
    }								#; print(palette)
    palette
}

graphics.resetColourPalette <- function (ctx) {
    userCtx <- ctx$rootCtx
    if ("userColourPalette" %in% names(userCtx)) {
        rm("userColourPalette", envir=userCtx)
    }
}

graphics.setColourPalette <- function (ctx, palette) {
    userCtx <- ctx$rootCtx
    userCtx$userColourPalette <- palette
}

graphics.makeTextPalette <- function (palette) {
    rgbs <- t(grDevices::col2rgb(palette))
    lum <- rgbs[,1]*0.299 + rgbs[,2]*0.587 + rgbs[,3]*0.114
    textCol <- rep("black",nrow(rgbs))
    textCol[which(lum<130)] <- "white"
    textCol
}
