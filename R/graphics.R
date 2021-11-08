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
graphics.loadAttributes <- function (ctx, attrId, attrField, attrFile, attrSheet) {		#; print(attrId)
    if (is.null(attrId)) {
        stop("A name was not specified for the attribute list, please provide one.")
    }												#; print(colnames(ctx$unfiltered$meta)); print(attrField)
    if (!(attrField %in% colnames(ctx$unfiltered$meta))) {
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
    userAttrList <- ctx$userGraphicAttributes
    if (is.null(userAttrList)) {
        userAttrList <- list()
    }
    userAttrList[[attrId]] <- userAttr
    ctx$userGraphicAttributes <- userAttrList
    ctx
}
#
###############################################################################
# Attributes Resolution
################################################################################
#
graphics.getGraphicalAttributes <- function (sampleMeta, userAttrNames) {
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

initializeGraphics <- function(graphicFilename, widthInch, heightInch, resolution=150, fontSize=14, background="white") {
    grDevices::png(filename=graphicFilename, width=widthInch, height=heightInch, units="in", pointsize=fontSize, bg=background, res=resolution)
}












###############################################################################
###############################################################################
# TO BE FACTORED OUT
###############################################################################
###############################################################################

applyGraphicalAttributes <- function(sampleMetadata, renderDefs) {
    
    numSamples <- nrow(sampleMetadata)
    sampleColour <- rep (analysis.defaultColour,  numSamples)
    samplePch    <- rep (analysis.defaultPch,     numSamples)
    sampleSize   <- rep (analysis.defaultSize,    numSamples)
    sampleLwd    <- rep (analysis.defaultLwd,     numSamples)
    sampleLcolour<- rep (analysis.defaultLcolour, numSamples)
    sampleOrder  <- rep (1000, numSamples)
    
    legendData <- data.frame(matrix(nrow=0, ncol=5))
    
    # Go through the list items in the "render" list one by one
    # This specify a set of graphical attributes to be applied depending on the samples; value in a specific metadata field
    for (renIdx in 1:length(renderDefs)) {
        renDef <- renderDefs[[renIdx]]				#; print(renIdx)
        renField <- renDef$field				#; print(paste("field", renField))
        
        # Get the properties to be applied.
        # This is a list of properties, each property specifying the attributes to be applied when the field value is in the specified set
        renProps <- renDef$properties
        for (propIdx in 1:length(renProps)) {
            renProp <- renProps[[propIdx]]
            
            # "Groups" is the list of values for which the property should be applied
            renGroups <- renProp$groups
            #print(paste("groups", renGroups))
            
            renSampleIdx <- which(sampleMetadata[,renField] %in% renGroups)
            numRendered <- length(renSampleIdx)
            #print(numRendered)
            if (numRendered > 0) {
              
                # We have some samples for which the property should be applied.
                # Apply only the graphical attributes specified.
                # A legend label is mandatory
                legendEntry <- c(renProp$label)
                # Symbol fill colour, if specified
                if (!is.null(renProp$colour)) {
                    sampleColour[renSampleIdx] <- renProp$colour
                    legendEntry <- c(legendEntry, renProp$colour)
                } else {
                    legendEntry <- c(legendEntry, analysis.defaultColour)
                }
                # Symbol shape (pch), if specified
                if (!is.null(renProp$pch)) {
                    samplePch[renSampleIdx] <- renProp$pch
                    legendEntry <- c(legendEntry, renProp$pch)
                } else {
                    legendEntry <- c(legendEntry, analysis.defaultPch)
                }
                # Symbol size factor (pex), if specified
                if (!is.null(renProp$size)) {
                    sampleSize[renSampleIdx] <- renProp$size
                    legendEntry <- c(legendEntry, renProp$size)
                } else {
                    legendEntry <- c(legendEntry, analysis.defaultSize)
                }
                # Symbol line width (lwd), if specified
                if (!is.null(renProp$lwd)) {
                    sampleLwd[renSampleIdx] <- renProp$lwd
                    legendEntry <- c(legendEntry, renProp$lwd)
                } else {
                    legendEntry <- c(legendEntry, analysis.defaultLwd)
                }
                # Symbol line colour, if specified
                if (!is.null(renProp$lcolour)) {
                    sampleLcolour[renSampleIdx] <- renProp$lcolour
                   legendEntry <- c(legendEntry, renProp$lcolour)
                } else {
                    legendEntry <- c(legendEntry, analysis.defaultLcolour)
                }
                # Order of plotting (reverse order actually, lower numbers are at the top of the stack)
                if (!is.null(renProp$order)) {
                    sampleOrder[renSampleIdx] <- renProp$order
                }
                # If the plot has a legend, add this category
                if (!is.null(renDef$legend) && renDef$legend) {
                    legendData <- rbind(legendData, legendEntry)
                }
            }
        }
    }
    headers <- colnames(sampleMetadata)
    sampleMetadata <- cbind(sampleMetadata, sampleColour, samplePch, sampleSize, sampleLwd, sampleLcolour, sampleOrder)
    colnames(sampleMetadata) <- c(headers, "plot__colour", "plot__pch", "plot__size", "plot__lwd", "plot__lcolour", "plot__order")
    colnames(legendData) <- c("label","colour","pch","size","lwd","lcolour")
    
    return (list(meta=sampleMetadata, legend=legendData))
}

###############################################################################
# Automatic attribute setup and rendering
###############################################################################
analysis.autoPalette <- c("red", "orange1", "yellow", "palegreen", 
                       "cyan", "blue", "magenta", "purple", 
                       "pink", "lightseagreen", "lightblue", "plum4", 
                       "darkred", "darkgoldenrod", "green4", "orchid")
analysis.autoPch <- c(rep(21,8),rep(24,8),rep(23,8),rep(22,8))

#
# Resolve all the automatic rendering in the plots
#
resolveAutomaticRenderingInPlots <- function(sampleMetadata, plotDefs) {
    for (plotIdx in 1:length(plotDefs)) {
        plotDef <- plotDefs[[plotIdx]]
        plotDef$render <- resolveAutomaticRendering (sampleMetadata, plotDef$render)
        plotDefs[[plotIdx]] <- plotDef
    }
    plotDefs
}

resolveAutomaticRendering <- function(sampleMetadata, renderDefs) {
    # Go through the list items in the "render" list one by one
    # This specify a set of graphical attributes to be applied depending on the samples; value in a specific metadata field
    for (renIdx in 1:length(renderDefs)) {
        renDef <- renderDefs[[renIdx]]
        # Get the graphical properties, and convert them if they are automatic
        renProps <- renDef$properties
        if ((length(renProps) == 1) && (renProps == "auto")) {
            renDef$properties <- makeAutomaticGraphicalProps (sampleMetadata, renDef)
        }
        renderDefs[[renIdx]] <- renDef
    }
    renderDefs
}

makeAutomaticGraphicalProps <- function(sampleMetadata, autoRenderDef) {
    
    # Parse the automatic rendering definition
    field      <- autoRenderDef$field
    autoParams <- autoRenderDef$autoParams
    numEntries <- autoParams$numEntries
    palette    <- autoParams$palette
    pch        <- autoParams$pch
    
    ignore    <- c()
    if (!is.null(autoParams$ignore)) {
        ignore <- autoParams$ignore
    }
    
    exclude    <- c()
    if (!is.null(autoRenderDef$exclude)) {
        exclude <- autoRenderDef$exclude
    }
    
    # Get the list of unique values in the field, sort them by how common they are
    fieldValues <- sampleMetadata[,field]
    if (length(exclude) > 0) {
        excludeIdx <- which(fieldValues %in% exclude)
        if (length(excludeIdx) > 0) {
            fieldValues <- fieldValues[-excludeIdx]
        }
    }
    if (length(ignore) > 0) {
        ignoreIdx <- which(fieldValues %in% ignore)
        if (length(ignoreIdx) > 0) {
            fieldValues <- fieldValues[-ignoreIdx]
        }
    }
    
    sortedGroups <- names(sort(table(fieldValues),decreasing=TRUE))
    groupCount <- length(sortedGroups)
    topGroupsCount <- min(numEntries, groupCount)
    topGroups <- sortedGroups[1:topGroupsCount]
    otherGroups <- NULL
    if (groupCount > topGroupsCount) {
        otherGroups <- sortedGroups[(topGroupsCount+1):groupCount]
    }
    
    # Create the palette/pch for this plot, if none specified
    if (is.null(palette)) palette <- analysis.autoPalette
    palette <- adjustVectorLength (palette, (topGroupsCount+1))
    
    if (is.null(pch)) pch <- analysis.autoPch
    pch <- adjustVectorLength (pch, (topGroupsCount+1))
    
    # Build up the list of properties now
    propList <- list()
    for (gIdx in 1:topGroupsCount) {
        gName <- topGroups[gIdx]
        prop <- list()
        prop$label <- gName
        prop$groups <- c(gName)
        prop$order <- gIdx

        if (!is.null(autoParams$size)) {
          prop$size <- autoParams$size
        }
        if (!is.null(autoParams$lwd)) {
          prop$lwd <- autoParams$lwd
        }
        if (!is.null(autoParams$lcolour)) {
          prop$lcolour <- autoParams$lcolour
        }
        prop$colour=palette[gIdx]
        prop$pch=pch[gIdx]
        
        propList[[gIdx]] <- prop
    }
    
    # Add an "Others" group if needed
    if (!is.null(otherGroups)) {
      
      prop <- list()
      prop$label  <- "Other"
      prop$groups <- otherGroups
      prop$order  <-  100

      if (!is.null(autoParams$size)) {
        prop$size <- autoParams$size
      }
      if (!is.null(autoParams$lwd)) {
        prop$lwd <- autoParams$lwd
      }
      if (!is.null(autoParams$lcolour)) {
        prop$lcolour <- autoParams$lcolour
      }
      prop$colour <- palette[topGroupsCount+1]
      prop$pch <- pch[topGroupsCount+1]
      
      propList[[topGroupsCount+1]] <- prop
    }
    propList
}
