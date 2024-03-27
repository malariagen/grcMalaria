###############################################################################
# Read a GRC data file from an Excel sheet
################################################################################
#
# Read sample metadata file
#
grcData.load <- function (grcDataFile, grcDataSheet, species, version) {
    if (!file.exists(grcDataFile)) {
        stop(paste("GRC data file", grcDataFile, "does not exist."))
    }
    if (is.null(version)) {
        stop(paste("You must specify a version number of the GRC data file."))
    }
    verParts <- unlist(strsplit(version,"\\."))
    majorVersion <- verParts[1]
    minorVersion <- verParts[2]
    #
    # The way in which this version checking is written is a bit odd, for historical reasons.
    # In V1, different species had different GRC format versions. Starting with V2, we have a GRC format version, used by all species.
    # Currently we still want to support V1.4 of Pf and 1.0 of Pv. This may change later.
    #
    grc <- NULL
    if (majorVersion >= 2) {
        #
        # Read the GRC file properties
        #
        grcProps <- getGrcProperties (grcDataFile)
        species <- grcProps$species
        version <- grcProps$version
        #
        # Read the GRC file features
        #
        grcFeatures <- readExcelData(grcDataFile, "Features")
        grcFeatures <- grcData.checkColumns (grcFeatures, c("FeatureName","ColumnName","Class","DataType","Reference","Alternative"), "Features", trim=TRUE)
        rownames(grcFeatures) <- grcFeatures$FeatureName
        grcDataSheet <- "Data"
        
    } else {
        if (is.null(species)) {
            stop("Species was not specified for the GRC data file")
        }
        if (species == "Pf") {
            if (version != "1.4") {
                stop(paste("Invalid species/version combination:",species, "/", version))
            }
            # The GRC file features are hardcoded
            grcFeatures <- setup.species.getFeatures.Pf.v1()
            
        } else if (species == "Pv") {
            if (version != "1.0") {
                stop(paste("Invalid species/version combination:",species, "/", version))
            }
            # The GRC file features are hardcoded
            grcFeatures <- setup.species.getFeatures.Pv.v1()
            
        } else {
            stop(paste("Invalid species:",species))
        }
    }
        
    #
    # Read the GRC Sample Data
    # Make sure all the feature columns are present
    #
    grcSampleData <- readExcelData(grcDataFile, grcDataSheet)
    grcSampleData <- grcData.checkColumns (grcSampleData, c("SampleId",grcSampleData$ColumnName), grcDataSheet)
    rownames(grcSampleData) <- grcSampleData$SampleId
    #
    # Process the features to create a configuation that works for this GRC
    #
    speciesConfig <- setup.getSpeciesConfig (species, version, grcFeatures)

    grc <- list(data=grcSampleData, speciesConfig=speciesConfig, version=version)
    grc
}
#
#
#
getGrcProperties <- function (grcDataFile) {
    propsData <- readExcelData(grcDataFile, "Properties")
    propsData <- grcData.checkColumns (propsData, c("Name","Value"), "Properties", trim=TRUE)
    props        <- propsData$Value
    names(props) <- propsData$Name
    #
    fileFormat <- getGrcPropertyValue ("Format", props, "GRC")
    majorVersion <- getGrcPropertyValue ("MajorVersion", props)
    minorVersion <- getGrcPropertyValue ("MinorVersion", props)
    version <- paste(majorVersion, minorVersion, sep=".")
    if (as.integer(majorVersion) > 2) {
        stop(paste0("Unsupported GRC version: ", version))
    }
    if (as.integer(minorVersion) > 0) {
        stop(paste0("Unsupported GRC version: ", version))
    }
    #
    supportedSpecies <- c(
        "Plasmodium falciparum"="Pf", 
        "Plasmodium vivax"="Pv"
    )
    speciesStr <- getGrcPropertyValue ("Species", props, names(supportedSpecies))
    species <- supportedSpecies[speciesStr]
    #
    grcProps <- list(species=species, version=version)
    grcProps
}
#
getGrcPropertyValue <- function (propName, props, expected=NULL) {
    if (!(propName %in% names(props))) {
        stop(paste("Missing GRC file property:", propName))
    }
    value <- props[propName]
    if (!is.null(expected)) {
        if (!(value %in% expected)) {
            stop(paste0("Invalid value spectified GRC file property '",propName,"': ", value))
        }
    }
    value
}
#
grcData.checkColumns <- function (data, columnNames, sheetName, trim=FALSE) {
    missing <- c()
    for (col in columnNames) {
        if (!(col %in% colnames(data))) {
            missing <- c(missing, col)
        } else {
            if (trim) {
                data[,col] <- trimws(data[,col])
            }
        }
    }
    if (length(missing) > 0) {
        missFields <- paste(missing, sep=",")
        stop(paste("Missing column(s) in GRC sheet '",sheetName,"': ",missFields))
    }
    data
}
#
################################################################################
# Merge two GRC datasets
################################################################################
#
# Merge two GRC datasets
#
grcData.merge <- function (inGrc, newGrc, overwrite=FALSE, extendColumns=FALSE) {
    inData  <- inGrc$data
    newData <- newGrc$data

    inSamples  <- rownames(inData);	inFlds  <- colnames(inData)
    newSamples <- rownames(newData);	newFlds <- colnames(newData)
    #
    # Add missing columns in the new data if necessary
    #
    missNewFlds <- inFlds[which(!(inFlds %in% newFlds))]
    if (length(missNewFlds) > 0) {
        missNewData <- data.frame(matrix("<NA>", nrow=length(newSamples), ncol=length(missNewFlds)))
        colnames(missNewData) <- missNewFlds
        rownames(missNewData) <- newSamples
        newData <- cbind(newData, missNewData)
    }
    #
    # Remove extra columns if needed, or add columns introduced by the new dataset
    #
    extraNewFlds <- newFlds[which(!(newFlds %in% inFlds))]
    if (length(extraNewFlds) > 0) {
        if (extendColumns) {
            extraInData <- data.frame(matrix("<NA>", nrow=nrow(inSamples), ncol=length(extraNewFlds)))
            colnames(extraInData) <- extraNewFlds
            rownames(extraInData) <- inSamples
            inData <- cbind(inData, extraInData)
        } else {
            newData <- newData[,-extraNewFlds]
        }
    }
    #
    # Align the columns (they are the same in both datasets)
    #
    inFlds  <- colnames(inData)    
    newData <- newData[,inFlds]
    #
    # Remove overlap in samples before merging
    #
    overlap <- inSamples[which(inSamples %in% newSamples)]
    if (length(overlap) > 0) {
        if (overwrite) {
            inData <- inData[-overlap,]
        } else {
            newData <- newData[-overlap,]
        }
    }
    #
    # Finally, merge the rows
    #
    mergedData <- rbind(inData, newData)
    mergedGrc <- list(data=mergedData, speciesConfig=inGrc$speciesConfig, version=inGrc$version)
    mergedGrc
}

