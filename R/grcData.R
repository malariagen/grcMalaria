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
    grcSampleData <- grcData.readExcel(grcDataFile, grcDataSheet)
    grc <- grcData.standardize(grcSampleData, species, version)
    grc
}

###############################################################################
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

###############################################################################
# Version checking and transformation
################################################################################
#
# Read sample metadata file
#
grcData.standardize <- function (grcSampleData, species, version) {
    speciesConfig <- setup.getSpeciesConfig (species, version)		#; print(species); print(str(speciesConfig))
    
    # TODO - in future, perform structure checks and homogenization, and handle versioning
    grc <- list(data=grcSampleData, speciesConfig=speciesConfig, version=version)
    grc
}

###############################################################################
# Data input routines
################################################################################
#
# Read sample metadata file
#
grcData.readExcel <- function (grcDataFile, grcDataSheet) {
    GRC_SAMPLE_ID_COL  <- "SampleId"				# TODO  This may need to be configured globally
    grcData <- data.frame(readxl::read_excel(grcDataFile, sheet=grcDataSheet, col_names=FALSE, col_types="text", .name_repair="minimal"))
    colnames(grcData) <- grcData[1,]				#; print(str(sampleMeta))
    grcData <- grcData[-1,]
    sampleNames <- grcData[, GRC_SAMPLE_ID_COL]			#; print (sampleNames)
    rownames(grcData) <- sampleNames				#; print(str(sampleMeta))
    grcData
}

