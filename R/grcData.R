
grcData.load <- function (grcDataFile, grcDataSheet, species, version) {
    if (!file.exists(grcDataFile)) {
        stop(paste("GRC data file", grcDataFile, "does not exist."))
    }
    grcSampleData <- grcData.readExcel(grcDataFile, grcDataSheet)
    grc <- grcData.standardize(grcSampleData, species, version)
    grc
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

    grcData <- data.frame(readxl::read_excel(grcDataFile, sheet=grcDataSheet, col_names=FALSE, col_types="text"))
    colnames(grcData) <- grcData[1,]				#; print(str(sampleMeta))
    grcData <- grcData[-1,]
    sampleNames <- grcData[, GRC_SAMPLE_ID_COL]			#; print (sampleNames)
    rownames(grcData) <- sampleNames				#; print(str(sampleMeta))
    grcData
}

