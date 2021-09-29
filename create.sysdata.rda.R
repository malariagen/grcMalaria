#
# Reconfigure when updating
#
getFilePathInfo <- function (filePath, checkFileExists=FALSE) {
    filePath <- normalizePath(filePath, mustWork=checkFileExists)
    fileDir  <- dirname(filePath)
    fileName <- basename(filePath)
    fileExt  <- tools::file_ext(fileName)
    fileBaseName <- tools::file_path_sans_ext(fileName)
    list(dir=fileDir, name=fileName, base=fileBaseName, ext=fileExt)
}
getDatedFilePath <- function(filePath) {
    fp <- getFilePathInfo(filePath)
    searchPattern <- paste(fp$base, "-\\d{8}.", fp$ext, sep="")
    filesFound <- dir(path=fp$dir, pattern=searchPattern)
    if (length(filesFound) == 0) {
        stop(paste("Dated file ",fp$name," not found in ",fp$dir,sep=""))
    }
    latestFile <- max(filesFound)  # The one with the latest date
    latestFile <- file.path(fp$dir, latestFile)
    latestFile
}
#
# ###########################################################
# Geographical setup data
# ###########################################################
#geoFile <- "data-raw/GeoData.xlsx"
convertUnicodeNames <- function (names) {
    converted <- gsub(">","", gsub("<U\\+", "\\\\u", names))
    result <- stringi::stri_unescape_unicode(converted)
    result
}
#
print("Loading Aggregation Units")
aggUnitsFile <- getDatedFilePath ("./data-raw/annotatedAggregationUnits.tab")
aggUnitData <- read.table(aggUnitsFile, as.is=TRUE, header=TRUE, quote="", sep="\t", encoding="ASCII")
#aggUnitData$GadmName <- convertUnicodeNames(aggUnitData$GadmName)
rownames(aggUnitData) <- aggUnitData$Id
#
print("Loading ISO-3166 Country Data")
countryFile <- "data-raw/ISO-3166.tab"
countryData <- read.table(countryFile, as.is=TRUE, header=TRUE, quote="", na.strings=c(), sep="\t")
rownames(countryData) <- countryData$iso2
#
map.geoTables <- list(gadmUnits=aggUnitData, countries=countryData)
#
# ###########################################################
# Barcode SNP structure
# ###########################################################
# Read and initialize the barcode metadata
# It must be a tab-separated file with one row per SNP, and at least four fields named "Order", "SnpName", "Ref", "Nonref"
readBarcodeMeta <- function(filename) {
    file <- paste("data-raw/barcodeMeta", filename, sep="/")
    meta <- utils::read.table(file, as.is=TRUE, header=TRUE, sep="\t")
    meta <- meta[,c("Order", "SnpName", "Ref", "Nonref")]
    meta <- meta[order(meta$Order),]
    rownames(meta) <- meta$SnpName
    meta
}
#
barcode.metadata.Pf <- list (
    "v1.0" = readBarcodeMeta("BarcodingSnps-Pf-1.0.tab")
)
barcode.metadata.Pv <- list (
    "v1.0" = readBarcodeMeta("BarcodingSnps-Pv-1.0.tab")
)
#
barcode.metadata <- list (
    Pf=barcode.metadata.Pf,
    Pv=barcode.metadata.Pv
)
#
# ###########################################################
# File Storage
# ###########################################################
#
usethis::use_data(map.geoTables, barcode.metadata, internal=TRUE, overwrite=TRUE)
