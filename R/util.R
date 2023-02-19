#
# Utilities to make a label corresponding to a minimum identity level (e.g. 0.55 -> "ge0.55")
# and vice versa
#
MIN_IDENTITY_PREFIX <- "ge"
getMinIdentityLabel <- function (minIdentity, prefix=MIN_IDENTITY_PREFIX) {
    minIdentityLabel <- paste(prefix, format(minIdentity, digits=2, nsmall=2), sep="")
    minIdentityLabel
}
getMinIdentityFromLabel <- function (minIdentityLabel, prefix=MIN_IDENTITY_PREFIX) {
    if (!startsWith(minIdentityLabel, prefix)) {
        return (NA)					#; print("Incorrect prefix")
    }
    minIdentity <- substring(minIdentityLabel,nchar(prefix)+1)		#; print(level)
    as.numeric(minIdentity)
}

# #####################################################################################
# Folder utilities
# #####################################################################################
#
getSubFolder <- function (parent, subnames, recursive=TRUE, create=TRUE) {
    dir <- parent
    for (i in 1 : length(subnames)) {
        name <- subnames[i]
        dir <- paste(dir, name, sep="/")
        if (create) {
            dir.create(dir, recursive=recursive, showWarnings=FALSE)
        }
    }
    dir
}
#
getOutFolder <- function (config, analysisName, subnames=NULL, create=TRUE) {
    sub <- getSubFolder (config$folder.out, analysisName, recursive=TRUE, create)
    if (!is.null(subnames)) {
        sub <- getSubFolder (sub, subnames, recursive=TRUE, create)
    }
    sub
}
#
#getCacheFolder <- function (config, subnames, create=TRUE) {
#    sub <- getSubFolder (config$folder.data, subnames, recursive=TRUE, create)
#    sub
#}
#
# Datafile Naming - prefix with context name so that datafiles are not overwritten
#
getContextCacheFile <- function (ctx, datasetName, subnames, filename) {
    dataFolder <- getContextCacheFolder (ctx, subnames)
    return(paste(dataFolder, "/", datasetName, ".", filename, sep=""))
}
#
getContextCacheFolder <- function (ctx, subnames, create=TRUE) {
    sub <- getSubFolder (ctx$config$folder.data, c(ctx$id, subnames), recursive=TRUE, create)
    sub
}
#
# #####################################################################################
# Vector/string Manipulation utilities
# #####################################################################################
#
# Utility to replicate a vector until it is a sufficent length, and rightsize it
#
adjustVectorLength <- function(vec, len) {
    if (length(vec) < len) {
        numRep <- ceiling(len / length(vec))
        vec <- rep (vec, numRep)
    }
    vec[1:len]
}
#
# A useful function: rev() for strings
#
strReverse <- function(x) {
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}
#
#
#
convertUnicodeNames <- function (names) {
    converted <- gsub(">","", gsub("<U\\+", "\\\\u", names))
    result <- stringi::stri_unescape_unicode(converted)
    result
}
#
# #####################################################################################
#  Data I/O Utilities - read Excel sheet as data frame
# #####################################################################################
#
readExcelData <- function (xFile, xSheet) {
    xData <- data.frame(readxl::read_excel(xFile, sheet=xSheet, col_names=FALSE, col_types="text", .name_repair="minimal"))
    colnames(xData) <- xData[1,]
    xData <- xData[-1,]
    xData
}
#
# #####################################################################################
#  Data I/O Utilities - read/write/filter samples tables data
# #####################################################################################
#
grcMalaria.stored.data <- NULL
#
fileExists <- function (dataFilename, ext) {
    file.exists(paste0(dataFilename,ext))
}
#
rdaFileExists <- function (dataFilename) {
    fileExists (dataFilename, ".rda") 
}
#
# Read sample data file (first column is an id)
#
readSampleData <- function (dataFilename, ext=".tab") {
    data <- utils::read.table(paste0(dataFilename,ext), as.is=TRUE, header=TRUE, quote="", sep="\t", check.names=FALSE)
    sampleNames <- data[,1]
    rownames(data) <- sampleNames
    data <- data[,-1]
    #str(data)
    data
}
#
readRdaSampleData <- function (dataFilename, ext=".rda") {
    load(file=paste0(dataFilename,ext))
    grcMalaria.stored.data
}
#
# Write sample data table file
#
writeSampleData <- function(sampleData, dataFilename, ext=".tab") {
    writeLabelledData (sampleData, "__Sample", paste0(dataFilename,ext))
}
#
writeLabelledData <- function(data, idLabel, dataFilename) {
    outData <- cbind(rownames(data), data)
    colnames(outData) <- c(idLabel,colnames(data))
    utils::write.table(outData, file=dataFilename, sep="\t", quote=FALSE, row.names=FALSE)
}
#
writeRdaSampleData <- function(sampleData, dataFilename, ext=".rda") {
    grcMalaria.stored.data <- sampleData
    save(grcMalaria.stored.data, file=paste0(dataFilename,ext), compress=TRUE, compression_level=9)
}
#
# #####################################################################################
#  Data I/O Utilities - read/write/filter samples matrices
# #####################################################################################
#
# Read distance matrix file
#
readMatrix <- function (matrixDataFile) {
    data <- utils::read.table(matrixDataFile, as.is=TRUE, header=TRUE, sep="\t")
    sampleNames <- data[,1]
    data <- data[,-1]
    rownames(data) <- sampleNames
    colnames(data) <- sampleNames
    #str(data)
    data
}
#
readRdaMatrix <- function (matrixDataFile) {
    load(file=paste0(matrixDataFile,".rda"))
    grcMalaria.stored.data
}
#
writeMatrix <- function (matrixData, matrixDataFile) {
  names <- colnames(matrixData)
  matrixData <- cbind(names, matrixData)
  colnames(matrixData) <- c ("__Sample", names)
  utils::write.table(matrixData, file=matrixDataFile, sep="\t", quote=FALSE, row.names=FALSE)
}
#
writeRdaMatrix <- function(matrixData, outFilename) {
    grcMalaria.stored.data <- matrixData
    save(grcMalaria.stored.data, file=paste0(outFilename,".rda"), compress=TRUE, compression_level=9)
}
#
# #####################################################################################
#  Date and interval calculation
# #####################################################################################
#
# 
#
parseTimeIntervals <- function (timePeriods) {
    intervals <- list()
    if (is.null(timePeriods)) {
         intervals[[1]] <- getDefaultTimeInterval()
    } else {
        for (tpIdx in 1 : length(timePeriods)) {
             tp <- timePeriods[[tpIdx]]
             if (tp$type == "year") {
                 startDate <- as.Date(tp$start, format="%d-%B-%Y")
                 endDate <- startDate + lubridate::years(1) - 1
                 intervals[[tpIdx]] <- list(name=tp$name, start=startDate, end=endDate)
             } else if (tp$type == "period") {
                 startDate <- as.Date(tp$start, format="%d-%B-%Y")
                 endDate <- as.Date(tp$end, format="%d-%B-%Y")
                 intervals[[tpIdx]] <- list(name=tp$name, start=startDate, end=endDate)
             } else {
                 stop(paste("Invalid type of time period:", tp$type))
             }
        }
    }
    intervals
}
#
getDefaultTimeInterval <- function () {
    list(name=NULL, start=NULL, end=NULL)
}
#
# #####################################################################################
#  ggplot2 support functions
# #####################################################################################
#
# This function replaces aes_strng() allowing the use of column names with dashes
#
fn_aesString <- get("aes_string", asNamespace("ggplot2"))
aes_string2 <- function(...){
    args <- lapply(list(...), function(x) sprintf("`%s`", x))
    do.call(fn_aesString, args)
}
