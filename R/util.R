#
# Utility to make a label corresponding to a minimum identity level (e.g. 0.55 -> "ge0.55")
#
getMinIdentityLabel <- function (minIdentity) {
    minIdentityLabel <- paste("ge", format(minIdentity, digits=2, nsmall=2), sep="")
    minIdentityLabel
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
getOutFolder <- function (ctx, analysisName, subnames=NULL, create=TRUE) {
    cfg <- ctx$config
    sub <- getSubFolder (cfg$folder.out, analysisName, recursive=TRUE, create)
    if (!is.null(subnames)) {
        sub <- getSubFolder (sub, subnames, recursive=TRUE, create)
    }
    sub
}
#
getCacheFolder <- function (ctx, subnames, create=TRUE) {
    sub <- getSubFolder (ctx$config$folder.data, subnames, recursive=TRUE, create)
    sub
}
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
# Read sample data file (first column is an id)
#
readSampleData <- function (dataFile) {
    data <- utils::read.table(dataFile, as.is=TRUE, header=TRUE, quote="", sep="\t", check.names=FALSE)
    sampleNames <- data[,1]
    rownames(data) <- sampleNames
    data <- data[,-1]
    #str(data)
    data
}
#
# Write sample data table file
#
writeSampleData <- function(sampleData, outFilename) {
    writeLabelledData (sampleData, "__Sample", outFilename)
}
#
writeLabelledData <- function(data, idLabel, outFilename) {
    outData <- cbind(rownames(data), data)
    colnames(outData) <- c(idLabel,colnames(data))
    utils::write.table(outData, file=outFilename, sep="\t", quote=FALSE, row.names=FALSE)
}
#
filterSampleData <- function (data, sNames) {
    data <- data[sNames,]
    data
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
writeMatrix <- function (matrixData, matrixDataFile) {
  names <- colnames(matrixData)
  matrixData <- cbind(names, matrixData)
  colnames(matrixData) <- c ("__Sample", names)
  utils::write.table(matrixData, file=matrixDataFile, sep="\t", quote=FALSE, row.names=FALSE)
}
#
filterMatrix <- function (matrixData, sampleNames) {
    matrixData <- matrixData[samplesNames,samplesNames]
    matrixData
}
