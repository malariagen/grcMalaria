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
getDataFolder <- function (subnames, create=TRUE) {
    sub <- getSubFolder (folder.data, subnames, recursive=TRUE, create)
    sub
}
#
getOutFolder <- function (analysisName, subnames=NULL, create=TRUE) {
    sub <- getSubFolder (folder.results, analysisName, recursive=TRUE, create)
    if (!is.null(subnames)) {
        sub <- getSubFolder (sub, subnames, recursive=TRUE, create)
    }
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
# #####################################################################################
# Palette for pie charts.
# #####################################################################################
#
# This one was submitted by Kevin Wright on https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
# Other similar palettes are shown in the same article
c25Palette <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
#
# #####################################################################################
#  Data I/O Utilities - read/write/filter samples tables data
# #####################################################################################
#
# Read sample data file (first column is an id)
#
readSampleData <- function (dataFile) {
    data <- read.table(dataFile, as.is=TRUE, header=TRUE, quote="", sep="\t", check.names=FALSE)
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
    write.table(outData, file=outFilename, sep="\t", quote=FALSE, row.names=FALSE)
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
    data <- read.table(matrixDataFile, as.is=TRUE, header=TRUE, sep="\t")
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
  write.table(matrixData, file=matrixDataFile, sep="\t", quote=FALSE, row.names=FALSE)
}
#
filterMatrix <- function (matrixData, sampleNames) {
    matrixData <- matrixData[samplesNames,samplesNames]
    matrixData
}
