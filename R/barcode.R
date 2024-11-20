###############################################################################
# Caching data files
################################################################################
barcode.getBarcodeDataFile <- function (ctx, datasetName) {
    dataFile <- getContextCacheFile(ctx, datasetName, "barcode", "barcodeAlleles")
    dataFile
}
barcode.getBarcodeSeqFile <- function (ctx, datasetName) {
    seqFile  <- getContextCacheFile(ctx, datasetName, "barcode", "barcodeSeqs")
    seqFile
}

barcode.getMetadata <- function (ctx, barcodeData) {
    barcodeSnps <- colnames(barcodeData)
    barcodeMeta <- ctx$config$barcodeMeta[barcodeSnps,]
    barcodeMeta
}

###############################################################################
# Barcode setting and caching
################################################################################
barcode.setDatasetBarcodes <- function (ctx, datasetName, barcodes, store=TRUE) {
    dataset <- ctx[[datasetName]]
    dataset$barcodes <- barcodes
    if (store) {
        barcodeDataFile <- barcode.getBarcodeDataFile (ctx, datasetName)
        writeRdaSampleData(barcodes, barcodeDataFile)
        barcodeSeqFile <- barcode.getBarcodeSeqFile (ctx, datasetName)
        barcode.writeFasta (barcodes, barcodeSeqFile)
    }
    #ctx[[datasetName]] <- dataset
    ctx
}

###############################################################################
# Barcode retrieval, validating and filtering
################################################################################
barcode.initializeBarcodes <- function (ctx, datasetName) {
    config <- ctx$config				#; print(config)
    dataset <- ctx[[datasetName]]
    barcodeDataFile <- barcode.getBarcodeDataFile (ctx, datasetName)

    if (rdaFileExists(barcodeDataFile)) {
        barcodeData <- readRdaSampleData (barcodeDataFile)
        barcode.setDatasetBarcodes (ctx, datasetName, barcodeData, store=FALSE)
    } else {
        #
        # Get barcode alleles, and discard samples that have too much missingness
        #
        barcodeMeta <- config$barcodeMeta		#; print (barcodeMeta)
        barcodeData <- barcode.getBarcodeAlleles (dataset$meta, barcodeMeta)
        print(paste("Barcode alleles - Samples:", nrow(barcodeData), "x SNPs:", ncol(barcodeData)))
        print("Validating barcodes")
        barcode.validateBarcodeAlleles (barcodeData, barcodeMeta)
        #
        # Filter the barcodes by typability, trying to throw away as little as possible
        #
        barcode.writeBarcodeStats (ctx, datasetName, barcodeData, "noFiltering")
        minSampleTypability <- config$minSampleTypability
        minSnpTypability    <- config$minSnpTypability
        print(paste("Filtering barcodes by typability (samples:", minSampleTypability, ", SNPs:", minSnpTypability, ")",sep=""))
        #
        # 1) remove all samples with <0.5 typability, so they affect less the removal of SNPs
        #
        filteredData <- barcode.filterByTypability (barcodeData, bySnp=FALSE, minTypability=0.5)
        #
        # 2) Refine further, using the thresholds specified
        #
        filteredData <- barcode.filterByTypability (filteredData, bySnp=TRUE,  minTypability=minSnpTypability)
        filteredData <- barcode.filterByTypability (filteredData, bySnp=FALSE, minTypability=minSampleTypability)
        barcodeData <- filteredData
        #
        barcode.writeBarcodeStats (ctx, datasetName, barcodeData,
                           paste("filtered-snps_", minSnpTypability, "-samples_", minSampleTypability, sep=""))
        print(paste("Barcode alleles after filtering - Samples:", nrow(barcodeData), "x SNPs:", ncol(barcodeData)))
        #
        barcode.setDatasetBarcodes (ctx, datasetName, barcodeData, store=TRUE)
    }

    # Report missingness
    totalCalls <- nrow(barcodeData) * ncol(barcodeData)
    callCounts <- table(unlist(barcodeData))
    missingCalls <- callCounts["X"]
    missing <- missingCalls/totalCalls
    het <- callCounts["N"]/(totalCalls-missingCalls)
    print(paste("Missing:", missing, "- Het:", het))

    ctx
}
#
# Convert barcodes into a dataframe of alleles
#
barcode.getBarcodeAlleles <- function(sampleMetadata, barcodeFeatures) {	#;print (head(barcodeFeatures))

    colNames <- barcodeFeatures$ColumnName
    alleleData <- data.frame(SampleId=rownames(sampleMetadata))
    for (i in 1:length(colNames)) {
        colName <- colNames[i]
        genos <- toupper(as.character(sampleMetadata[,colName]))
        genos[which(genos=="*")] <- "N"
        genos[which(genos=="-")] <- "X"
        genos[which(genos=="<NA>")] <- "X"
        alleleData <- cbind(alleleData, genos)
    }
    alleleData <- alleleData[,2:ncol(alleleData)]
    rownames(alleleData) <- rownames(sampleMetadata);
    colnames(alleleData) <- colNames;
    
    missingCount <- rowSums(alleleData=="X")
    validIdx <- which(missingCount < ncol(alleleData))
    alleleData <- alleleData[validIdx,]

    alleleData
}
#
# Verify all alleles extracted are valid
#
barcode.validateBarcodeAlleles <- function (barcodeData, barcodeMeta) {
    snpCount <- nrow(barcodeMeta)				#; print(snpCount)
    if (ncol(barcodeData) != snpCount) {
        stop (paste("Number of barcode SNPs in the metadata (",snpCount,") does not match the length of the barcodes (",ncol(barcodeData),")",sep=""))
    }
    validAlleles <- c("A","C","G","T","X","N")
    errorCount <- 0
    for (sIdx in 1:snpCount) {
        calls <- barcodeData[,sIdx]				#; print(calls)
        snpMeta <- barcodeMeta[sIdx,]				#; print(snpMeta)
        snpAlleles <- c(snpMeta$Reference, snpMeta$Alternative,
        		"X","N")				#; print(snpAlleles)
        badIdx <- which(!(calls %in% snpAlleles))		#; print(badIdx)
        if (length(badIdx) > 0) {
            for (bIdx in badIdx) {
                badAllele <- calls[bIdx]
                badSample <- rownames(barcodeData)[bIdx]
                isNt <- (badAllele %in% validAlleles)
                if (isNt) {
                    problem <- "Info: Unexpected allele"
                } else {
                    errorCount <- errorCount+1
                    problem <- "Error: Invalid symbol"
                }
                cat (paste0(problem," found at SNP #",sIdx," in sample ",badSample,": found ",badAllele), fill=TRUE)
            }
        }
    }
    if (errorCount > 0) {
        stop ("Validation errors were detected in the barcodes, the process will be stopped.")
    }
}
#
###############################################################################
# Barcode Sample/SNP filtering
################################################################################
barcode.filterByTypability <- function(barcodeData, bySnp=FALSE, minTypability=0.75) {
    stats <- barcode.computeBarcodeStats (barcodeData, bySnp)
    result <- barcode.filterByTypabilityStats (barcodeData, stats, bySnp, minTypability)
    result
}

barcode.filterByHeterozygosity <- function(barcodeData, bySnp=FALSE, maxHeterozygosity=0.0) {
    stats <- barcode.computeBarcodeStats (barcodeData, bySnp)
    result <- barcode.filterByHeterozygosityStats (barcodeData, stats, bySnp, maxHeterozygosity)
    result
}

barcode.filterByStats <- function(barcodeData, bySnp=FALSE, minTypability=0.0, maxHeterozygosity=1.0) {
    stats <- barcode.computeBarcodeStats (barcodeData, bySnp)
    result <- barcode.filterByTypabilityStats (barcodeData, stats, bySnp, minTypability)
    result <- barcode.filterByHeterozygosityStats (result, stats, bySnp, maxHeterozygosity)
    result
}

barcode.filterByTypabilityStats <- function(barcodeData, stats, bySnp, minTypability) {
    selectIdx <- which(stats$typable >= minTypability)
    result <- if (bySnp) barcodeData[,selectIdx] else barcodeData[selectIdx,]
    result
}

barcode.filterByHeterozygosityStats <- function(barcodeData, stats, bySnp, maxHeterozygosity) {
    selectIdx <- which(stats$het <= maxHeterozygosity)
    result <- if (bySnp) barcodeData[,selectIdx] else barcodeData[selectIdx,]
    result
}

barcode.computeBarcodeStats <- function (bcodes, bySnp) {
    margin <- if (bySnp) 2 else 1
    occurrences <- if (bySnp) nrow(bcodes) else ncol(bcodes)
    rnames <- if (bySnp) colnames(bcodes) else rownames(bcodes)

    missCounts <- apply(bcodes, margin, function(x) length(which(x=="X")))
    hetCounts <- apply(bcodes, margin, function(x) length(which(x=="N")))

    missing <- (missCounts / occurrences)
    typabile <- 1 - missing;
    valid   <- (occurrences - missCounts)
    het     <- (hetCounts / valid)

    stats <- data.frame(missing, typabile, het)
    colnames(stats) <- c("missing", "typable", "het")
    rownames(stats) <- rnames
    stats
}

barcode.writeBarcodeStats <- function(ctx, datasetName, barcodeData, suffix="") {
    if (nchar(suffix) > 0) {
        suffix <- paste(".", suffix, sep="")
    }
    stats <- barcode.computeBarcodeStats (barcodeData, bySnp=TRUE)
    statsFilename <- paste("stats-snps", suffix, ".tab",  sep="")
    statsFile  <- getContextCacheFile(ctx, datasetName, "barcode", statsFilename)
    writeLabelledData (stats, "Snp", statsFile)

    stats <- barcode.computeBarcodeStats (barcodeData, bySnp=FALSE)
    statsFilename  <- paste("stats-samples", suffix, ".tab",  sep="")
    statsFile  <- getContextCacheFile(ctx, datasetName, "barcode", statsFilename)
    writeLabelledData (stats, "Sample", statsFile)
}

###############################################################################
# Sequence Output
################################################################################
barcode.writeFasta <- function(allelesData, genosFilename) {
  strData <- data.frame(lapply(allelesData, as.character), stringsAsFactors=FALSE)
  txt <- c()
  sampleNames <- rownames(allelesData)
  for (mIdx in 1:length(sampleNames)) {
    header <- paste (">",sampleNames[mIdx],sep='')
    seq <- paste (strData[mIdx,], collapse='')
    txt <- c(txt, header, seq)
  }
  writeLines(txt, paste0(genosFilename,".fasta"))
}
