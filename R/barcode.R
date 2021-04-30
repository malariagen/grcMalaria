###############################################################################
# Declarations: sample data and metadata
################################################################################
folder.data.barcode <- getDataFolder ("barcodes")

barcodeFname     <- "barcodeAlleles.tab"
barcodeMetaFname <- "barcodeMeta.tab"
barcodeSeqFname  <- "barcodeSeqs.fasta"

###############################################################################
# Barcode retrieval, validating and filtering
################################################################################
initializeBarcodes <- function (context, loadFromCache=TRUE, outFolder=folder.data.barcode) {
    barcodeDataFile <- analysis.getDataFile(context, folder.data.barcode, barcodeFname)
    barcodeMetaFile <- analysis.getDataFile(context, folder.data.barcode, barcodeMetaFname)
    
    if (loadFromCache & file.exists(barcodeDataFile)) {
        barcodeMeta <- read.table(barcodeMetaFile, as.is=TRUE, header=TRUE, sep="\t")
        barcodeData <- readSampleData (barcodeDataFile)
        context <- setContextBarcodes (context, barcodeData, barcodeMeta, store=FALSE)
    } else {
        # Read and initialize the barcode metadata
        # It must be a tab-separated file with one row per SNP, and at least four fields named "Order", "SnpName", "Ref", "Nonref"
        barcodeMetaFile <- paste (barcodeMeta.folder, barcodeMeta.file, sep="/")
        barcodeMeta <- read.table(barcodeMetaFile, as.is=TRUE, stringsAsFactors=FALSE, quote="", header=TRUE, sep="\t")
        barcodeMeta <- barcodeMeta[,c("Order", "SnpName", "Ref", "Nonref")]
        colnames(barcodeMeta) <- c("Order", "SnpName", "Ref", "Nonref")
        barcodeMeta <- barcodeMeta[order(barcodeMeta$Order),]
        rownames(barcodeMeta) <- barcodeMeta$SnpName
    
        # Get barcode alleles, and discard samples that have too much missingness
        barcodeData <- getAllelesFromBarcodes (context$meta, meta.colName.barcode, barcodeMeta)
        
        print(paste("Barcode alleles - Samples:", nrow(barcodeData), "x SNPs:", ncol(barcodeData)))
        print("Validating barcodes")
        validateBarcodeAlleles (barcodeData, barcodeMeta)
        
        # Filter the barcodes by typability, trying to throw away as little as possible
        writeBarcodeStats (barcodeData, outFolder=folder.data.barcode, prefix=context$name, suffix="noFiltering")
        print(paste("Filtering barcodes by typability (samples:", barcode.minTypability.sample, ", SNPs:", barcode.minTypability.snp, ")",sep=""))
        #
        # 1) remove all samples with <0.5 typability, so they affect less the removal of SNPs
	filteredData <- filterByTypability (barcodeData, bySnp=FALSE, minTypability=0.5)
	#
	# 2) Refine further, using the thresholds specified
	filteredData <- filterByTypability (filteredData, bySnp=TRUE,  minTypability=barcode.minTypability.snp)
	filteredData <- filterByTypability (filteredData, bySnp=FALSE, minTypability=barcode.minTypability.sample)
	barcodeData <- filteredData
        writeBarcodeStats (barcodeData, outFolder=folder.data.barcode, prefix=context$name, 
                           suffix=paste("filtered-snps_", barcode.minTypability.snp, "-samples_", barcode.minTypability.sample, sep=""))
        print(paste("Barcode alleles after filtering - Samples:", nrow(barcodeData), "x SNPs:", ncol(barcodeData)))
        #
        barcodeMeta <- trimBarcodeMeta (barcodeMeta, barcodeData)
        context <- setContextBarcodes (context, barcodeData, barcodeMeta, store=TRUE, outFolder=outFolder)
    }

    # Report missingness 
    totalCalls <- nrow(barcodeData) * ncol(barcodeData)
    callCounts <- table(unlist(barcodeData))
    missingCalls <- callCounts["X"]
    missing <- missingCalls/totalCalls
    het <- callCounts["N"]/(totalCalls-missingCalls)
    print(paste("Missing:", missing, "- Het:", het))

    context
}

setContextBarcodes <- function (context, newBarcodes, newBarcodeMeta, store=TRUE, outFolder=folder.data.barcode) {
    context$barcodes    <- newBarcodes
    context$barcodeMeta <- newBarcodeMeta
    if (store) {
        barcodeDataFile <- analysis.getDataFile(context, outFolder, barcodeFname)
        barcodeMetaFile <- analysis.getDataFile(context, outFolder, barcodeMetaFname)
        barcodeSeqFile  <- analysis.getDataFile(context, outFolder, barcodeSeqFname)
        writeSampleData(newBarcodes, barcodeDataFile)
        write.table(newBarcodeMeta, file=barcodeMetaFile, sep="\t", quote=FALSE, row.names=FALSE)
        writeFasta (newBarcodes, barcodeSeqFile)
    }
    context
}

#
# Verify all alleles extracted are valid
#
validateBarcodeAlleles <- function (barcodeData, barcodeMeta) {
    snpCount <- nrow(barcodeMeta)				#; print(snpCount)
    if (ncol(barcodeData) != snpCount) {
        stop (paste("Number of barcode SNPs in the metadata (",snpCount,") does not match the length of the barcodes (",ncol(barcodeData),")",sep=""))
    }
    for (sIdx in 1:snpCount) {
        snpMeta <- barcodeMeta[sIdx,]				#; print(snpMeta)
        alleles <- c(snpMeta$Ref,snpMeta$Nonref,"X","N")	#; print(alleles)
        calls <- barcodeData[,sIdx]				#; print(calls)
        badIdx <- which(!(calls %in% alleles))			#; print(badIdx)
        if (length(badIdx) > 0) {
            bad <- badIdx[1]
            stop (paste("Bad allele found in SNP #",sIdx," in sample ",rownames(barcodeData)[bad],": found ",calls[bad,sIdx],sep=""))
        }
    }
}
#
# Convert barcodes into a dataframe of alleles, filtering both samples and barcode SNPs by typability
#
getAllelesFromBarcodes <- function(sampleMetadata, barcodeColumnName, barcodeMeta) {
    barcodes <- as.character(sampleMetadata[,barcodeColumnName])
    names(barcodes) <- rownames(sampleMetadata)

    # Eliminate all samples without barcode
    barcodes <- barcodes[which(barcodes != "-")]
    
    # Split the barcodes into constituent alleles
    alleleMat <- extractBarcodeAlleles (barcodes, rownames(barcodeMeta))
    alleleData <- data.frame(alleleMat)
    rownames(alleleData) = names(barcodes);
    alleleData
}
cppFunction('
    StringMatrix extractBarcodeAlleles (StringVector barcodes, StringVector snpNames) {
        unsigned int nsnps = snpNames.length();
        unsigned int nsamples = barcodes.length();
        StringVector sampleNames = barcodes.names();
        //CharacterVector sampleNames = barcodes.attr("names");
        StringMatrix alleles(nsamples, nsnps);
        for (unsigned int s1 = 0; s1 < nsamples; s1++) {
            for (unsigned int s2 = 0; s2 < nsnps; s2++) {
                std::string s;
                s.push_back(barcodes[s1][s2]);
                alleles(s1, s2) = s;
            }
        }
        rownames(alleles) = sampleNames;
	colnames(alleles) = snpNames;
        return (alleles);
    }
')
#
# Adjust barcode SNP columns after filtering (e.g. by typability)
#
trimBarcodeMeta <- function(barcodeMeta, barcodeData) {
    snpIds <- colnames(barcodeData)
    barcodeMeta <- barcodeMeta[snpIds,]
}
#
###############################################################################
# Barcode Sample/SNP filtering
################################################################################
filterByTypability <- function(barcodeData, bySnp=FALSE, minTypability=0.75) {
    stats <- computeBarcodeStats (barcodeData, bySnp)
    result <- filterByTypabilityStats (barcodeData, stats, bySnp, minTypability)
    result
}

filterByHeterozygosity <- function(barcodeData, bySnp=FALSE, maxHeterozygosity=0.0) {
    stats <- computeBarcodeStats (barcodeData, bySnp)
    result <- filterByHeterozygosityStats (barcodeData, stats, bySnp, maxHeterozygosity)
    result
}

filterByStats <- function(barcodeData, bySnp=FALSE, minTypability=0.0, maxHeterozygosity=1.0) {
    stats <- computeBarcodeStats (barcodeData, bySnp)
    result <- filterByTypabilityStats (barcodeData, stats, bySnp, minTypability)
    result <- filterByHeterozygosityStats (result, stats, bySnp, maxHeterozygosity)
    result
}

filterByTypabilityStats <- function(barcodeData, stats, bySnp, minTypability) {
    selectIdx <- which(stats$typable >= minTypability)
    result <- if (bySnp) barcodeData[,selectIdx] else barcodeData[selectIdx,]
    result
}

filterByHeterozygosityStats <- function(barcodeData, stats, bySnp, maxHeterozygosity) {
    selectIdx <- which(stats$het <= maxHeterozygosity)
    result <- if (bySnp) barcodeData[,selectIdx] else barcodeData[selectIdx,]
    result
}

computeBarcodeStats <- function (bcodes, bySnp) {
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

writeBarcodeStats <- function(barcodeData, outFolder=folder.data.barcode, prefix="", suffix="") {
    if (nchar(prefix) > 0) {
        prefix <- paste(prefix, ".", sep="")
    }
    if (nchar(suffix) > 0) {
        suffix <- paste(".", suffix, sep="")
    }
    stats <- computeBarcodeStats (barcodeData, bySnp=TRUE)    
    statsFile  <- paste(outFolder, "/", prefix, "stats-snps", suffix, ".tab",  sep="")
    writeLabelledData (stats, "Snp", statsFile)

    stats <- computeBarcodeStats (barcodeData, bySnp=FALSE)
    statsFile  <- paste(outFolder, "/", prefix, "stats-samples", suffix, ".tab",  sep="")
    writeLabelledData (stats, "Sample", statsFile)
}

###############################################################################
# Sequence Output
################################################################################
writeFasta <- function(allelesData, genosFilename) {
  strData <- data.frame(lapply(allelesData, as.character), stringsAsFactors=FALSE)
  txt <- c()
  sampleNames <- rownames(allelesData)
  for (mIdx in 1:length(sampleNames)) {
    header <- paste (">",sampleNames[mIdx],sep='')
    seq <- paste (strData[mIdx,], collapse='')
    txt <- c(txt, header, seq)
  }
  writeLines(txt, genosFilename)
}
