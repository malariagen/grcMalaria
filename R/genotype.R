###############################################################################
# Declarations: sample data and metadata
################################################################################
geno.getGenoFile <- function (ctx, datasetName) {
    genoFile  <- getContextCacheFile(ctx, datasetName, "genotypes", "sampleGenotypes.tab")
    genoFile
}

geno.initialize <- function (ctx, datasetName, loadFromCache=TRUE, store=TRUE) {
    config <- ctx$config
    dataset <- ctx[[datasetName]]

    genoDataFile <- geno.getGenoFile (ctx, datasetName)
    
    if (loadFromCache & file.exists(genoDataFile)) {
        genoData <- readSampleData (genoDataFile)
        geno.setDatasetGenotypes (ctx, datasetName, genoData, store=FALSE)
    } else {
        genoData <- geno.convertAllelesToGenos (ctx, dataset$barcodes)
        geno.setDatasetGenotypes (ctx, datasetName, genoData, store=store)
    }
}

geno.setDatasetGenotypes <- function (ctx, datasetName, genos, store=TRUE) {
    config <- ctx$config
    dataset <- ctx[[datasetName]]
    dataset$genos <- genos
    if (store) {
        genoDataFile <- geno.getGenoFile (ctx, datasetName)
        writeSampleData(genos, genoDataFile)
    }
}

#
# Convert a dataframe of alleles, to a dataframe of genotypes (values between 0 and 1)
# Missing values are NA and het calls are 0.5
#
geno.convertAllelesToGenos <- function(ctx, alleleData) {
    alleleMeta <- barcode.getMetadata (ctx, alleleData)

    # Change to 0-1 genotypes
    sampleCount <- nrow(alleleData)
    snpCount    <- ncol(alleleData)
  
    genoData <- data.frame(matrix(nrow=sampleCount, ncol=0))
    for (pIdx in 1:snpCount) {
        ref <- alleleMeta$Ref[pIdx]
        nref <- alleleMeta$Nonref[pIdx]
        
        posNt <- alleleData[,pIdx]
        badIdx <- which(!(posNt %in% c(ref, nref, "N","X")))
        if (length(badIdx)>0) {
            bad <- badIdx[1]
            stop (paste("Bad allele found in SNP #",pIdx," in sample ",rownames(alleleData)[bad],": found ",posNt[bad],sep=""))
        }

        genos <- rep(1.0, sampleCount)
        genos[which(posNt == ref)] <- 0.0
        genos[which(posNt == "X")] <- NA
        genos[which(posNt == "N")] <- 0.5
        genoData <- cbind(genoData, genos)
    }
  
    # Name the rows and columns
    colnames(genoData) <- colnames(alleleData) 
    rownames(genoData) <- rownames(alleleData)
    genoData
}
