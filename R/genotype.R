###############################################################################
# Declarations: sample data and metadata
################################################################################
geno.getGenoFile <- function (ctx, datasetName) {
    genoFile  <- getContextCacheFile(ctx, datasetName, "genotypes", "sampleGenotypes")
    genoFile
}

geno.initialize <- function (ctx, datasetName, loadFromCache=TRUE, store=TRUE) {
    config <- ctx$config
    dataset <- ctx[[datasetName]]

    genoDataFile <- geno.getGenoFile (ctx, datasetName)
    
    if (loadFromCache & rdaFileExists(genoDataFile)) {
        genoData <- readRdaSampleData (genoDataFile)
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
        writeRdaSampleData(genos, genoDataFile)
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
  
    validAlleles <- c("A","C","G","T","X","N")
    genoData <- data.frame(matrix(nrow=sampleCount, ncol=0))
    for (pIdx in 1:snpCount) {
        ref <- alleleMeta$Ref[pIdx]
        
        posNt <- alleleData[,pIdx]
        badIdx <- which(!(posNt %in% validAlleles))
        if (length(badIdx) > 0) {
            for (bIdx in 1:length(badIdx)) {
                badAllele <- posNt[bIdx]
                badSampleIdx <- badIdx[bIdx]
                badSample <- rownames(alleleData)[badSampleIdx]
                cat (paste0("Error: Unexpected allele found at SNP #",pIdx," in sample ",badSample,": found ",badAllele), fill=TRUE)
            }
            stop ("Errors found in barcoding alleles - exiting.")
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
