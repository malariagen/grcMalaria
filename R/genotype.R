###############################################################################
# Declarations: sample data and metadata
################################################################################
folder.data.geno <- getDataFolder ("genotypes")

genoFname      <- "sampleGenotypes.tab"
genoMetaFname  <- "sampleGenotypeMeta.tab"

geno.initialize <- function (context, loadFromCache=TRUE, store=TRUE, outFolder=folder.data.geno) {
    genoDataFile <- analysis.getDataFile(context, folder.data.geno, genoFname)
    genoMetaFile <- analysis.getDataFile(context, folder.data.geno, genoMetaFname)
    if (loadFromCache & file.exists(genoDataFile)) {
        genoMeta <- read.table(genoMetaFile, as.is=TRUE, header=TRUE, sep="\t")
        genoMeta <- genoMeta[,c("Order", "SnpName", "Ref", "Nonref")]
	colnames(genoMeta) <- c("Order", "SnpName", "Ref", "Nonref")
        genoData <- readSampleData (genoDataFile)
        context <- geno.setContextGenotypes (context, genoData, genoMeta, store=FALSE)
    } else {
        genoMeta <- context$barcodeMeta
        genoData <- geno.convertAllelesToGenos (context$barcodes, genoMeta)
        context <- geno.setContextGenotypes (context, genoData, genoMeta, store=store, outFolder=outFolder)
    }
    context
}

geno.setContextGenotypes <- function (context, newGenos, newGenoMeta, store=TRUE, outFolder=folder.data.geno) {
    context$genos    <- newGenos
    context$genoMeta <- newGenoMeta
    if (store) {
        genoDataFile <- analysis.getDataFile(context, outFolder, genoFname)
        genoMetaFile <- analysis.getDataFile(context, outFolder, genoMetaFname)
        writeSampleData(newGenos, genoDataFile)
        write.table(newGenoMeta, file=genoMetaFile, sep="\t", quote=FALSE, row.names=FALSE)
    }
    context
}

#
# Convert a dataframe of alleles, to a dataframe of genotypes (values between 0 and 1)
# Missing values are NA and het calls are 0.5
#
geno.convertAllelesToGenos <- function(alleleData, alleleMeta) {
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
