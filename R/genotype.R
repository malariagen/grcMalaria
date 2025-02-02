#
# This module contains the code that translates genotype strings in the GRC data file into genotype data internally used by the package.
# The code works on columns of data=; please note that the same column may be used by more than one feature (e.g. a multi-allele variant and a specific mutation)
#
GENO.HOM  <- 1
GENO.HET  <- 2
GENO.HET_NO_PROPS <- 3
GENO.MISS <- 9
#
# Process all data in the allele-based columns referenced by features. The names of these columns are listed in config$alleleColumns.
#
# Args:   grcData - A data frame containing the GRC data, read from file
#         genoCols - A vector with the names of the columns to be genotyped
# Return: GenotypeData: a list containing the following elements:
#         - <for each barcoding feature>: a VariantGenotypeData list containing the genotype results for that variant.
#         - samples: the vector of sample IDs
#         - sampleMissing: a vector of numerics containing the proportion of variants that were missing for each sample
#         - sampleHet: a vector of numerics containing the proportion of non-missing variant calls that were het for each sample
#         Each element of the list is indexed by the variant feature name
#
genotype.processGenotypes <- function (grcData, genoCols) {
    missingCols <- which(!(genoCols %in% colnames(grcData)))
    if (length(missingCols) > 0) {
        stop(paste0("Missing columns of genotype data: [", paste0(genoCols[missingCols], collapse=","), "]"))
    }
    genotypeData <- geno_processGenotypes(grcData, genoCols)  # NOW USES RCPP CODE
    genotypeData
}
#
# ###########################################################################################
# Genotype Data Filtering
# ###########################################################################################
#
# Build a GenotypeData list that contains a subset of samples of a larger GenotypeData list
#
genotype.filterGenotypeDataBySample <- function (genotypeData, selSamples) {
    #
    # Check that all the samples to be selected are in the source barcode set
    #
    sCount <- length(selSamples)			#; print(paste("DBG filterGenotypeDataBySample new sample count:", sCount))
    srcSamples <- genotypeData$samples			#; print(paste("DBG filterGenotypeDataBySample prev sample count:", length(srcSamples)))
    selIdx <- which(selSamples %in% srcSamples)
    if (length(selIdx) != sCount) {
        xIdx <- which(!(selSamples %in% srcSamples))
        print(selSamples[xIdx])
        stop("One or more selected samples are not in the original sample list")
    }
    #
    columnGenoData <- genotypeData$columnGenoData
    colNames <- names(columnGenoData)
    colCount <- length(colNames)			#; print(paste("DBG filterGenotypeDataBySample col count:", colCount))
    newColumnGenoData <- list()
    #
    for (cIdx in 1:colCount) {
        colName <- colNames[cIdx]			#; print(paste("DBG filterGenotypeDataBySample", cIdx, colName))
        colGenotypes <- columnGenoData[[colName]]
        #
        sampleGenotypes <- colGenotypes$sampleGenotypes
        newSampleGenotypes <- sampleGenotypes[selSamples]
        #
        missingCount <- length(which(newSampleGenotypes==GENO.MISS))
        homCount     <- length(which(newSampleGenotypes==GENO.HOM))
        hetCount     <- sCount - (missingCount + homCount)
        #
        sampleAlleles <- colGenotypes$sampleAlleles
        keep <- selSamples[which(selSamples %in% names(sampleAlleles))]
        newSampleAlleles <- sampleAlleles[keep]
        #
        newColGenotypes <- list(samples=selSamples,
                                sampleGenotypes=newSampleGenotypes,
                                sampleAlleles=newSampleAlleles,
                                totalCount=sCount,
                                missingCount=missingCount,
                                homCount=homCount,
                                hetCount=hetCount)
        newColumnGenoData[[colName]] <- newColGenotypes
    }
    #
    # Estimate sample missingness and heterozyous call proportion
    #
    sProp <- geno_estimateSampleProperties (selSamples, colNames, newColumnGenoData)
    cProp <- geno_estimateColumnProperties (colNames, newColumnGenoData)
    #
    newGenotypeData <- list(samples=selSamples, 
                            columns=colNames,
                            columnGenoData=newColumnGenoData, 
                            #
                            columnCount=sProp$columnCount, 
                            sampleMissingCounts=sProp$missingCounts, 
                            sampleHomCounts=sProp$homCounts,
                            sampleHetCounts=sProp$hetCounts,
                            #
                            sampleCount=cProp$sampleCount, 
                            columnMissingCounts=cProp$missingCounts, 
                            columnHomCounts=cProp$homCounts, 
                            columnHetCounts=cProp$hetCounts
                            )
    newGenotypeData
}
#
# Build a GenotypeData list that contains a subset of columns of a larger GenotypeData list
#
genotype.filterGenotypeDataByColumn <- function (genotypeData, selColumns) {
    #
    # Check that all the columns to be selected are in the source barcode set
    #
    columnGenoData <- genotypeData$columnGenoData
    colNames <- names(columnGenoData)
    selIdx <- which(colNames %in% selColumns)
    if (length(selIdx) != length(selColumns)) {
        stop("One or more selected columns are not in the original list")
    }
    #
    newColumnGenoData<- columnGenoData[selColumns]
    #
    # Estimate sample missingness and heterozyous call proportion
    #
    sProp <- geno_estimateSampleProperties (genotypeData$samples, selColumns, newColumnGenoData)
    cProp <- geno_estimateColumnProperties (selColumns, newColumnGenoData)
    #
    newGenotypeData <- list(samples=genotypeData$samples, 
                            columns=selColumns,
                            columnGenoData=newColumnGenoData,
                            #
                            columnCount=sProp$columnCount, 
                            sampleMissingCounts=sProp$missingCounts, 
                            sampleHomCounts=sProp$homCounts,
                            sampleHetCounts=sProp$hetCounts,
                            #
                            sampleCount=cProp$sampleCount, 
                            columnMissingCounts=cProp$missingCounts, 
                            columnHomCounts=cProp$homCounts, 
                            columnHetCounts=cProp$hetCounts
                            )
    newGenotypeData
}
