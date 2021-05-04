###############################################################################
# Declarations: sample data and metadata
################################################################################

#
###################################################################
# Imputed Data, with filled-in missingness in filtered barcodes 
###################################################################
#
impute.createImputedDataset <- function (ctx, loadFromCache=TRUE) {
    print("Initializing Imputed Dataset")
    ctx$imputed <- list(name="imputed")
    filteredDs <- ctx$filtered
    config     <- ctx$config
    
    ctx <- meta.setDatasetMeta (ctx, "imputed", filteredDs$meta)
    impBarcodeDataFile <- barcode.getBarcodeDataFile (ctx, "imputed")
    if (loadFromCache & file.exists(impBarcodeDataFile)) {
        impBarcodeData <- readSampleData (impBarcodeDataFile)
        print(paste("Loaded imputed barcodes - Samples:", nrow(impBarcodeData), "x SNPs:", ncol(impBarcodeData)))
    } else {
        barcodeData <- filteredDs$barcodes
        barcodeMeta <- barcode.getMetadata(ctx, barcodeData)
        impBarcodeData <- impute.imputeBarcodes (barcodeData, barcodeMeta, filteredDs$distance)
    }
    ctx <- setDatasetBarcodes (ctx, "imputed", impBarcodeData)

    # Get the genotypes and distance matrix from imputed data
    ctx <- geno.initialize(ctx, "imputed")
    ctx <- distance.initialize(ctx, "imputed")
    ctx
}
#
impute.imputeBarcodes <- function (barcodeData, barcodeMeta, distData) {
    sampleNames <-  rownames(barcodeData)
    distData <- distData[sampleNames, sampleNames]
    impMat <- imputeBarcodesCpp (barcodeData, barcodeMeta, distData)
    impBarcodeData <- data.frame(impMat)
    rownames(impBarcodeData) <- sampleNames
    colnames(impBarcodeData) <- colnames(barcodeData)
    #print(barcodeData[1:10,1:30])
    #print(impBarcodeData[1:10,1:30])
    impBarcodeData
}

