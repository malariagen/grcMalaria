###################################################################
# Filtered Data, with high-missingness barcodes removed
###################################################################
#
filter.createFilteredDataset <- function (ctx, loadFromCache=TRUE) {
    print("Initializing Filtered Dataset")
    filteredDs <- context.createContextDataset (ctx, "filtered")
    #ctx$filtered <- list(name="filtered")
    #filteredDs   <- ctx$filtered
    
    unfilteredDs <- ctx$unfiltered
    config       <- ctx$config

    filteredMetaFile        <- meta.getMetaDataFile(ctx, "filtered")
    filteredBarcodeFile     <- barcode.getBarcodeDataFile(ctx, "filtered")

    if (loadFromCache & file.exists(filteredMetaFile) & file.exists(filteredBarcodeFile)) {
        meta <- readSampleData (filteredMetaFile)		#; print(colnames(meta))
        meta.setDatasetMeta (ctx, "filtered", meta, store=FALSE)
        barcodeData <- readSampleData (filteredBarcodeFile)
        barcode.setDatasetBarcodes (ctx, "filtered", barcodeData, store=FALSE)
        print(paste("Loaded filtered barcodes - Samples:", nrow(barcodeData), "x SNPs:", ncol(barcodeData)))
    } else {
        meta.setDatasetMeta (ctx, "filtered", ctx$unfiltered$meta, store=FALSE)
        barcode.initializeBarcodes (ctx, "filtered")
        # Trim the metadata to cover the barcodes selected
        sampleNames   <- rownames(ctx$filtered$barcodes)
        filteredMeta  <- filter.filterSampleData (ctx$unfiltered$meta, sampleNames)
        meta.setDatasetMeta (ctx, "filtered", filteredMeta)
    }

    # Get the genotypes, distance matrix and execute the pop structure analysis
    geno.initialize(ctx, "filtered")
    distance.initialize(ctx, "filtered")
}
#
filter.filterSampleData <- function (data, sNames) {
    data <- data[sNames,]
    data
}
#
#filter.filterMatrix <- function (matrixData, sampleNames) {
#    matrixData <- matrixData[sampleNames,sampleNames]
#    matrixData
#}
