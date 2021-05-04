## code to prepare `DATASET` dataset goes here

# ###########################################################
# Geographical setup data
# ###########################################################
geoFile <- "data-raw/GeoData.xlsx"

print("Loading GADM province names")
gadmProvTable <- data.frame(readxl::read_excel(geoFile, sheet="GADM-provinces", col_types="text"))

print("Loading Country Data")
countryDataTable <- data.frame(readxl::read_excel(geoFile, sheet="ISO-3166", col_types="text"))
rownames(countryDataTable) <- countryDataTable$iso2

map.geoTables <- list(gadmProv=gadmProvTable, countries=countryDataTable)

# ###########################################################
# Barcode SNP structure
# ###########################################################
# Read and initialize the barcode metadata
# It must be a tab-separated file with one row per SNP, and at least four fields named "Order", "SnpName", "Ref", "Nonref"
readBarcodeMeta <- function(filename) {
    file <- paste("data-raw/barcodeMeta", filename, sep="/")
    meta <- utils::read.table(file, as.is=TRUE, header=TRUE, sep="\t")
    meta <- meta[,c("Order", "SnpName", "Ref", "Nonref")]
    meta <- meta[order(meta$Order),]
    rownames(meta) <- meta$SnpName
    meta
}

barcode.metadata.Pf <- list (
    "v1.0" = readBarcodeMeta("BarcodingSnps-Pf-1.0.tab")
)
barcode.metadata.Pv <- list (
    "v1.0" = readBarcodeMeta("BarcodingSnps-Pv-1.0.tab")
)

barcode.metadata <- list (
    Pf=barcode.metadata.Pf,
    Pv=barcode.metadata.Pv
)

# ###########################################################
# File Storage
# ###########################################################

usethis::use_data(map.geoTables, barcode.metadata, internal=TRUE, overwrite=TRUE)
