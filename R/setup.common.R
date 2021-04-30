folder.root <<- getwd()
print(paste("Root folder:", folder.root))

# Load all libraries, quietly
print("Loading libraries")
PKGs <- c("Rcpp", 
          "readxl", "openxlsx", 
          "ape", "pegas", "pcaMethods", #"umap",
          "igraph",
          "GADMTools", "OpenStreetMap", 
          "RColorBrewer", 
          "ggmap", "ggplot2", "ggrepel", "ggforce")
if (exists("showLibraryLoadingMessages") && (showLibraryLoadingMessages==TRUE)) {
    lapply(PKGs, library, character.only=TRUE)
} else {
    loaded <- suppressMessages(lapply(PKGs, library, character.only=TRUE, quietly=TRUE))
}
#
source(paste(folder.code, "util.R", sep="/"), local=FALSE)
#
folder.config  <- getSubFolder (folder.root, "config", create=FALSE)
folder.data    <- getSubFolder (folder.root, "data")
folder.results <- getSubFolder (folder.root, "out")
print(paste("Result folder:", folder.results))

print("Loading configuration")
setwd(folder.config)
source("config.R", local=FALSE)

print("Loading Tools")
setwd(folder.code)
source("metadata.R", local=FALSE)
source("barcode.R", local=FALSE)
source("impute.R", local=FALSE)
source("genotype.R", local=FALSE)
source("distance.R", local=FALSE)
source("analysis.R", local=FALSE)
#
# Options
#
setwd(folder.root)
options(scipen=10)
options(stringsAsFactors=FALSE)
options(error=dump.frames)

