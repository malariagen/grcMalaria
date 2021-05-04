# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Test method
#'
#' @export
testRcppFunctions <- function() {
  sampleNames <- c("S1","S2","S3")
  snps <- c("BC1","BC2","BC3","BC4")
  bcodes <- c("AGCT","AXAT","GGAG")
  genos <- matrix(data=c(0,0,1,0, 0,NA,0,0, 1,0,0,1), nrow=3, ncol=4, byrow=TRUE)
  bcodeMeta  <- data.frame (Order=c(1,2,3,4), SnpName=snps, 
                            Ref=c("A","G","A","T"),	Nonref=c("G","T","C","G"))
  names(bcodes) <- sampleNames
  rownames(genos) <- sampleNames
  colnames(genos) <- snps
  
  bcodeMat <- extractBarcodeAlleles (bcodes, snps)
  print (bcodeMat)

  
  distMat <- computeDistances (genos)
  print (distMat)

  
  bcodesData <- data.frame (bcodeMat)
  distData   <- data.frame (distMat)
  mat <- imputeBarcodesCpp (bcodesData, bcodeMeta, distData)
  print (mat)
}
