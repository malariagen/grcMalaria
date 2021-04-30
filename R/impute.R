###############################################################################
# Declarations: sample data and metadata
################################################################################

#
###################################################################
# Imputed Data, with filled-in missingness in filtered barcodes 
###################################################################
#
impute.createContext <- function (filteredContext, loadFromCache=TRUE) {
    print("Initializing Imputed Dataset")
    imputedContext <- list(name="imputed")
    imputedContext <- meta.setContextSampleMeta (imputedContext, filteredContext$meta)

    impBarcodeDataFile <- analysis.getDataFile(imputedContext, folder.data.barcode, barcodeFname)
    impBarcodeMeta <- filteredContext$barcodeMeta
    if (loadFromCache & file.exists(impBarcodeDataFile)) {
        impBarcodeData <- readSampleData (impBarcodeDataFile)
        print(paste("Loaded", imputedContext$name, "barcodes - Samples:", nrow(impBarcodeData), "x SNPs:", ncol(impBarcodeData)))
    } else {
        impBarcodeData <- impute.imputeBarcodes (filteredContext$barcodes, impBarcodeMeta, filteredContext$distance)
    }
    imputedContext <- setContextBarcodes (imputedContext, impBarcodeData, impBarcodeMeta)

    # Get the genotypes and distance matrix from imputed data
    imputedContext <- geno.initialize(imputedContext)
    imputedContext <- distance.initialize(imputedContext)
    imputedContext
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

cppFunction('
    StringMatrix imputeBarcodesCpp (DataFrame barcodeData, DataFrame barcodeMeta, DataFrame distData) {
        int nsamples = barcodeData.nrow();
        int nsnps    = barcodeData.ncol();
        CharacterVector sampleNames = barcodeData.attr("row.names");
        CharacterVector snpNames = barcodeData.attr("names");

	Function order_R("order");
	
        CharacterVector refAlleles = barcodeMeta["Ref"];
	CharacterVector nrefAlleles = barcodeMeta["Nonref"];
	
        int nTopAlleles = 100;
        StringVector topAlleles(nTopAlleles);
        NumericVector topScores(nTopAlleles);
        
        StringMatrix imp(nsamples, nsnps);
        for (int sIdx = 0; sIdx < nsamples; sIdx++) {
            String sName = sampleNames[sIdx];
            for (int vIdx = 0; vIdx < nsnps; vIdx++) {
                StringVector posCalls = barcodeData[vIdx];
                String call = posCalls[sIdx];
                if ((call == "X") || (call == "N")) {
                    // We need to impute at this position.
                    // Get the calls for the most similar samples
                    NumericVector dists = distData[sName];
                    StringVector alleles = barcodeData[vIdx];
		    IntegerVector perm = order_R(dists);
		    int topIdx = 0;
                    for (int i = 0; ((topIdx < nTopAlleles) && (i < nsamples)); i++) {
		         int idx = (perm[i]-1);
		         if (idx == sIdx) {
		             continue;
		         }
		         String a = alleles[idx];
		         if ((a == "X") || (a == "N")) {
		             continue;
		         }
		         topAlleles[topIdx] = a;
		         topScores[topIdx] = 1.0 - dists[idx];
		         topIdx++;
                    }
                    //Rcout << topAlleles << "\\n";
                    //Rcout << topScores << "\\n";
                    // Now score ref and nonref alleles
                    String refAllele = refAlleles[vIdx];
                    String nrefAllele = nrefAlleles[vIdx];
                    double refScore = 0.0;
                    double nrefScore = 0.0;
                    for (int i = 0; i < nTopAlleles; i++) {
                         String allele = topAlleles[i];
		         if(allele == refAllele) {
		             refScore += topScores[i];
		         } else if(allele == nrefAllele){
		             nrefScore += topScores[i];
		         } else {
		             String snpName = snpNames[vIdx];
		             stop("Unknown allele for sample %s at position %s: %s", sName.get_cstring(), snpName.get_cstring(), allele.get_cstring());
		         }
                    }
                    String impAllele = (refScore > nrefScore) ? refAllele : nrefAllele;
                    imp(sIdx, vIdx) = impAllele;
                } else {
                    imp(sIdx,vIdx) = call;
                }
            }
        }
        return (imp);
    }
')

#impBarcodes <- initializeImputation(loadFromCache=FALSE)


##########################################################
exampleRearrange <- function() {
    d <- data.frame(c("Napoli","Milano","Roma"),c(3,1,2))
    colnames(d) <- c("City","Num")
    exampleRearrangeCpp(d)
}
cppFunction('
    StringVector exampleRearrangeCpp (DataFrame df) {
        int nrec = df.nrows();
        IntegerVector num = df["Num"];     // Second column
        StringVector cities = df["City"];  // First column
        Function srt("order");   
        IntegerVector perm = srt(num);
        StringVector result(nrec);
        for (int i = 0; i < nrec; i++) {
             int idx = (perm[i]-1);
             String s = cities[idx];
             result[i] = s;
        }
        return (result);
    }
')

cppFunction('
    String test1 (String s1, String s2) {
        if(s1 == s2) {
            return ("Yes");
        }
        return ("No");
    }
')
#test1 ("Hello", "there")
#test1 ("there", "there")

exampleRowNames <- function() {
    d <- data.frame(c("Napoli","Milano","Roma"),c(3,1,2))
    colnames(d) <- c("City","Num")
    rownames(d) <- c("A","B","C")
    exampleRowNamesCpp(d)
}
cppFunction('
    StringVector exampleRowNamesCpp (DataFrame df) {
        int nrec = df.nrows();
        CharacterVector sNames = df.attr("row.names");
        StringVector snpNames = colnames(df);
        for (int sIdx = 0; sIdx < nrec; sIdx++) {
            String sName = sNames[sIdx];
            Rcout << sName.get_cstring() << "\\n";
        }
        return (sNames);
    }
')
#exampleRowNames()


##########################################################


