
#include <Rcpp.h>
#include <cstdio>
#include <cstdlib>

using namespace Rcpp;

// [[Rcpp::export]]
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
    rownames(imp) = sampleNames;
    colnames(imp) = snpNames;
    //Rcpp::Rcout << "Created Matrix: " << nsamples<< " x " << nsnps << std::endl;
    for (int sIdx = 0; sIdx < nsamples; sIdx++) {
        Rcpp::String sName = sampleNames[sIdx];
        //Rcpp::Rcout << "Sample: " << sName.get_cstring() << std::endl;
        for (int vIdx = 0; vIdx < nsnps; vIdx++) {
            StringVector posCalls = barcodeData[vIdx];
            String call = posCalls[sIdx];
            //Rcpp::Rcout << "Snp #" << vIdx << ": " << call.get_cstring() << std::endl;
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
                        // Rcout << "Found : " << a ;
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
                for (int i = 0; i < topIdx; i++) {
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
