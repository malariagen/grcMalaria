
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
StringMatrix extractBarcodeAlleles (StringVector barcodes, StringVector snpNames) {
    unsigned int nsnps = snpNames.length();
    unsigned int nsamples = barcodes.length();
    StringVector sampleNames = barcodes.names();
    StringMatrix alleles(nsamples, nsnps);
    Rcpp::Rcout << "Created Matrix: " << nsamples<< " x " << nsnps << std::endl;
    rownames(alleles) = sampleNames;
    colnames(alleles) = snpNames;
    for (unsigned int s1 = 0; s1 < nsamples; s1++) {
        for (unsigned int s2 = 0; s2 < nsnps; s2++) {
             std::string s;
             s.push_back(barcodes[s1][s2]);
             alleles(s1, s2) = s;
        }
    }
    return (alleles);
}
