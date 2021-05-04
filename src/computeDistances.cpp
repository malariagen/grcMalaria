
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix computeDistances (NumericMatrix genos) {
  unsigned int nsnps    = genos.ncol();
  unsigned int nsamples = genos.nrow();
  CharacterVector sNames = rownames(genos);

  NumericMatrix dist(nsamples, nsamples);
  rownames(dist) = sNames;
  colnames(dist) = sNames;

  for (unsigned int s1 = 0; s1 < (nsamples-1); s1++) {
    for (unsigned int s2 = (s1+1); s2 < nsamples; s2++) {
      double d = 0.0;
      double cnt = 0.0;
      for (unsigned int v = 0; v < nsnps; v++) {
        double g1 = genos(s1, v);
        if (NumericVector::is_na(g1)) continue;
        double g2 = genos(s2, v);
        if (NumericVector::is_na(g2)) continue;
        d += (g1*(1.0-g2) + g2*(1.0-g1));
        cnt++;
      }
      d = (cnt > 0) ? (d / cnt) : NA_REAL;
      dist(s1,s2) = d;
      dist(s2,s1) = d;
    }
  }
  return (dist);
}

