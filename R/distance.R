###############################################################################
# Declarations: sample data and metadata
################################################################################
folder.data.distance <- getDataFolder ("distance")

distFname      <- "sampleDistance.tab"

distance.initialize <- function (context, loadFromCache=TRUE, store=TRUE, outFolder=folder.data.distance) {
    distDataFile <- analysis.getDataFile(context, folder.data.distance, distFname)
    if (loadFromCache & file.exists(distDataFile)) {
        distData <- readMatrix(distDataFile)
        context <- distance.setContextDistance(context, distData, store=FALSE)
        print(paste("Loaded", context$name, "distance matrix - Samples:", nrow(distData)))
    } else {
        genoData <- context$genos
        print(paste("Computing pairwise distances for",nrow(genoData),"samples using",ncol(genoData),"SNPs"))
        distData <- computeDistances (as.matrix(genoData))
        context <- distance.setContextDistance(context, distData, store=store, outFolder=outFolder)
    }
    context
}

distance.setContextDistance <- function (context, newDistance, store=TRUE, outFolder=folder.data.distance) {
    context$distance <- newDistance
    if (store) {
        distDataFile <- analysis.getDataFile(context, outFolder, distFname)
        writeMatrix(newDistance, distDataFile)
    }
    context
}

cppFunction('
    NumericMatrix computeDistances (NumericMatrix genos) {
        unsigned int nsnps    = genos.ncol();
        unsigned int nsamples = genos.nrow();
        CharacterVector sNames = rownames(genos);

        NumericMatrix dist(nsamples, nsamples);
        rownames(dist) = sNames;
        colnames(dist) = sNames;
        
        for (int s1 = 0; s1 < (nsamples-1); s1++) {
            for (int s2 = (s1+1); s2 < nsamples; s2++) {
                double d = 0.0;
                double cnt = 0.0;
                for (int v = 0; v < nsnps; v++) {
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
')

