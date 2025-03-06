#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <Rcpp.h>
using namespace Rcpp;

#define GENO_HOM 1
#define GENO_HET 2
#define GENO_HET_NO_PROPS 3
#define GENO_MISS 9

// [[Rcpp::export]]
std::vector<std::string> getColumnAlleles(unsigned int cIdx, List colGenotypes) {
    std::vector<unsigned int> sGenotypes = colGenotypes["sampleGenotypes"];
    unsigned int sCount = sGenotypes.size();
    List sAlleles = colGenotypes["sampleAlleles"];
    std::unordered_set<std::string> uniqueStrings;  // Set to store unique strings
    for (unsigned int sIdx = 0; sIdx < sCount; sIdx++) {
        //Rcpp::Rcout << cIdx << " - " << sIdx << "\n";
        switch (sGenotypes[sIdx]) {
        case GENO_HOM:
        case GENO_HET:
            List alleleProps = sAlleles[sIdx];
            //CharacterVector ch = wrap(alleleProps.attributeNames());
            //Rcout << ch << " - " << alleleProps.length() << "\n";
            StringVector alleleNames = alleleProps.names();
            unsigned int aCount = alleleNames.length();
            for (unsigned int aIdx = 0; aIdx < aCount; aIdx++) {
                uniqueStrings.insert(std::string(alleleNames[aIdx]));
		    }
            break;
        }
    }
    // Convert set to vector for returning
    std::vector<std::string> result(uniqueStrings.begin(), uniqueStrings.end());
    return result;
}

// [[Rcpp::export]]
std::vector<std::vector<double>> getAllelePropTable (List colGenotypes, std::vector<std::string> colAlleles) {
	//
	// Make an index table to facilitate index lookups
    std::unordered_map<std::string, int> indexMap;
    for (size_t i = 0; i < colAlleles.size(); ++i) {
        indexMap[colAlleles[i]] = i;
    }
    //
    std::vector<unsigned int> sGenotypes = colGenotypes["sampleGenotypes"];
    unsigned int sCount = sGenotypes.size();
    unsigned int aCount = colAlleles.size();
    List sAlleles = colGenotypes["sampleAlleles"];
    //
    // Create a 2D vector (sCount rows and aCount columns), initialized to 0.0
    std::vector<std::vector<double>> propTable(sCount, std::vector<double>(aCount, 0.0));
    for (unsigned int sIdx = 0; sIdx < sCount; sIdx++) {
        switch (sGenotypes[sIdx]) {
        case GENO_HOM:
        case GENO_HET:
            List alleleProps = sAlleles[sIdx];
            StringVector alleleNames = alleleProps.names();
            unsigned int aCount = alleleNames.length();
            for (unsigned int aIdx = 0; aIdx < aCount; aIdx++) {
				double aProp = alleleProps[aIdx];
				std::string aName = std::string(alleleNames[aIdx]);
                int alleleIndex = indexMap[aName];
				propTable[sIdx][alleleIndex] = aProp;
			}
            break;
        }
    }
    return (propTable);
}

// [[Rcpp::export]]
std::vector<int> getHomAlleleIndexes (List colGenotypes, std::vector<std::string> colAlleles) {
	//
	// Make an index table to facilitate index lookups
    std::unordered_map<std::string, int> indexMap;
    for (size_t i = 0; i < colAlleles.size(); ++i) {
        indexMap[colAlleles[i]] = i;
    }
    //
    std::vector<unsigned int> sGenotypes = colGenotypes["sampleGenotypes"];
    unsigned int sCount = sGenotypes.size();
    List sAlleles = colGenotypes["sampleAlleles"];
    std::vector<int> homIndexes (sCount, -1);
    for (unsigned int sIdx = 0; sIdx < sCount; sIdx++) {
        switch (sGenotypes[sIdx]) {
        case GENO_HOM:
            List alleleProps = sAlleles[sIdx];
            StringVector alleleNames = alleleProps.names();
            std::string allele (alleleNames[0]);
            homIndexes[sIdx] = indexMap[allele];
            break;
        }
    }
    return (homIndexes);
}

// [[Rcpp::export]]
NumericMatrix dist_calculateDistanceMatrix (List genotypeData) {
    //Rcpp::Rcout << "Estimating Distance Matrix\n";
    StringVector samples = genotypeData["samples"];
    unsigned int sCount = samples.size();
    List colGenoData = genotypeData["columnGenoData"];
    StringVector cols = colGenoData.names();
    unsigned int cCount = cols.size();
    //
    NumericMatrix distMat  = NumericMatrix(sCount, sCount); rownames(distMat)  = samples; colnames(distMat)  = samples;
    NumericMatrix countMat = NumericMatrix(sCount, sCount); rownames(countMat) = samples; colnames(countMat) = samples;
    //
    for (unsigned int cIdx = 0; cIdx < cCount; cIdx++) {
        //Rcpp::Rcout << "Barcoding column: " << cIdx << "\n";
        String colName = cols[cIdx];
        List colGenotypes = colGenoData[colName];
        //for (unsigned int idx = 0; idx < colGenotypes.names().size(); ++idx) { Rcpp::Rcout << colGenoNames[idx] << idx << "\n"; }
        IntegerVector sGenoVec = colGenotypes["sampleGenotypes"];
        std::vector<unsigned int> sGenotypes = as<std::vector<unsigned int>>(sGenoVec);
        List sAlleles = colGenotypes["sampleAlleles"];
        //
        // Get all the alleles present in this column, and construct a hashmap to find their index
        std::vector<std::string> colAlleles = getColumnAlleles(cIdx, colGenotypes);
        std::unordered_map<std::string, int> indexMap;
        for (unsigned int i = 0; i < colAlleles.size(); ++i) {
            indexMap[colAlleles[i]] = i;
        }
        //
        // Create a table with the proportions for each allele/sample combination, to make lookup faster
        std::vector<std::vector<double>> propTable = getAllelePropTable (colGenotypes, colAlleles);
        //
        // For each hom call, get the index of the allele, to facilitate fast computation
        std::vector<int> homIndexes = getHomAlleleIndexes (colGenotypes, colAlleles);
        //
        // And now calculate pairwise distances
        //
        //
        for (unsigned int s1Idx = 0; s1Idx < (sCount-1); s1Idx++) {
            //Rcpp::Rcout << "Sample1: " << s1Idx << "\n";
            unsigned int g1 = sGenotypes[s1Idx];
            if (g1 == GENO_MISS) continue;
            //String s1Name = samples[s1Idx];
            //
            for (unsigned int s2Idx = (s1Idx+1); s2Idx < sCount; s2Idx++) {
                //Rcpp::Rcout << "Sample2: " << s2Idx << "\n";
                unsigned int g2 = sGenotypes[s2Idx];
                if (g2 == GENO_MISS) continue;
                //String s2Name = samples[s2Idx];
                //
                double similarity = 0;
                //
                if ((g1 == GENO_HET_NO_PROPS) || (g2 == GENO_HET_NO_PROPS)) {
					// At least one of the sample is heterozygous but we don't know the proportion, assume a distance of 0.5
					// TODO - We might do this in a more sophisticated way later on: e.g. getAllelePropTable() could assign 0.5
					// proportion to each of the two most frequent alleles at this position, and treat as regular HET.
					// But for now this will do, since we still do not have much in terms of allele proportions.
                    similarity = 0.5;

                } else if ((g1==GENO_HOM) && (g2==GENO_HOM)) {
					// There is only one allele in each StringVector, so get the first element of each
                    if (homIndexes[s1Idx]==homIndexes[s2Idx]) {
                        similarity = 1;
                    }

				} else if (g1==GENO_HOM) {
					// There is only one allele in sample 1, so look up the proportion in sample 2 for the same allele
					int a1Index = homIndexes[s1Idx];
					similarity = propTable[s2Idx][a1Index];

				} else if (g2==GENO_HOM) {
					// There is only one allele in sample 2, so look up the proportion in sample 1 for the same allele
					int a2Index = homIndexes[s2Idx];
					similarity = propTable[s1Idx][a2Index];

                } else {
					//
                    // Het+Het
                    for (unsigned int i = 0; i < colAlleles.size(); i++) {
						if (propTable[s1Idx][i] > 0) {
							similarity += (propTable[s1Idx][i] * propTable[s2Idx][i]);
						}
				    }
                }
				double distance = 1 - similarity;
			    //Rcpp::Rcout << "Distance: " << distance << "\n";
                distMat(s1Idx,s2Idx) += distance;
                distMat(s2Idx,s1Idx) += distance;
                countMat(s1Idx,s2Idx)++;
                countMat(s2Idx,s1Idx)++;
			}
		}
	}
	for (unsigned int s1Idx = 0; s1Idx < (sCount-1); s1Idx++) {
        for (unsigned int s2Idx = (s1Idx+1); s2Idx < sCount; s2Idx++) {
			distMat(s1Idx,s2Idx) /= countMat(s1Idx,s2Idx);
            distMat(s2Idx,s1Idx) /= countMat(s2Idx,s1Idx);
		}
	}
    return (distMat);
}
