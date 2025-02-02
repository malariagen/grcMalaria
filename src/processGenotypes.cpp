#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <Rcpp.h>

using namespace Rcpp;

#define GENO_HOM 1
#define GENO_HET 2
#define GENO_HET_NO_PROPS 3
#define GENO_MISS 9

/*

*/
// [[Rcpp::export]]
StringVector geno_splitString (String input, char delimiter) {
    // Read tokens from a string stream separated by the delimiter
    std::istringstream stream(input);
    std::string token;
    StringVector split;
    while (std::getline(stream, token, delimiter)) {
        split.push_back(token);
    }
    return (split);
}
/*
# Test calls, run in R
testSplitString <- function () {
    geno_splitString("Hello there", " ")
    geno_splitString("Hello there", "e")
}
*/

/*
 Tests whether a given character (of C++ type char) is present within a string.
*/
// [[Rcpp::export]]
bool geno_containsChar(std::string input, char sub) {
    std::size_t subPos = input.find(sub);
    return (subPos != std::string::npos);
}

/*
# Test calls, run in R
testContainsString <- function () {
    geno_containsChar ("CVI-T", '-')
    geno_containsChar ("CVIET", '-')
    geno_containsChar ("CVIE-", '-')
    geno_containsChar ("CV---", '-')
    geno_containsChar ("Hello there", ' ')
    geno_containsChar ("<NA>", '2')
    geno_containsChar ("Hello there", '2')
    geno_containsChar ("Hello there", 'e')
    geno_containsChar ("Hello there", 'z')
}
*/

/*
 Tests whether a genotype sting represents a missing value (a dash, a string containing a dash, or the string "<NA>")
*/
// [[Rcpp::export]]
bool geno_isMissingGenotype(String genoStr) {
	std::string input = std::string(genoStr);
	return ((genoStr == "-")
	     || (genoStr == "<NA>")
	     || (geno_containsChar(genoStr, '-')));
}

/*
# Test calls, run in R
testMissingGenotype <- function () {
    geno_isMissingGenotype ("CVI-T")
    geno_isMissingGenotype ("CVIET")
    geno_isMissingGenotype ("CVIE-")
    geno_isMissingGenotype ("<NA>")
    geno_isMissingGenotype ("Hello there")
}
*/

/*
 Parse allele string of the form <allele>[:<proportion>]. Proportion is set to 0 if unspecified.
 Args: alleleStr - An allele string of the form <allele>[:<proportion>]
 Return: allele: a one-element list with a name (the allele) and a numeric (the proportion)
*/
// [[Rcpp::export]]
List geno_parseAllele (String alleleStr) {
    StringVector p1 = geno_splitString (alleleStr, ':');
    unsigned int cnt = 0;
    if (p1.length() > 1) {
        cnt = stoi(std::string(p1[1]));
    }
    String n = p1[0];
    List allele = List::create(Named(n)=cnt);
    return (allele);
}

/*
 Parses a comma-separated sequences of allele strings, representing a sample genotype at a given locus.
 If the input string is "-" (missing), returns NULL.

 Args:   genotypeStr - a string containing a comma-separated sequence of allele strings of the form <allele>[:<proportion>]
 Return: Alleles     - a list in which each element has a name (the allele) and a numeric (the proportion).
                       The sum of proportions for all alleles is 1.0.
*/
// [[Rcpp::export]]
List geno_parseGenotypeString (String genotypeStr) {
    List alleles = List::create();
    if (geno_isMissingGenotype(genotypeStr)) {
        return (NULL);
    }
    StringVector p1 = geno_splitString (genotypeStr, ',');
    unsigned int p1cnt = p1.length();
    if (p1cnt == 1) {
        alleles = geno_parseAllele (genotypeStr);
        alleles[0] = 1.0;
    } else {
        bool noProps = false;
        StringVector  alleleNames;
        NumericVector props;
        for (unsigned int i = 0; i < p1cnt; i++) {
            String pString = p1[i];
            List a = geno_parseAllele (pString);
            double prop = a[0];
            if (prop == 0.0) {
                noProps = true;
            }
            CharacterVector an = a.names();
            alleleNames.push_back(an[0]);
            props.push_back(prop);
        }
        double alleleCount = (double)alleleNames.length();
        if (noProps) {
            double prop = 1.0 / alleleCount;
            for (unsigned int i = 0; i < alleleCount; i++) {
                props[i] = prop;
            }
        } else {
            double total = 0.0;
            for (int i = 0; i < alleleCount; i++) {
                total += props[i];
            }
            for (int i = 0; i < alleleCount; i++) {
                props[i] /= total;
            }
        }
        for (int i = 0; i < alleleCount; i++) {
             String n = alleleNames[i];
             double p = props[i];
             alleles[n] = p;
        }
    }
    return (alleles);
}
/*
# Test calls, run in R
testParseGenotypeString <- function (gs) {
    print(paste0("Testing genotype string: [",gs,"]"))
    print(geno_parseGenotypeString(gs))
}
testParseGenotypeStrings <- function () {
    testParseGenotypeString("-")
    testParseGenotypeString("A")
    testParseGenotypeString("A,G")
    testParseGenotypeString("A:10,T:20")
    testParseGenotypeString("A:10,T:20,C:120")
    testParseGenotypeString("C:10,T")
}
testParseGenotypeStrings()
*/

/*
 Parses a single column, each sample genotype encoded as a comma-separated sequence of
 allele strings of the form <allele>[:<proportion>]

 Args:   genoStrs - A vector of string, each containing a sample genotype encoded as a comma-separated sequence
                    of allele strings of the form <allele>[:<proportion>]; the element name is the sample ID.
         samples - A vector of string, each containing the sample ID for the corresponding genotype string.
 Return: VariantGenotypes: a list containing the following elements:
         - genotypes: a vector of integers indicating for each sample the genotyping poutcome (hom, het, missing etc.)
         - miss: a numeric indicating the proportion of samples with a missing call
         - het: a numeric indicating the proportion of non-missing sample calls that are heterozygous (multiple alleles)
         - hom: a numeric indicating the proportion of non-missing sample calls that are homozygous (one alleles only)
         - <for each sample with non-missing genotype>: a SampleAlleles list containing the alleles and their proportions
                                                        for a given variant/sample combination.
*/
// [[Rcpp::export]]
List geno_parseColumnGenotypes (StringVector genoStrs, StringVector samples) {
    unsigned int sCount = samples.length();
    IntegerVector genotypes (sCount, 0);
    List sampleAlleles(sCount,R_NilValue);
    //
    unsigned int hom  = 0;
    unsigned int miss = 0;
    for (unsigned int sIdx = 0; sIdx < sCount; sIdx++) {
        unsigned int sGenotype;
        List sAlleles = NULL;
        String sample = samples[sIdx];
        String sGenotypeStr = genoStrs[sIdx];
        if (geno_isMissingGenotype(sGenotypeStr)) {
            sGenotype = GENO_MISS;
            miss++;
        } else if (((int)std::string(sGenotypeStr).find("*")) >= 0) {
            sGenotype = GENO_HET_NO_PROPS;
        } else {
            //
            // Parse the genotype string to extract an array of named values, corresponding to
            // the relative proportion of the allele; the names are the names of the alleles
            //
            sAlleles = geno_parseGenotypeString (sGenotypeStr);
            unsigned int alleleCnt = sAlleles.length();
            if (alleleCnt == 1) {
                sGenotype = GENO_HOM;
                hom++;
            } else {
                sGenotype = GENO_HET;
            }
        }
        genotypes[sIdx] = sGenotype;
        if ((sGenotype == GENO_HOM) || (sGenotype == GENO_HET)) {
            sampleAlleles[sIdx] = sAlleles;
        }
    }
    genotypes.attr("names") = samples;
    sampleAlleles.attr("names") = samples;
    //
    double het = sCount - (miss + hom);
    //
    List colGenotypes = List::create(
        Named("samples")         = samples,
        Named("sampleGenotypes") = genotypes,
        Named("sampleAlleles")   = sampleAlleles,
        Named("totalCount")      = sCount,
        Named("missingCount")    = miss,
        Named("homCount")        = hom,
        Named("hetCount")        = het
    );
    return (colGenotypes);
}

/*
# Test calls, run in R
testParseGenotypeStrings <- function () {
    genos <- c("-", "A", "CVIET", "CVI-T", "A,G", "A:10,T:20", "-", "A:10,T:20,C:120", "C,T", "C:10,T")
    samples <- paste0("S", 1:10)
    print(geno_parseColumnGenotypes (genos, samples))
}
testParseGenotypeStrings()
*/

/*
 Works out the proportions of missing and heterozygous calls for all samples, acrosss all barcoding variants
 Args: samples - The vector of sample IDs
       genoCols - The name of the variant columns
       columnGenoData - a list containing the parsed genotyping data for all variants and samples
 Return: a list contianing two vectors (each vector has names, which are the sample IDs):
       missing: a vector of the proportions of calls that are missing, for each sample
       het: a vector of the proportions of non-missing calls calls that are het, for each sample
*/
// [[Rcpp::export]]
List geno_estimateSampleProperties (StringVector samples, StringVector genoCols, List columnGenoData) {
    unsigned int sCount = samples.length();
    unsigned int colCount = genoCols.length();

    IntegerVector missingCounts(sCount, 0);
    IntegerVector homCounts(sCount, 0);
    IntegerVector hetCounts(sCount, 0);

    for (unsigned int colIdx = 0; colIdx < colCount; colIdx++) {
        String genoCol = genoCols[colIdx];
        List columnGenos = columnGenoData[genoCol];
        IntegerVector sGenotypes = columnGenos["sampleGenotypes"];
        for (unsigned int sIdx = 0; sIdx < sCount; sIdx++) {
             unsigned int sGeno = sGenotypes[sIdx];
             if (sGeno == GENO_MISS) {
                 missingCounts[sIdx]++;
             } else if (sGeno == GENO_HOM) {
                 homCounts[sIdx]++;
             } else {
                 hetCounts[sIdx]++;
			 }
         }
    }
    missingCounts.attr("names") = samples;
    homCounts.attr("names") = samples;
    hetCounts.attr("names") = samples;
    //
    List props = List::create(
        Named("columnCount")   = colCount,
        Named("missingCounts") = missingCounts,
        Named("homCounts")     = homCounts,
        Named("hetCounts")     = hetCounts
    );
    return (props);
}

/*
 Works out the proportions of missing and heterozygous calls for all variant columns, acrosss all samples
 Args: genoCols - A data frame of GRC Features relevant to barcoding variants (one row per feature)
       columnGenoData: a  list containing the parsed genotyping data for all variants and samples
 Return: a list containing two vectors (each vector has names, which are the sample IDs):
         missing: a vector of the proportions of calls that are missing, for each variant
         het: a vector of the proportions of non-missing calls calls that are het, for each variant
*/
// [[Rcpp::export]]
List geno_estimateColumnProperties (StringVector genoCols, List columnGenoData) {
    unsigned int colCount = genoCols.length();

    IntegerVector missingCounts(colCount, 0);
    IntegerVector homCounts(colCount, 0);
    IntegerVector hetCounts(colCount, 0);

    for (unsigned int colIdx = 0; colIdx < colCount; colIdx++) {
        String genoCol = genoCols[colIdx];
        List columnGenos = columnGenoData[genoCol];
        missingCounts[colIdx] = columnGenos["missingCount"];
        homCounts[colIdx]     = columnGenos["homCount"];
        hetCounts[colIdx]     = columnGenos["hetCount"];
    }
    missingCounts.attr("names") = genoCols;
    homCounts.attr("names") = genoCols;
    hetCounts.attr("names") = genoCols;
    unsigned int sampleCount = missingCounts[0] + homCounts[0] + hetCounts[0];

    List props = List::create(
        Named("sampleCount")   = sampleCount,
        Named("missingCounts") = missingCounts,
        Named("homCounts")     = homCounts,
        Named("hetCounts")     = hetCounts
    );
    return (props);
}

/*
 Process all data in the allele-based columns referenced by features.
 The names of these columns are listed in grc$config$alleleColumns.

 Args:   grcData - A data frame containing the GRC data, read from file
         genoCols - A vector with the names of the columns to be genotyped
 Return: GenotypeData: a list containing the following elements:
         - <for each barcoding feature>: a VariantGenotypeData list containing the genotype results for that variant.
         - samples: the vector of sample IDs
         - sampleMissing: a vector of numerics containing the proportion of variants that were missing for each sample
         - sampleHet: a vector of numerics containing the proportion of non-missing variant calls that were het for each sample
         Each element of the list is indexed by the variant feature name
*/
// [[Rcpp::export]]
List geno_processGenotypes (DataFrame grcData, StringVector genoCols) {
    unsigned int colCount = genoCols.length();
    StringVector samples = grcData["SampleId"];
    unsigned int sCount = samples.length();
    //
    // Parse the text-encoded variant genotype data, one variant at a time.
    //
    List columnGenoData = List::create();
    for (unsigned int vIdx = 0; vIdx < colCount; vIdx++) {
        String genoCol = genoCols[vIdx];
        StringVector genoStrs = as<StringVector>(grcData[genoCol]);
        List colGenos = geno_parseColumnGenotypes(genoStrs, samples);
        columnGenoData[genoCol] = colGenos;
    }
    //
    // Estimate missingness and heterozyous call proportion
    //
    List sProp = geno_estimateSampleProperties (samples, genoCols, columnGenoData);
    List cProp = geno_estimateColumnProperties (genoCols, columnGenoData);
    //
    List genotypeData = List::create(
        Named("samples")        = samples,
        Named("columns")        = genoCols,
        Named("columnGenoData") = columnGenoData,
        //
        Named("columnCount")         = sProp["columnCount"],
        Named("sampleMissingCounts") = sProp["missingCounts"],
        Named("sampleHomCounts")     = sProp["homCounts"],
        Named("sampleHetCounts")     = sProp["hetCounts"],
        //
        Named("sampleCount")         = cProp["sampleCount"],
        Named("columnMissingCounts") = cProp["missingCounts"],
        Named("columnHomCounts")     = cProp["homCounts"],
        Named("columnHetCounts")     = cProp["hetCounts"]
    );
    return (genotypeData);
}

