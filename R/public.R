#############################################################
#
#' grcMalaria: Processing and analyzing SpotMalaria Genetic Report Cards
#'
#'\if{html}{\figure{logo.jpg}{options: align='right' alt='logo' width='120'}}
#' Various functions for producing analysis results such as maps and reports from the Genetic Report Cards datasets
#' delivered by SpotMalaria projects, such as GenRe-Mekong.
#' Typical outputs include drug resistance prevalence maps, population genetic analyses such as PCA and NJ trees,
#' and analyses of genetic barcodes such as barcode sharing networks, etc.
#' The data frame returned contains a data table whose values and structure may have been manipulated after being read from the file,
#' and may have been updated to the latest version fo the format.
#'
#' @section Author:
#' Olivo Miotto \email{olivo@tropmedres.ac}
#'
#' @seealso Useful links:
#' \itemize{
#'  \item user guide {\url{https://olivo-miotto-genre.notion.site/grcMalaria-R-package-User-Guide-dfaac7ea30b5430b8522e8b531fd2e4b}}
#'  \item source code {\url{https://github.com/malariagen/grcMalaria}}
#'  \item reporting bugs {\url{https://github.com/malariagen/grcMalaria/issues}}
#' }
#'
#' @docType package
#' @name grcMalaria
NULL

#############################################################
#
#' Load a Genetic Report Cards data file
#'
#' Load a GRC data file, ready for analysis.
#' The data frame returned contains a data table whose values and structure may have been manipulated after being read from the file,
#' and may have been updated to the latest version fo the format.
#'
#' @param file The path to the GRC data Microsoft Excel file (use forward slashes in the path)
#' @param sheet The name of the datasheet or tab name within the Excel file containing the GRC data
#' @param species The species being analyzed ("Pf" (i.e. Plasmodium falciparum) or "Pv" (i.e. Plasmodium vivax))
#' @param version The version number of the GRC data file format. This is in the documentation when you receive or download the file.
#'
#' @return A list containing a data frame with the data ready to be analyzed, plus some configuration metadata
#' @export
#'
#' @seealso Information on the Genetic Report Card (GRC) structure input file {\url{https://olivo-miotto-genre.notion.site/Genetic-Report-Card-GRC-Structure-40c4f6cd8a8144eab5f1a4151ce9127a}}
#'
#' @examples
#' # Load data file
#' # Change the path to where your file is located before running the code
#' Data <- loadGrc("D:/.../name_file.xlsx",
#'                sheet = "GenRe-Mekong",
#'                species = "Pf", version = "1.0")
#
loadGrc <- function (file, sheet="GenRe-Mekong", species="Pf", version="1.0") {
  grcData <- grcData.load (file, sheet, species, version)
  grcData
}

#############################################################
#
#' Combine two Genetic Report Cards datasets
#'
#' Merge two GRC datasets which have been loaded from data files. The new dataset is added to the source dataset.
#' If the new dataset contains new columns, they can be either ignored, or the same columns are added to the source dataset,
#' with null values. If the new dataset contains sample data that is already in the source dataset, they can be either ignored,
#' or the source data rows can be overwritten.
#'
#' @param srcGrc The source GRC dataset (loaded by loadGrc())
#' @param newGrc The new GRC dataset to be added(loaded by loadGrc())
#' @param overwrite Boolean: if TRUE, sample rows in newGrc will replace rows in srcGrc that have the same sampleId. If FALSE, they will be disregarded.
#' @param extendColumns Boolean: if TRUE, new columns in newGrc will also be added to srcGrc, filled with <NA> values.  If FALSE, they will be disregarded.
#'
#' @return A list containing a data frame with the merged data ready to be analyzed, plus some configuration metadata
#' @export
#'
#' @examples
#' # Load spreadsheet 1
#' Sheet1 <- loadGrc("D:/.../sheet1.xlsx",
#'                  sheet="GRC", species="Pf", version="1.0")
#' # Load spreadsheet 2
#' Sheet2 <- loadGrc("D:/.../sheet2.xlsx",
#'                  sheet="GRC", species="Pf", version="1.0")
#' # Merge datasets
#' mergeData <- mergeGrc(Sheet1, Sheet2)

#
mergeGrc <- function (srcGrc, newGrc, overwrite=FALSE, extendColumns=FALSE) {
    grcData <- grcData.merge (srcGrc, newGrc, overwrite, extendColumns)
    grcData
}

#############################################################
#
#' Initialize a Genetics Report Cards dataset
#'
#' Initialize a GRC dataset (obtained from function loadGrc()) so that it is ready for analysis.
#' It performs barcode imputation, computes genetic distances, and initializes data that will subsequently be used in analyses.
#’ It creates a folder where output of future analyses will be saved. It will take a while to create a context object which will be used for the subsequent analysis tasks (can take as long as 5 mins to run).
#'
#' @param grcData The data obtained from reading the GRC Excel data.
#' @param dir The folder where the outputs from this and subsequent analyses will be stored.
#' @param minSnpTypability The minimum proportion of non-missing samples for a barcode position to be retained in analysis. The default is 0.8. Meaning, SNPs where the allele is missing in more than 20 percent of samples are removed.
#' @param minSampleTypability The minimum proportion of non-missing positions for a sample to be retained in analysis. The default is 0.75. Meaning, samples where more than 25 percent of the remaining SNPs are missing are removed.
#'
#' @return A analysis context object, which is a list that contains all the data for analysis, which will be passed to subsequent analysis tasks.
#' @export
#'
#' @examples
#' #Change the path to where you want output file to be
#' ctx <- initializeContext(Data,
#'                          dir="D:/...")
#
initializeContext <- function (grcData, dir=".",
                               minSnpTypability=0.8,
                               minSampleTypability=0.75) {
  options(scipen=10)
  options(stringsAsFactors=FALSE)

  config <- setup.getConfig (grcData, dir, minSnpTypability, minSampleTypability)
  ctx <- analysis.createContext (grcData$data, config)
  ctx
}

#############################################################
#
#' Select a set of samples using metadata
#'
#' Selects a set of samples for analysis, based on the provided metadata values.
#' A given analysis context can contain multiple sample sets. These can be labelled with different names.
#' The selection criteria are specified as a list containing a sequence of lists with two elements each:
#' "field" which is the column name to be checked, and "values" which is an array of possible values that can be matched for a sample to be selected.
#' A selected sample must match all the criteria.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSetName The name of the sample set, to be used to identify it when calling analysis tasks. This will create a folder with the same name e.g. “Laos”, or "Test170522".
#' @param select  The criteria for sample selection. Selected samples will match all specified criteria. The values must be a list. "field" corresponds to a column name in the datafile. "values" correspond to an array of possible values that can be matched for a sample to be selected. Use "," to add additional parameters inside a quotation mark " ".
#'
#' @return The analysis context, augmented with the new sample set
#' @export
#'
#' @examples
#' # Create a SampleSet named "Laos", based on "Timepoint" and "Study"
#' ctx <- selectSampleSet(ctx, sampleSet="Laos", select=list(
#' list(field="TimePoint", values=c("D00H00","-")),
#' list(field="Study", values="1208-PF-LA-CMPE-GENRE") ))
#'
#' # Create a SampleSet named "Test170522", containing only isolates from 2018, and 2019 from Vietnam
#' ctx <- selectSampleSet(ctx, sampleSet="Test170522", select=list(
#' list(field="Year", values=c("2018","2019")),
#' list(field="Country", values="VN") ))
#'
#' # Create a SampleSet named "C580Y", containing only isolates with the C580Y genotype
#' ctx <- selectSampleSet(ctx, sampleSet="C580Y", select=list(
#' list(field="Kelch", values="C580Y") ))
#
selectSampleSet <- function (ctx, sampleSetName, select) {
    ctx <- analysis.selectSampleSet (ctx, sampleSetName, select)
    ctx
}

#############################################################
#
#' Map of Sample Counts
#'
#' Creates maps showing the numbers of samples collected at different aggregation levels.
#' For each aggregation unit, we place a marker on the map, coloured according to the administrative
#' division it is in.
#' Additionally data files are created. Created .tab files can be opened with for example Excel.
#' Maps and data are located in .../output/out/sampleSetName/map-sampleCount/.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet().
#' @param timePeriods A separately given list of time period object for partitioning samples into time-interval plots. Each list specifies a time period, with three parameters: name, type and, start. Parameter type can currently only be year. Parameter start must be a date in the format “dd-MMM-yyyy”, and parameter name must be provided. The string passed in the name parameter is added to the end of the file name of the produced files.
#' Currently, the public GRC data only contains years, so only the year in the provided “start” date is used.
#' @param aggregate The administrative level at which we aggregate (Province or District). Separate maps are created for each administrative level.
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District).
#' @param colourBy Shows the aggregation level to be used to color the markers (Country or Province).
#' @param markerSize Allows adjustment of the size of markers on the map. If only one value is passed,
#'                   all markers will be that size; if two values are passed, they will be used as the min
#'                   and max size of the marker, whose size will reflect the number of samples.
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @export
#'
#' @examples
#' #Given lists of time periods of interest, in this case calendar years
#' periods <- list(
#' list(name="2018", type="year", start="1-Jan-2018"),
#' list(name="2019", type="year", start="1-Jan-2019"),
#' list(name="2020", type="year", start="1-Jan-2020")
#' )
#'
#' #Maps showing the numbers of samples collected at Province and District aggregation levels.
#' mapSampleCounts (ctx, sampleSet="Laos", timePeriods=periods,
#' aggregate=c("Province","District"), minAggregateCount=1,
#' markerSize=c(10,40), colourBy="Province", showNames=TRUE,
#' width=15, height=15)
#'
#' #Given lists of time periods of interest, in this case the Southern Laos malaria season
#' periods2 <- list(
#' list(name="2018", type="year", start="1-Sep-2018"),
#' list(name="2019", type="year", start="1-Sep-2019"),
#' list(name="2020", type="year", start="1-Sep-2020")
#' )
#'
#' #Maps showing the numbers of samples collected at Province and District aggregation levels.
#' mapSampleCounts (ctx, sampleSet="Laos", timePeriods=periods2,
#' aggregate=c("Province","District"), minAggregateCount=1,
#' markerSize=c(10,40), colourBy="Province", showNames=TRUE,
#' width=15, height=15)
#
mapSampleCounts <- function (ctx, sampleSet, timePeriods=NULL,
                   aggregate="Province", minAggregateCount=1,
                   markerSize=c(4,40), showNames=TRUE, colourBy="Province",
                   width=15, height=15) {

    if (length(colourBy) > 1) {
        stop ("colourBy parameter can only accet a single value (\"Country\" or \"Province\")")
    }
    colourAggLevel <- map.getAggregationLevelsFromLabels (colourBy)
    if (colourAggLevel > 1) {
        stop ("colourBy parameter can only accet values \"Country\" or \"Province\"")
    }
    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    timeIntervals <- parseTimeIntervals(timePeriods)
    params <- list(
        analysis.timeIntervals=timeIntervals,
        aggegation.levels=aggLevels,
        map.aggregateCountMin=minAggregateCount,
        map.markerColourAggLevel=colourAggLevel,
        map.markerSize=markerSize,
        map.markerNames=showNames,
        plot.size=list(width=width,height=height)
    )
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/sampleCount", params=params)
}
#
#############################################################
#
#' Map prevalence of Drug resistance (predicted phenotype)
#'
#' Creates a map showing the levels of estimated resistance to a particular antimalarial drug. Based on published genetic markers, each sample has a predicted phenotype for different types of antimalarial drugs.
#' If multiple drugs are specified, then different maps will be created for different drugs.
#' The predictions of resistance for any given drug are aggregated at the desired administrative level: Province (level 1), or District (level 2). Separate maps are created for each drug-administrative-combination.
#' For each aggregation unit, we place a marker on the map, coloured according to the level of resistance to the drug, with a label indicating the prevalence.
#' To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' Additionally data files are created. Created .tab files can be opened with for example Excel.
#' Maps and data are located in .../output/out/sampleSetName/map-drug/.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet().
#' @param timePeriods A separately given list of time period object for partitioning samples into time-interval plots. Each list specifies a time period, with three parameters: name, type and, start. Parameter type can currently only be year. Parameter start must be a date in the format “dd-MMM-yyyy”, and parameter name must be provided. The string passed in the name parameter is added to the end of the file name of the produced files.
#' Currently, the public GRC data only contains years, so only the year in the provided “start” date is used.
#' @param drugs The antimalarial drugs for which prevalence of phenotypic resistance will be estimated; "ALL" creates maps for all the drugs for which phenotypic resistance predictions are available,
#' which include "Artemisinin", "Chloroquine", "Piperaquine", "DHA-PPQ" (i.e. Dihydroartemisinin/piperaquine), "Sulfadoxine", "Pyrimethamine", "S-P" (i.e. Sulfadoxine-Pyrimethamine), and "S-P-IPTp" (i.e. Sulfadoxine-Pyrimethamine).
#' To specify a drug put the drug name in between quotation marks e.g. "Artemisinin", or  c(“Artemisinin”, “Chloroquine”, “S-P”) to select several specific drugs.
#' @param aggregate The administrative level at which we aggregate.
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param markerSize Allows adjustment of the size of markers on the map.
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @export
#'
#' @examples
#' #Given lists of time periods of interest, in this case calendar years
#' periods <- list(
#' list(name="2018", type="year", start="1-Jan-2018"),
#' list(name="2019", type="year", start="1-Jan-2019"),
#' list(name="2020", type="year", start="1-Jan-2020")
#' )
#'
#' # Map estimated Artemisinin and Chloroquine phenotypic resistance prevalence for sampleSet "Laos" for both Province and District level for calendar years 2018-2020.
#' mapDrugResistancePrevalence (ctx, sampleSet="Laos", timePeriods=periods,
#' drugs=c("Artemisinin", "Chloroquine"), aggregate=c("Province","District"),
#' minAggregateCount=10, showNames=TRUE, markerSize=16,
#' width=15, height=15)
#'
#' @seealso Useful links:
#' \itemize{
#' \item Details of how drug resistance phenotypes are predicted based on genetic markers {\url{http://ngs.sanger.ac.uk/production/malaria/Resource/29/20200705-GenRe-05-PhenotypeRules-0.39.pdf}}
#' \item Other useful background information about the SpotMalaria platform  {\url{https://www.malariagen.net/resource/29}}
#' }
#
mapDrugResistancePrevalence <- function (ctx, sampleSet, timePeriods=NULL,
                   drugs="ALL",
                   aggregate="Province", minAggregateCount=10,
                   showNames=TRUE, markerSize=16,
                   width=15, height=15) {

    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    timeIntervals <- parseTimeIntervals(timePeriods)
    params <- list(
        analysis.timeIntervals=timeIntervals,
        analysis.measures=drugs,
        aggegation.levels=aggLevels,
        map.aggregateCountMin=minAggregateCount,
        map.markerSize=markerSize,			# Use this for constant marker size
        #map.markerSize=c(4,40),			# Use this for count-proportional marker size
        map.markerNames=showNames,
        plot.size=list(width=width,height=height)
    )
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/drug", params=params)
}
#
#############################################################
#
#' Map prevalence of genetic mutations of drug resistance
#'
#' Creates maps showing the prevalence of genetic mutations types known to be associated with antimalarial drug resistance.
#' If multiple mutation types are specified, then a map will be created for each mutation type.
#' The prevalence of genetic mutation types are aggregated at the desired administrative level: Province (level 1), or District (level 2). Separate maps are created for each mutation-administrative-combination.
#' For each aggregation unit, we place a marker on the map, coloured according to the prevalence of the genetic mutation type, with a label indicating the prevalence.
#' To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' Additionally data files are created. Created .tab files can be opened with for example Excel.
#' Maps and data are located in .../output/out/sampleSetName/map-mutation/.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet().
#' @param timePeriods A separately given list of time period object for partitioning samples into time-interval plots. Each list specifies a time period, with three parameters: name, type and, start. Parameter type can currently only be year. Parameter start must be a date in the format “dd-MMM-yyyy”, and parameter name must be provided. The string passed in the name parameter is added to the end of the file name of the produced files.
#' Currently, the public GRC data only contains years, so only the year in the provided “start” date is used.
#' @param mutations The genetic mutation types for which prevalence will be calculated; "ALL" creates maps for all 42 available genetic mutation types.
#' Available mutation types include: "crt_C72S",	"crt_M74I",	"crt_N75E",	"crt_N75D",	"crt_K76T",	"crt_T93S",	"crt_H97Y",	"crt_H97L",	"crt_I218F",	"crt_A220S",
#' "crt_Q271E",	"crt_N326S",	"crt_N326D",	"crt_T333S",	"crt_I356T",	"crt_I356L",	"crt_R371I",	"dhfr_N51I",	"dhfr_C59R",	"dhfr_C59H",	"dhfr_S108N",
#' "dhfr_S108T",	"dhfr_I164L",	"dhps_S436A",	"dhps_S436F",	"dhps_A437G",	"dhps_K540E",	"dhps_K540N",	"dhps_A581G",	"dhps_A613T",	"dhps_A613S",	"mdr1_N86Y",
#' mdr1_Y184F",	"mdr1_S1034I",	"mdr1_F1226Y",	"mdr1_D1246Y",	"arps10_V127M",	"arps10_D128H",	"arps10_D128Y",	"fd_D193Y",	"mdr2_T484I",	"exo_E415G".
#' To specify a mutation put the mutation type in between quotation marks e.g. "crt_N75E", or  c("arps10_D128Y","fd_D193Y","crt_I218F") to select several genetic mutation types.
#' @param aggregate The administrative level at which we aggregate.
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param markerSize Allows adjustment of the size of markers on the map.
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @export
#'
#' @seealso Useful links:
#' \itemize{
#'  \item user guide for information on functionality of specific genetic mutation types {\url{https://olivo-miotto-genre.notion.site/grcMalaria-R-package-User-Guide-dfaac7ea30b5430b8522e8b531fd2e4b}}
#' }
#'
#' @examples
#' #Given lists of time periods of interest, in this case calendar years
#' periods <- list(
#' list(name="2018", type="year", start="1-Jan-2018"),
#' list(name="2019", type="year", start="1-Jan-2019"),
#' list(name="2020", type="year", start="1-Jan-2020")
#' )
#'
#' # Produce prevalence maps for all available genetic mutation types, for sampleSet "Laos", for District level, and for calendar years 2018-2020.
#' mapMutationPrevalence (ctx, sampleSet="Laos", timePeriods=periods,
#' mutations="ALL", aggregate="District",
#' minAggregateCount=10, showNames=TRUE, markerSize=16,
#' width=15, height=15)
#
mapMutationPrevalence <- function (ctx, sampleSet, timePeriods=NULL,
                   mutations="ALL",
                   aggregate="Province", minAggregateCount=10,
                   showNames=TRUE, markerSize=16,
                   width=15, height=15) {

    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    timeIntervals <- parseTimeIntervals(timePeriods)
    params <- list(
        analysis.timeIntervals=timeIntervals,
        analysis.measures=mutations,
        aggegation.levels=aggLevels,
        map.aggregateCountMin=minAggregateCount,
        map.markerSize=markerSize,			# Use this for constant marker size
        #map.markerSize=c(4,40),			# Use this for count-proportional marker size
        map.markerNames=showNames,
        plot.size=list(width=width,height=height)
    )
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/mutation", params=params)
}

#############################################################
#
#' Map genetic Diversity
#'
#' Creates a map showing the calculated genetic diversity of the collected Plasmodium samples in the selected sampleSet.
#' Genetic diversity can be calculated using four different measures: "maxHaploFreq","haploHet", "meanSnpHet", and "medianDistance".
#' Calculations are based on the 'genetic barcode', which consists of 101 SNPs determined for each sample. These SNPs are distributed across all nuclear chromosomes.
#' SNP positions were chosen based on their ability to differentiate populations and their power to recapitulate genetic distance.
#' If multiple measures are specified, then a map will be created for each measure.
#' The genetic diversity is aggregated at the desired administrative level: Province (level 1), or District (level 2).
#' For each aggregation unit, we place a marker on the map, coloured according to the genetic diversity, with a label indicating the genetic diversity.
#' To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' Additionally data files are created. Created .tab files can be opened with for example Excel.
#' Maps and data are located in .../output/out/sampleSetName/map-diversity/.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet().
#' @param timePeriods A separately given list of time period object for partitioning samples into time-interval plots. Each list specifies a time period, with three parameters: name, type and, start. Parameter type can currently only be year. Parameter start must be a date in the format “dd-MMM-yyyy”, and parameter name must be provided. The string passed in the name parameter is added to the end of the file name of the produced files.
#' Currently, the public GRC data only contains years, so only the year in the provided “start” date is used.
#' @param measures This can be "ALL", or any vector containing one or more of c("maxHaploFreq","haploHet", "meanSnpHet","medianDistance").
#' The method "maxHaploFreq" is a measure of loss of diversity. This measure gives the proportion of samples carrying the most common haplotype (defined as samples with identical barcodes).
#' The output ranges from 0-1, with a low value corresponding to a low proportion of samples with the most common haplotype, and a high value corresponding to a
#' high proportion of samples with the most common haplotype.
#' The method "haploHet" is a measure of heterozygosity based on the complete barcode.
#' This measure gives the probability of two randomly selected samples carrying a different barcode, and is useful to detect large changes in a population structure.
#' The output ranges from 0-1, with a low value corresponding to low diversity (low probability of carrying a different barcode), while a high value corresponds to high diversity (high probability of a different barcode).
#' The method "meanSnpHet" is also known as the expected heterozygosity or gene diversity of a locus.
#' For each barcode SNP, expected heterozygosity is calculated using: HE= n/(n-1) [1-∑ipi^2 ],
#' where n = the number of samples and pi = the allele frequency of the ith SNP in the barcode.
#' The final value shows the mean of heterozygosity across all the loci.
#' This measure is able to detect smaller changes in a population structure, and stable as it uses the mean of all SNPs.
#' The output ranges from 0-1, with a low value corresponding to low diversity (low probability of a different allele), and a high value corresponding to high diversity (high probability of a different allele).
#' The method "medianDistance" is a measure that gives the median genetic distance in a assessed group.
#' The is calculated using pairwise genetic distance as the proportion of SNPs differing between two samples. Then taking the median of these in the group of samples of interest.
#' @param aggregate The administrative level at which we aggregate.
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District).
#' @param markerColours The colour to indicate the level of genetic diversity. Default: "red3", other examples: "forestgreen", "cornflowerblue","darkgoldenrod3".
#' @param markerSize Allows adjustment of the size of markers on the map.
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @export
#'
#' @examples
#' #' #Given lists of time periods of interest, in this case calendar years
#' periods <- list(
#' list(name="2018", type="year", start="1-Jan-2018"),
#' list(name="2019", type="year", start="1-Jan-2019"),
#' list(name="2020", type="year", start="1-Jan-2020")
#' )
#'
#' # Produce genetic diversity maps for all available measures, for sampleSet "Laos", for District level, and for calendar years 2018-2020.
#' mapDiversity (ctx, sampleSet="Laos", timePeriods=periods,
#' measures="ALL", aggregate="District", markerColours="red3",
#' minAggregateCount=10, showNames=TRUE, markerSize=16,
#' width=15, height=15)
#

mapDiversity <- function (ctx, sampleSet, timePeriods=NULL,
                   measures="ALL",
                   aggregate="Province", minAggregateCount=10,
                   showNames=TRUE, markerSize=16, markerColours="red3",
                   width=15, height=15) {

    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    timeIntervals <- parseTimeIntervals(timePeriods)
    params <- list(
        analysis.timeIntervals=timeIntervals,
        analysis.measures=measures,
        aggegation.levels=aggLevels,
        map.aggregateCountMin=minAggregateCount,
        map.markerSize=markerSize,			# Use this for constant marker size
        #map.markerSize=c(4,40),			# Use this for count-proportional marker size
        map.markerNames=showNames,
        map.diversity.markerColours=markerColours,
        plot.size=list(width=width,height=height)
    )
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/diversity", params=params)
}

#############################################################
#
#' Map Connectedness between sites
#'
#' Creates a map of connectedness of haplotypes between different administrative divisions (i.e. sites) in the selected sampleSet.
#' The administrative divisions are compared in a pairwise fashion, within this comparison output is produced that
#' gives either i) the proportion of sample pairs that have at least the assigned percentage of similarity (default: 1), or
#' ii) the proportion of sample pairs with a mean genetic distance of at least the meanDistanceLevels (default: 0.5).
#' Subsequently, maps are produced that show connections between the administrative divisions,
#' in which the thickness of the connection line corresponds to the proportion of either i), or ii).
#' Both similarity and mean genetic distance are calculated based on the 'genetic barcode', which consists of 101 SNPs determined for each sample.
#' These SNPs are distributed across all nuclear chromosomes.
#' SNP positions were chosen based on their ability to differentiate populations and their power to recapitulate genetic distance.
#' If multiple levels are specified at minIdentity, and meanDistanceLevels, then multiple maps will be produced.
#' The connectedness is aggregated at the desired administrative level: Province (level 1), or District (level 2).
#' To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' Additionally data files are created. Created .tab files can be opened with for example Excel.
#' Maps and data are located in .../output/out/sampleSetName/map-connect/.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet().
#' @param measures default is "ALL", other options are "meanDistance", or "similarity".
#' "meanDistance" produces a map with pairwise comparisons between administrative divisions (e.g. sites) in which the thickness of the connection line corresponds to
#' the proportion of sample pairs with a mean genetic distance of at least the meanDistanceLevels
#' "similarity" produces a map with pairwise comparisons between administrative divisions (e.g. sites) in which the thickness of the connection line corresponds to
#' the proportion of sample pairs with at least the percentage of similarity given in minIdentity.
#' @param minIdentity value between 0 and 1, default is 1
#' @param meanDistanceLevels value between 0 and 1, default is 0.5
#' @param aggregate The administrative level at which we perform pairwise comparisons
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District).
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @export
#'
#' @examples
#' #Map connectedness of haplotypes between districts for a similarity of 1.0, and 0.95, and a mean genetic distance of 0.5, 0.6, and 0.7
#' mapConnections (ctx, sampleSet="Laos",
#' measures="ALL", aggregate="District", minIdentity=c(1.0,0.95),
#' meanDistanceLevels=c(0.5,0.6,0.7),
#' minAggregateCount=10, showNames=TRUE,
#' width=15, height=15)
#
mapConnections <- function (ctx, sampleSet,
                   measures="ALL",
                   minIdentity=1.0,
                   meanDistanceLevels=0.5,
                   aggregate="Province", minAggregateCount=10,
                   showNames=TRUE, markerSize=16,
                   width=15, height=15) {

    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    params <- list(
        analysis.measures=measures,
        aggegation.levels=aggLevels,
        map.aggregateCountMin=minAggregateCount,
        map.markerSize=markerSize,			# Use this for constant marker size
        #map.markerSize=c(4,40),			# Use this for count-proportional marker size
        map.markerNames=showNames,
	map.connect.identity.min=minIdentity,
	map.connect.meanDistance.min=meanDistanceLevels,
        plot.size=list(width=width,height=height)
    )
    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/connect", params=params)
}

#############################################################
#
#' Map barcode group frequencies
#'
#' Creates a map showing the frequency of groups of Plasmodium samples with identical a 'genetic barcode' within administrative divisions.
#' The 'genetic barcode' consists of 101 SNPs determined for each sample. These SNPs are distributed across all nuclear chromosomes.
#' SNP positions were chosen based on their ability to differentiate populations and their power to recapitulate genetic distance.
#' The samples are grouped according to identical barcode within each administrative division, at the desired administrative level: Province, or District.
#' Barcode groups are not compared between administrative divisions. To compare barcode groups between administrative divisions please use mapClusterSharing(). Use command ?mapClusterSharing for more information.
#' To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' Additionally, data files are created. Created .tab files can be opened with for example Excel.
#' Maps and data are located in .../output/out/sampleSetName/map-barcodeFrequency/.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet().
#' @param type Frequency of barcode groups can be visualized either as a "bar" and/or a "pie".
#' @param aggregate The administrative level at which we perform pairwise comparisons, either "Province" and/or "District".
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District).
#' @param markerScale Allows adjustment of the size of markers on the map, default: 0.8.
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @export
#'
#' @examples
#' #Map barcode group frequencies in both bar and pie form, for sampleSet "Laos", on both a Province and District level
#' mapBarcodeFrequencies (ctx, sampleSet="Laos",
#' type=c("bar","pie"),
#' aggregate=c("Province","District"),
#' minAggregateCount=10, showNames=TRUE, markerScale=0.8,
#' width=15, height=15)
#
mapBarcodeFrequencies <- function (ctx, sampleSet,
                   type=c("bar","pie"),
                   aggregate="Province", minAggregateCount=10,
                   showNames=TRUE, markerScale=0.8,
                   width=15, height=15) {

    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    mapParams <- list(
        aggegation.levels=aggLevels,
        map.aggregateCountMin=minAggregateCount,
        map.markerNames=showNames,
        map.cluster.visualizations=type,
        map.cluster.markerScale=markerScale,
        plot.size=list(width=width,height=height)
    )
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/barcodeFrequency", params=mapParams)
}

#############################################################
#
#' Clustering (Updating the analysis context)
#'
#' This function partitions the Plasmodium samples into clusters with a similar genetic background, based on the 'genetic barcode' similarity.
#' The 'genetic barcode' consists of 101 SNPs determined for each sample. These SNPs are distributed across all nuclear chromosomes.
#' SNP positions were chosen based on their ability to differentiate populations and their power to recapitulate genetic distance.
#' The analysis context will be updated by findClusters(), and will be saved in the RStudio environment to be used in future cluster analysis, as seen in examples.
#' Functions providing cluster analyses include: mapClusterSharing(), and mapClusterPrevalence(). Use commands ?mapClusterSharing, and ?mapClusterPrevalence for more information.
#' findClusters() also produces three .tab files that can be opened with for example Excel.
#' Output files are located in .../output/out/(sampleSetName)/cluster/data/(clusterSetname)/ge(minIdentity).
#' The output files include: i) clusterMembers.tab: A list of samples and the cluster they belong to. Samples that are missing from the list are the one that do not belong to a cluster.
#' ii) clusters.tab: A summary of cluster size and members. iii) clusterStats.tab: Frequencies of antimalarial drug resistance predictions and genetic mutation types for each cluster.
#'
#' @param ctx The analysis context, updated with intializeContext().
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet().
#' @param clusterSet The name of the clustering set.
#' @param minIdentity The minimal similarity level set for a pair of samples to be in a cluster. For example, "0.95" corresponds to at least 95 percent genetic barcode similarity.
#' The default is 1.
#' @param clusteringMethod The clustering method. Two methods are available: "allNeighbours" and "louvain". The default is "allNeighbours".
#' The "allNeighbours" method clusters samples together that are above the set "minIdentity" threshold.
#' This method is less informative at low similarity levels, because each sample will be assigned to a single cluster.
#' The "louvain" method is the preferred method. Also known as Louvain Community-based clustering,
#' this method uses an algorithm to identify clusters within a network that are strongly connected to each other, and more weakly connected to other clusters.
#' This method is superior to "allNeighbours", particularly in sample sets that have low genetic similarity.
#' @param minClusterSize To avoid creating very small clusters, one can set a minimum cluster size. The default is 10 samples.
#'
#' @export
#'
#' @examples
#' # Replaces the analysis context to also contain information about clustering
#' # and clusterSet output named "Laos_clusters", containing clusters of at least 5 samples
#' # based on sampleSet "Laos", using the Louvain method at genetic similarity levels of 100 and 95 percent
#' ctx <- findClusters(ctx, sampleSet="Laos", clusterSet = "Laos_clusters",
#' minIdentity = c(1, 0.95),
#' clusteringMethod = "louvain", minClusterSize = 5)
#
findClusters <- function (ctx, sampleSet, clusterSet,
                   minIdentity=1.0,
                   clusteringMethod="allNeighbours",
                   minClusterSize=10) {

    mapParams <- list(
        cluster.clusterSet.name=clusterSet,
        cluster.identity.min=minIdentity,
        cluster.method=clusteringMethod,
        cluster.minSize=minClusterSize
    )
    ctx <- cluster.findClusters (ctx, sampleSetName=sampleSet, params=mapParams)
    ctx
}

#############################################################
#
#' Plot cluster networks based on genetic barcode
#'
#' Creates graphs of cluster networks within a sampleSet of Plasmodium samples.
#' The graphs are based on the cluster analysis context, which needs to be produced with the findClusters function in advance.
#' The findClusters function partitions the Plasmodium samples into clusters with a similar genetic background, based on the 'genetic barcode' similarity.
#' Use command ?findClusters for more information.
#' The 'genetic barcode' consists of 101 SNPs determined for each sample. These SNPs are distributed across all nuclear chromosomes.
#' SNP positions were chosen based on their ability to differentiate populations and their power to recapitulate genetic distance.
#' Additionally, data files are created. Created .tab files can be opened with for example Excel.
#' Maps and data are located in .../output/out/(sampleSetName)/graph/.
#'
#' @param ctx The analysis context, updated with findClusters().
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet().
#' @param clusterSet The name of the clustering set, defined in findClusters().
#' @param graphLayout ?, the default="fr".
#' @param weightPower ?, the default is 2.
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @return Produces network graphs of clusters assigned in findClusters function
#' @export
#'
#' @examples
#' #Plot cluster network graph for sampleSet "Laos", and clusterSet "LAclust"
#' # based on analysis context called ctx updated with findClusters
#' plotClusterGraph (ctx, sampleSet="Laos", clusterSet="LAclust",
#' graphLayout="fr",weightPower=2,
#' width=15,height=15)
#
plotClusterGraph <- function (ctx, sampleSet, clusterSet,
                   graphLayout="fr", weightPower=2,
                   width=15, height=15) {
    mapParams <- list(
        cluster.clusterSet.name=clusterSet,
        graph.layoutAlgorithm=graphLayout,
        graph.weightPower=weightPower,
        plot.size=list(width=width,height=height)
    )
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="graph", params=mapParams)
}

#############################################################
#
#' Map cluster frequencies
#'
#' Creates maps showing the frequency of clusters, allowing for comparison across administrative divisions.
#' Each colour in the bar/pie represents a different cluster.
#' White represents samples that do not belong to a cluster because they did not meet the specified parameter criteria.
#' The maps are based on the cluster analysis context, which needs to be produced with the findClusters function in advance.
#' The findClusters function partitions the Plasmodium samples into clusters with a similar genetic background, based on the 'genetic barcode' similarity.
#' Use command ?findClusters for more information.
#' The 'genetic barcode' consists of 101 SNPs determined for each sample. These SNPs are distributed across all nuclear chromosomes.
#' SNP positions were chosen based on their ability to differentiate populations and their power to recapitulate genetic distance.
#' Additionally, data files are created. Created .tab files can be opened with for example Excel.
#' Maps and data are located in .../output/out/(sampleSetName)/map-clusterSharing/.
#'
#' @param ctx The analysis context, updated with findClusters().
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet().
#' @param clusterSet The name of the clustering set, defined in findClusters().
#' @param type Frequency of clusters can be visualized either as a "bar" and/or a "pie". The default is c("bar", "pie").
#' @param aggregate The administrative level at which we perform pairwise comparisons, either "Province" and/or "District".
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples,
#' below which the marker is not shown.The default is 5.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District).
#' @param markerScale Allows adjustment of the size of markers on the map, default: 0.8.
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @return Maps with cluster frequencies across administrative divisions
#' @export
#'
#' @examples
#' # Map cluster frequencies in both bar and pie form, for sampleSet "Laos", and clusterSet "LAclust" on both a Province and District level,
#' # based on analysis context called ctx updated with findClusters
#' mapClusterSharing (ctx, sampleSet="Laos", clusterSet = "LAclust",
#' type=c("bar", "pie"),
#' aggregate=c("Province","District"),
#' minAggregateCount=5, showNames=TRUE, markerScale=0.8,
#' width=15, height=15)
#'
#
mapClusterSharing <- function (ctx, sampleSet, clusterSet,
                   type=c("bar","pie"),
                   aggregate="Province", minAggregateCount=5,
                   showNames=TRUE, markerScale=0.8,
                   width=15, height=15) {

    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    mapParams <- list(
        cluster.clusterSet.name=clusterSet,
        aggegation.levels=aggLevels,
        map.aggregateCountMin=minAggregateCount,
        map.markerNames=showNames,
        map.cluster.visualizations=type,
        map.cluster.markerScale=markerScale,
        plot.size=list(width=width,height=height)
    )
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/clusterSharing", params=mapParams)
}

#############################################################
#
#' Map cluster networks between sites based on prevalence
#'
#' For each assigned cluster a map is created, showing the prevalence and overlap in prevalence of that cluster across administrative divisions.
#' Each map also contains information on sample count, drug-resistant and mutation frequencies.
#' The maps are based on the cluster analysis context, which needs to be produced with the findClusters function in advance.
#' The findClusters function partitions the Plasmodium samples into clusters with a similar genetic background, based on the 'genetic barcode' similarity.
#' Use command ?findClusters for more information.
#' The 'genetic barcode' consists of 101 SNPs determined for each sample. These SNPs are distributed across all nuclear chromosomes.
#' SNP positions were chosen based on their ability to differentiate populations and their power to recapitulate genetic distance.
#' The prevalence of the assessed cluster is displayed in the nodes.
#' The overlap in prevalence is displayed in the edges (i.e. lines).
#' The definition of prevalence overlap is the proportion of sample pairs that have at least the assigned "minIdentity" of genetic barcode similarity in the findCluster function.
#' Additionally, data files are created. Created .tab files can be opened with for example Excel.
#' Maps and data are located in .../output/out/(sampleSetName)/map-clusterPrevalence/.
#'
#' @param ctx The analysis context, updated with findClusters().
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet().
#' @param clusterSet The name of the clustering set, defined in findClusters().
#' @param aggregate The administrative level at which we perform pairwise comparisons, either "Province" and/or "District". The default is "Province".
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples,
#' below which the marker is not shown. The default is 5.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District).
#' @param markerScale Allows adjustment of the size of markers on the map. The default is 0.8.
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @return For each assigned cluster a map with cluster prevalence and overlap in cluster prevalence across administrative divisions
#' @export
#'
#' @examples
#' # Produce a network and prevalence map for each assigned cluster,
#' # using analysis context called ctx, updated with findClusters()
#' # for sampleSet "Laos", and clusterSet "LAclust" on both a Province and District level,
#' mapClusterPrevalence (ctx, sampleSet="Laos", clusterSet = "LAclust",
#' aggregate=c("Province","District"),
#' minAggregateCount=5, showNames=TRUE, markerScale=0.8,
#' width=15, height=15)
#
mapClusterPrevalence <- function (ctx, sampleSet, clusterSet,
                   aggregate="Province", minAggregateCount=5,
                   showNames=TRUE, markerScale=0.8,
                   width=15, height=15) {

    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    mapParams <- list(
        cluster.clusterSet.name=clusterSet,
        aggegation.levels=aggLevels,
        map.aggregateCountMin=minAggregateCount,
        map.markerNames=showNames,
        map.cluster.visualizations="cluster",
        map.cluster.markerScale=markerScale,
        plot.size=list(width=width,height=height)
    )
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/clusterPrevalence", params=mapParams)
}

#############################################################
#
#' Graphical Attributes for principal component analysis and neighbour-joining tree
#'
#' This function allows for specification of graphical display of samples based on their metadata,
#' that can subsequently be used to plot the population structure with a principal component analysis (plotPrincipalComponents()),
#' or creation of a neighbour-joining tree (plotTree()).
#' Give command ?plotPrincipalComponents, or ?plotTree for more information.
#' Specification based on metadata (i.e. attributes) is done by loading an Excel spreadsheet into R.
#' The spreadsheet and instructions can be found in the grcMalaria user-guide, see link in 'See also' section.
#'
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param name The name of the graphical attribute, this can be any metadata category of interest.
#' Once defined, graphical attributes can be applied by referencing their name.
#' Examples of names are: 'Countries', 'Provinces', 'Kelch13', etc.
#' @param field The corresponding column-name in the original data file of the specified graphical attributes.
#' @param file The file name of the Excel spreadsheet where all graphical attributes are specified.
#' @param sheet The sheet name of the Excel spreadsheet where the to parameter 'name' corresponding graphical attribute is specified.
#'
#' @export
#'
#' @seealso grcMalaria user-guide:
#' \itemize{
#'  \item {\url{https://olivo-miotto-genre.notion.site/grcMalaria-R-package-User-Guide-dfaac7ea30b5430b8522e8b531fd2e4b}}
#' }
#'
#' @examples
#'  #Load Excel spreadsheet into R with specified graphical attributes
#'  gaFile <- "/path/to/my/file/GenRe-GraphicAttributes.xlsx"
#'
#'  #Load graphical attributes for province, Kelch13, and Plasmepsin 2,3
#' ctx <- loadGraphicAttributes (ctx, name="province",      field="AdmDiv1",       file=gaFile, sheet="Provinces")
#' ctx <- loadGraphicAttributes (ctx, name="k13",           field="Kelch",         file=gaFile, sheet="Kelch13")
#' ctx <- loadGraphicAttributes (ctx, name="pm23",          field="Plasmepsin2/3", file=gaFile, sheet="Plasmepsin23")
#' ctx <- loadGraphicAttributes (ctx, name="pm23-noColour", field="Plasmepsin2/3", file=gaFile, sheet="Plasmepsin23-noColour")
#'
loadGraphicAttributes <- function (ctx, name, field, file, sheet) {
    ctx <- graphics.loadAttributes (ctx, name, field, file, sheet)
    ctx
}

#############################################################
#
#' Plot population structure (PCA)
#'
#' This function produces principal component analysis (PCA) plots, in which the population (i.e. samples in the sampleSet),
#' is structured in principal components based on their 'genetic barcode'.
#' The 'genetic barcode' consists of 101 SNPs determined for each sample. These SNPs are distributed across all nuclear chromosomes.
#' SNP positions were chosen based on their ability to differentiate populations and their power to recapitulate genetic distance.
#' One plot is created for each specified attribute in the parameter 'plots', in which the samples are coloured according to user input.
#' Graphical display of attributes needs to be specified before using this function.
#' Give command ?loadGraphicAttributes for more information.
#' Additionally, data files are created. Created .tab files can be opened with for example Excel.
#' Plots and data are located in .../output/out/(sampleSetName)/(PCA-method, see parameter 'type')/.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet().
#' @param type The type of method used for the principal component analysis.
#' "PCoA", performs multidimensional scaling on a distance matrix, in which each sample is a variable.
#' "nipals", "bpca" are methods applied to the set of genetic barcode genotypes (each barcode position is a variable).
#' The default is "PCoA".
#' @param plots The list of attributes for which plots will be created. See example below.
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @return Principal component plots and .tab data files.
#' @export
#'
#' @examples
#' # Perform a principal component analysis for province, kelch13, plasmepsin 2,3, and kelch13 in combination with plasmepsin 2,3,
#' # using the "PCoA" method on sampleSet "Laos"
#' plotPrincipalComponents (ctx, sampleSet="Laos", type="PCoA",
#' plots=list(
#'  list(name="ByProvince",       attributes="province"),
#'  list(name="ByKelch13",        attributes="k13"),
#'  list(name="ByPm23",           attributes="pm23"),
#'  list(name="ByKelch13AndPm23", attributes=c("k13","pm23-noColour"))
#' ), width=15, height=10)
#'
plotPrincipalComponents <- function (ctx, sampleSet,
                                     type="PCoA", plots,
                                     width=15, height=15) {
    plotParams <- list(
        plot.plotList=plots,
        plot.size=list(width=width,height=height)
    )
    task <- paste("pca", type, sep="/")
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks=task, params=plotParams)
}

#############################################################
#
#' Plot Trees
#'
#' This function creates phylogenetic tree plots,
#' in which the evolution of the population (i.e. samples in the sampleSet) is reconstructed based on the 'genetic barcode'.
#' The 'genetic barcode' consists of 101 SNPs determined for each sample. These SNPs are distributed across all nuclear chromosomes.
#' SNP positions were chosen based on their ability to differentiate populations and their power to recapitulate genetic distance.
#' For each specified attribute, a plot is created for the phylogenetic tree, in which the nodes (i.e. samples) are coloured according to that attribute.
#' Graphical display of attributes needs to be specified before using this function.
#' Give command ?loadGraphicAttributes for more information.
#' Additionally, separate .npg legend files are created.
#' Furthermore, data files are created. Created .tab files can be opened with for example Excel.
#' Lastly, a newick file is created that can be used for phylogenetic tree visualization tools.
#' Plots and data are located in .../output/out/(sampleSetName)/njt/.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet().
#' @param type The method used to create a phylogenetic tree.
#' Currently the only option "njt" which creates a neighbour-joining tree.
#' @param plots The list of attributes for which plots will be created. See example below.
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @return Phylogenetic tree plots, .tab data files, and a .newick tree file.
#' @export
#'
#' @examples
#' #Produce a neighbour-joining tree for sampleSet "Laos"
#' #with plots in which the samples (i.e. nodes) are coloured according to i) province, ii) kelch13, iii) plasmepsin 2,3, iv) kelch13 and plasmepsin 2,3
#' plotTree (ctx, sampleSet="Laos", type="njt",
#' plots=list(
#'  list(name="ByProvince",       attributes="province"),
#'  list(name="ByKelch13",        attributes="k13"),
#'  list(name="ByPm23",           attributes="pm23"),
#'  list(name="ByKelch13AndPm23", attributes=c("k13","pm23-noColour"))
#' ), width=15, height=10)
#'
plotTree <- function (ctx, sampleSet,
                      type="njt", plots,
                      width=15, height=15) {
    plotParams <- list(
        plot.plotList=plots,
        plot.size=list(width=width,height=height)
    )
    task <- paste("tree", type, sep="/")
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks=task, params=plotParams)
}
