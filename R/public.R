#############################################################
#
#' grcMalaria: Processing and analyzing SpotMalaria Genetic Report Cards
#'
#'\if{html}{\figure{logo.jpg}{options: align='right' alt='logo' width='120'}}
#' Various functions for producing analysis results such as maps and reports from the Genetic Report Cards datasets
#' delivered by SpotMalaria projects, such as GenRe-Mekong.
#' Typical outputs include drug resistance prevalence maps, population genetic analyses such as PCA (principal component analysis) and NJ (neighbor joining) trees,
#' and analyses of genetic barcodes such as barcode sharing networks, etc.
#' The data frame returned contains a data table whose values and structure may have been manipulated after being read from the file,
#' and may have been updated to the latest version fo the format.
#'
#' @section Author:
#' Olivo Miotto \email{olivo@tropmedres.ac}
#'
#' @seealso Useful links:
#' \itemize{
#'  \item user guide {\url{https://cutt.ly/grcMalaria-R-UserGuide}}
#'  \item source code {\url{https://github.com/malariagen/grcMalaria}}
#'  \item reporting bugs {\url{https://github.com/malariagen/grcMalaria/issues}}
#' }
#'
#' @docType package
#' @name grcMalaria
NULL

#############################################################
#
#' Load a Genetic Report Card (GRC) data file
#'
#' Load a GRC data file, ready for analysis.
#' The data frame returned contains a data table whose values and structure may have been manipulated after being read from the file,
#' and may have been updated to the latest version fo the format.
#'
#' @param file The path to the GRC data Microsoft Excel file (use forward slashes in the path)
#' @param sheet The name of the datasheet or tab name within the Excel file containing the GRC data
#' @param species The species being analyzed ("Pf" (i.e. Plasmodium falciparum) or "Pv" (i.e. Plasmodium vivax))
#' @param version The version number of the GRC data file format. This found is in the documentation when you receive or download the file.
#'
#' @return A list containing a data frame with the data ready to be analyzed, plus some configuration metadata
#' @export
#'
#' @seealso Information on the GRC input file {\url{https://cutt.ly/grcMalaria-R-UserGuide}}
#'
#' @examples \dontrun{
#' # Load data file
#' # Change the path to where your file is located before running the code
#' Data <- loadGrc("D:/myFiles/name_file.xlsx",
#'                sheet = "GRC",
#'                species = "Pf", version = "1.2")
#'}
#
loadGrc <- function (file, sheet="GRC", species="Pf", version="1.2") {
    grcData <- grcData.load (file, sheet, species, version)
    grcData
}

#############################################################
#
#' Combine two Genetic Report Cards (GRC) datasets
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
#' @examples \dontrun{
#' # Load spreadsheet 1
#' Sheet1 <- loadGrc("D:/myFiles/sheet1.xlsx",
#'                  sheet="GRC", species="Pf", version="1.2")
#' # Load spreadsheet 2
#' Sheet2 <- loadGrc("D:/myFiles/sheet2.xlsx",
#'                  sheet="GRC", species="Pf", version="1.2")
#' # Merge datasets
#' mergeData <- mergeGrc(Sheet1, Sheet2)
#'}
#
mergeGrc <- function (srcGrc, newGrc, overwrite=FALSE, extendColumns=FALSE) {
    grcData <- grcData.merge (srcGrc, newGrc, overwrite, extendColumns)
    grcData
}

#############################################################
#
#' Initialize a Genetics Report Card (GRC) dataset
#'
#' Initialize a GRC dataset (obtained from function loadGrc()) so that it is ready for analysis.
#' It performs barcode imputation, computes genetic distances, and initializes data that will subsequently be used in analyses.
#' It creates a folder where output of future analyses will be saved. It will take a while to create a context object which will be used for the subsequent analysis tasks (can take as long as 5 mins to run).
#'
#' @param grcData The data obtained from reading the GRC Excel data.
#' @param dir The folder where the outputs from this and subsequent analyses will be stored.
#' @param minSnpTypability The minimum proportion of non-missing samples for a barcode position to be retained in analysis.
#' @param minSampleTypability The minimum proportion of non-missing positions for a sample to be retained in analysis.
#'
#' @return An analysis context object, which is a list that contains all the data for analysis, which will be passed to subsequent analysis tasks.
#' @export
#'
#' @examples \dontrun{
#' Data <- loadGrc("D:/myFiles/name_file.xlsx", sheet = "GRC",
#'                species = "Pf", version = "1.2")
#' 
#' # Change the path to where you want output file to be
#' ctx <- initializeContext(Data, dir="/path/to/my/folder")
#'}
#
initializeContext <- function (grcData, dir=".", 
                               minSnpTypability=0.8, 
                               minSampleTypability=0.75) {
    options(scipen=10)
    options(stringsAsFactors=FALSE)

    config <- setup.getConfig (grcData, dir, minSnpTypability, minSampleTypability)
    ctx <- context.createRootContext (grcData$data, config)
    ctx
}

#############################################################
#
#' Select a set of samples using metadata
#'
#' Selects a set of samples for analysis, based on their metadata values.
#' A given analysis context can contain multiple sample sets. These will be labelled with different names.
#' The selection criteria are specified as a list containing a sequence of lists wth two elements each:
#' "field" which is the column name to be checked, and "values" which is an array of possible values that can be matched for a sample to be selected.
#' A selected sample must match all the criteria.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSetName The name of the sample set, to be used to identifty it when calling analysis tasks. This will create a folder with the same name e.g. ?Laos?, or "Test170522".
#' @param select The criteria for sample selection. Selected samples will match all specified criteria. The values must be a list. "field" corresponds to a column name in the datafile. "values" correspond to an array of possible values that can be matched for a sample to be selected. Use "," to add additional parameters inside a quotation mark " ".
#'
#' @return The analysis context, augmented with the new sample set
#' @export
#'
#' @examples \dontrun{
#' Data <- loadGrc("D:/myFiles/name_file.xlsx", sheet = "GRC",
#'                species = "Pf", version = "1.2")
#' ctx <- initializeContext(Data, dir="/path/to/my/folder")
#'
#' # Create a SampleSet named "Laos", based on "Timepoint" and "Study"
#' selectSampleSet(ctx, sampleSet="Laos", 
#'                 select=list(
#'                     list(field="TimePoint", values=c("D00H00","-")),
#'                     list(field="Study", values="1208-PF-LA-CMPE-GENRE") 
#'                  ))
#'
#' # Create a SampleSet named "Test170522", containing only isolates from 2018, and 2019 from Vietnam
#' selectSampleSet(ctx, sampleSet="Test170522",
#'                 select=list(
#'                     list(field="Year", values=c("2018","2019")),
#'                     list(field="Country", values="VN")
#'                 ))
#'
#' # Create a SampleSet named "C580Y", containing only isolates with the C580Y genotype
#' selectSampleSet(ctx, sampleSet="C580Y", 
#'                 select=list(
#'                     list(field="Kelch", values="C580Y") 
#'                  ))
#'}
#
selectSampleSet <- function (ctx, sampleSetName, select) {
    sampleSet.selectSampleSet (ctx, sampleSetName, select)
}

#############################################################
#
#' Map of Locations
#'
#' Creates maps showing markers for all Provinces/Districts/Sites where samples were collected.
#' For each location, we place a marker on the map, coloured according to the administrative
#' division it is in.
#' Additionally data files are created. Created .tab files can be opened with for example Excel.
#' Maps and data are located in .../output/out/sampleSetName/map-location/.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param aggregate The administrative level at which we aggregate (Province/District/Site). Separate maps are created for each administrative level.
#' @param markerSize Allows adjustment of the size of markers on the map. If only one value is passed, 
#'                   all markers will be that size; if two values are passed, they will be used as the min 
#'                   and max size of the marker, whose size will reflect the number of samples.
#' @param colourBy Shows the aggregation level to be used to colour the markers (Country or Province)
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... any of plot parameters, including: width, height, units, dpi, legendPosition, legendWidth. 
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default:"NULL"
#'
#' @export
#'
#' @examples \dontrun{
#' # Create 2 maps showing the locations for sampleSet Cambodia on a Province, and a District level
#' mapLocations (ctx, sampleSet="Cambodia", aggregate=c("Province", "District"), 
#'              markerSize=5, colourBy="Country", showNames=FALSE, 
#'              width=18, height=10, units="in", format="pdf")
#'}
#
mapLocations <- function (ctx, sampleSet,
                   aggregate="Province",
                   markerSize=c(4,40), colourBy="Province",
                   showNames=TRUE, nameFontSize=5, 
                   ...) {

    task <- "map/location"
    params <- param.makeParameterList (ctx, task, 
                  aggregate=aggregate, markerSize=markerSize, colourBy=colourBy, showNames=showNames, nameFontSize=nameFontSize, 
                  ...)
    execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
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
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param timePeriods A separately given list of time period objects for partitioning samples into time-interval plots. Each list specifies a time period, with three parameters: name, type and, start. 
#' Parameter type can currently only be year. Parameter start must be a date in the format ?dd-MMM-yyyy?, and parameter name must be provided. The string passed in the name parameter is added to the end of the file name of the produced files.
#' Currently, the public GRC data only contains years, so only the year in the provided ?start? date is used.
#' @param aggregate The administrative level at which we aggregate (Province or District). Separate maps are created for each administrative level.
#' @param minAggregateCount The minimum count of aggregated samples, below which a marker is not shown.
#' @param markerSize Set the size of markers on the map. If only one value is passed, 
#'                   all markers will be that size; if two values are passed, they will be used as the min 
#'                   and max size of the marker, whose size will reflect the number of samples.
#' @param markerFontSize Set the font size of the value labels inside the markers on the map, default=6
#' @param colourBy Shows the aggregation level to be used to colour the markers (Country or Province)
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... any of plot parameters, including: width, height, units, dpi, legendPosition, legendWidth. 
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default:"NULL"
#'
#' @export
#'
#' @examples \dontrun{
#' # Provide lists of time periods of interest, in this case calendar years
#' periods <- list(
#'                list(name="2018", type="year", start="1-Jan-2018"),
#'                list(name="2019", type="year", start="1-Jan-2019"),
#'                list(name="2020", type="year", start="1-Jan-2020")
#'             )
#'
#' # Create maps showing the numbers of samples collected at Province and District aggregation levels.
#' mapSampleCounts (ctx, sampleSet="Laos", timePeriods=periods,
#'                  aggregate=c("Province","District"), minAggregateCount=1,
#'                  markerSize=c(10,40), colourBy="Province", showNames=TRUE)
#'
#' # Provide lists of time periods of interest, in this case the Southern Laos malaria season
#' periods2 <- list(
#'                 list(name="2018", type="year", start="1-Sep-2018"),
#'                 list(name="2019", type="year", start="1-Sep-2019"),
#'                 list(name="2020", type="year", start="1-Sep-2020")
#'             )
#'
#' #Maps showing the numbers of samples collected at Province and District aggregation levels.
#' mapSampleCounts (ctx, sampleSet="Laos", timePeriods=periods2,
#'                  aggregate=c("Province","District"), minAggregateCount=1,
#'                  markerSize=c(10,40), colourBy="Province", showNames=TRUE,
#'                  width=15, height=15, units="in")
#'}
#
mapSampleCounts <- function (ctx, sampleSet, timePeriods=NULL,
                   aggregate="Province", minAggregateCount=1,
                   markerSize=c(4,40), markerFontSize=6, colourBy="Province",
                   showNames=TRUE, nameFontSize=5, 
                   ...) {

    task <- "map/sampleCount"
    params <- param.makeParameterList (ctx, task, 
                  timePeriods=timePeriods, aggregate=aggregate, minAggregateCount=minAggregateCount, markerSize=markerSize, 
                  markerFontSize=markerFontSize, colourBy=colourBy, showNames=showNames, nameFontSize=nameFontSize, 
                  ...)
    execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Map prevalence of Drug resistance (Predicted Phenotype)
#'
#' Creates a map showing the levels of resistance to a particular antimalarial drug. Based on published genetic markers, each sample has a predicted phenotype for different types of antimalarial drugs.
#' If multiple drugs are specified, then different maps will be created for different drugs.
#' The predictions of resistance for any given drugs are aggregated at the desired administrative level: Province (level 1), or District (level2). Separate maps are created for each drug-administrative-combination.
#' For each aggregation unit, we place a marker on the map, coloured according to the level of resistance to the drug, with a label indicating the prevalence.
#' To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' Additionally data files are created. Created .tab files can be opened with for example Excel.
#' Maps and data are located in .../output/out/sampleSetName/map-drug/.
#'
#' @param ctx The analysis context, created by intializeContext()
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet()
#' @param timePeriods A separately given list of time period objects for partitioning samples into time-interval plots. Each list specifies a time period, with three parameters: name, type and, start. 
#' Parameter type can currently only be year. Parameter start must be a date in the format ?dd-MMM-yyyy?, and parameter name must be provided. The string passed in the name parameter is added to the end of the file name of the produced files.
#' Currently, the public GRC data only contains years, so only the year in the provided ?start? date is used
#' @param drugs The antimalarial drugs for which prevalence of phenotypic resistance will be estimated; "ALL" creates maps for all the drugs for which phenotypic resistance predictions are available,
#' which include "Artemisinin", "Chloroquine", "Piperaquine", "DHA-PPQ" (i.e. Dihydroartemisinin/piperaquine), "Mefloquine", "Sulfadoxine", "Pyrimethamine", "S-P" (i.e. Sulphadoxine-Pyrimethamine),
#' "AS-MQ" (i.e. Artesunate-Melfoquine), and "S-P-IPTp" (i.e. intermittent preventive treatment of malaria in pregnancy using Sulfadoxine-Pyrimethamine)
#' To specify a drug put the drug name in between quotation marks e.g. "Artemisinin", or  c(?Artemisinin?, ?Chloroquine?, ?S-P?) to select several specific drugs
#' @param aggregate The administrative level at which we aggregate. Separate maps are created for each administrative level
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown
#' @param markerSize Allows adjustment of the size of markers on the map
#' @param markerFontSize Set the font size of the value labels inside the markers on the map, default=6
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... any of plot parameters, including: width, height, units, dpi, legendPosition, legendWidth. 
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default:"NULL"
#'
#' @export
#'
#' @examples \dontrun{
#' # Provide lists of time periods of interest, in this case calendar years
#' periods <- list(
#'                list(name="2018", type="year", start="1-Jan-2018"),
#'                list(name="2019", type="year", start="1-Jan-2019"),
#'                list(name="2020", type="year", start="1-Jan-2020")
#'             )
#'
#' # Create maps for artemisinin and chloroquine phenotypic resistance prevalence for sampleSet "Laos" 
#' # for both Province and District level for calendar years 2018-2020.
#' mapDrugResistancePrevalence (ctx, sampleSet="Laos", timePeriods=periods,
#'                              drugs=c("Artemisinin", "Chloroquine"), 
#'                              aggregate=c("Province","District"),
#'                              minAggregateCount=10, 
#'                                markerSize=16, showNames=TRUE)
#'}
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
                   markerSize=16, markerFontSize=6,
                   showNames=TRUE, nameFontSize=5, 
                   ...) {

    task <- "map/drug"
    params <- param.makeParameterList (ctx, task, 
                  timePeriods=timePeriods, drugs=drugs, aggregate=aggregate, minAggregateCount=minAggregateCount, markerSize=markerSize, 
                  markerFontSize=markerFontSize, showNames=showNames, nameFontSize=nameFontSize, 
                  ...)
    execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Map the prevalence of genetic mutations of drug resistance
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
#' @param timePeriods A separately given list of time period objects for partitioning samples into time-interval plots. Each list specifies a time period, with three parameters: name, type and, start. 
#' Parameter type can currently only be year. Parameter start must be a date in the format ?dd-MMM-yyyy?, and parameter name must be provided. The string passed in the name parameter is added to the end of the file name of the produced files.
#' Currently, the public GRC data only contains years, so only the year in the provided ?start? date is used.
#' @param mutations The genetic mutation types for which prevalence will be calculated; "ALL" creates maps for all 42 available genetic mutation types.
#' Available mutation types include: "crt_C72S",	"crt_M74I",	"crt_N75E",	"crt_N75D",	"crt_K76T",	"crt_T93S",	"crt_H97Y",	"crt_H97L",	"crt_I218F",	"crt_A220S",
#' "crt_Q271E",	"crt_N326S",	"crt_N326D",	"crt_T333S",	"crt_I356T",	"crt_I356L",	"crt_R371I",	"dhfr_N51I",	"dhfr_C59R",	"dhfr_C59H",	"dhfr_S108N",
#' "dhfr_S108T",	"dhfr_I164L",	"dhps_S436A",	"dhps_S436F",	"dhps_A437G",	"dhps_K540E",	"dhps_K540N",	"dhps_A581G",	"dhps_A613T",	"dhps_A613S",	"mdr1_N86Y",
#' mdr1_Y184F",	"mdr1_S1034I",	"mdr1_F1226Y",	"mdr1_D1246Y",	"arps10_V127M",	"arps10_D128H",	"arps10_D128Y",	"fd_D193Y",	"mdr2_T484I",	"exo_E415G".
#' To specify a mutation put the mutation type in between quotation marks e.g. "crt_N75E", or  c("arps10_D128Y","fd_D193Y","crt_I218F") to select several genetic mutation types.
#' @param aggregate The administrative level at which we aggregate.
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param markerSize Allows adjustment of the size of markers on the map.
#' @param markerFontSize Set the font size of the value labels inside the markers on the map, default=6
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District).
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... any of plot parameters, including: width, height, units, dpi, legendPosition, legendWidth. 
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default:"NULL"
#'
#' @export
#'
#' @seealso Useful links:
#' \itemize{
#'  \item user guide for information on functionality of specific genetic mutation types {\url{https://cutt.ly/grcMalaria-R-UserGuide}}
#' }
#'
#' @examples \dontrun{
#' # Provide lists of time periods of interest, in this case calendar years
#' periods <- list(
#'                list(name="2018", type="year", start="1-Jan-2018"),
#'                list(name="2019", type="year", start="1-Jan-2019"),
#'                list(name="2020", type="year", start="1-Jan-2020")
#'            )
#'
#' # Produce prevalence maps for all available genetic mutation types, for sampleSet "Laos", 
#' # at District level, and for calendar years 2018-2020.
#' mapMutationPrevalence (ctx, sampleSet="Laos", timePeriods=periods,
#'                        mutations="ALL", aggregate="District",
#'                        minAggregateCount=10, showNames=TRUE, markerSize=16)
#'}
#
mapMutationPrevalence <- function (ctx, sampleSet, timePeriods=NULL,
                   mutations="ALL",
                   aggregate="Province", minAggregateCount=10,
                   markerSize=16, markerFontSize=6,
                   showNames=TRUE, nameFontSize=5, 
                   ...) {

    task <- "map/mutation"
    params <- param.makeParameterList (ctx, task, 
                  timePeriods=timePeriods, mutations=mutations, aggregate=aggregate, minAggregateCount=minAggregateCount, 
                  markerSize=markerSize, markerFontSize=markerFontSize, showNames=showNames, nameFontSize=nameFontSize, 
                  ...)
    execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}
                  
#############################################################
#
#' Map the proportion of alleles for a given mutation (including kelch13 and amplifications).
#' 
#' Creates maps showing the proportions of observed allele variants of a given genetic mutation.
#' If multiple genetic mutations are specified, then a map will be created for each mutation.
#' Missing genotypes and heterozygous genotypes are excluded in the calculations.
#' The proportions of observed allele variants are aggregated at the desired administrative level: Province (level 1), or District (level 2). Separate maps are created for each allele-administrative-combination.
#' For each aggregation unit, a pie-chart is created, coloured according to the observed proportions of observed allele variants of the mutation of interest.
#' To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' Additionally data files are created. Created .tab files can be opened with for example Excel.
#' Maps and data are located in .../output/out/sampleSetName/map-alleleProp/.
#'
#' @param ctx The analysis context, created by intializeContext()
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet()
#' @param timePeriods A separately given list of time period objects for partitioning samples into time-interval plots. Each list specifies a time period, with three parameters: name, type and, start. 
#' Parameter type can currently only be year. Parameter start must be a date in the format ?dd-MMM-yyyy?, and parameter name must be provided. The string passed in the name parameter is added to the end of the file name of the produced files.
#' Currently, the public GRC data only contains years, so only the year in the provided ?start? date is used
#' @param mutations The genetic mutation types for which prevalence will be calculated; "ALL" creates maps for all 11 available genetic mutations
#' Available genetic mutations include: "Pfkelch13", "pm23-Amp", "mdr1-Amp", "crt", "dhfr", "dhps", "mdr1", "mdr2", "arps10", "fd", "exo"
#' To specify a mutation put it in between quotation marks e.g. "Pfkelch13", or  c("Pfkelch13","dhfr","fd") to select several genetic mutation types
#' @param aggregate The administrative level at which we aggregate
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown
#' @param markerSize Allows adjustment of the size of markers on the map
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... any of plot parameters, including: width, height, units, dpi, legendPosition, legendWidth. 
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default:"NULL"
#'
#' @export
#' @seealso Useful links:
#' \itemize{
#'  \item user guide for information on functionality of specific genetic mutation types {\url{https://cutt.ly/grcMalaria-R-UserGuide}}
#' }
#'
#' @examples \dontrun{
#' # Produce maps for all available genetic mutations for sample set called Cambodia, 
#' # at province level, with a minimum of 10 samples
#' mapAlleleProportions (ctx, sampleSet="Cambodia", 
#'                       mutations="ALL",
#'                       aggregate="Province", minAggregateCount=10,
#'                       showNames=TRUE, markerSize=c(10,25))
#'}
#
mapAlleleProportions <- function (ctx, sampleSet, timePeriods=NULL,
                   mutations="ALL",
                   aggregate="Province", minAggregateCount=10, 
                   markerSize=16,
                   showNames=TRUE, nameFontSize=5,
                   ...) {

    task <- "map/alleleProp"
    params <- param.makeParameterList (ctx, task, 
                  timePeriods, mutations, aggregate, minAggregateCount, markerSize, showNames, nameFontSize, 
                  timePeriods=timePeriods, mutations=mutations, aggregate=aggregate, minAggregateCount=minAggregateCount, 
                  markerSize=markerSize, showNames=showNames, nameFontSize=nameFontSize, 
                  ...)
    execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
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
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param timePeriods A separately given list of time period objects for partitioning samples into time-interval plots. Each list specifies a time period, with three parameters: name, type and, start. 
#' Parameter type can currently only be year. Parameter start must be a date in the format ?dd-MMM-yyyy?, and parameter name must be provided. The string passed in the name parameter is added to the end of the file name of the produced files.
#' Currently, the public GRC data only contains years, so only the year in the provided ?start? date is used.
#' @param measures This can be "ALL", or any vector containing one or more of c("maxHaploFreq","haploHet", "meanSnpHet","medianDistance").
#' The method "maxHaploFreq" is a measure of loss of diversity. This measure gives the proportion of samples carrying the most common haplotype (defined as samples with identical barcodes).
#' The output ranges from 0-1, with a low value corresponding to a low proportion of samples with the most common haplotype, and a high value corresponding to a
#' high proportion of samples with the most common haplotype.
#' The method "haploHet" is a measure of heterozygosity based on the complete barcode.
#' This measure gives the probability of two randomly selected samples carrying a different barcode, and is useful to detect large changes in a population structure.
#' The output ranges from 0-1, with a low value corresponding to low diversity (low probability of carrying a different barcode), while a high value corresponds to high diversity (high probability of a different barcode).
#' The method "meanSnpHet" is also known as the expected heterozygosity or gene diversity of a locus.
#' For each barcode SNP, expected heterozygosity is calculated using: HE=(n/(n-1))(1-sum(pi^2)),
#' where n = the number of samples and pi = the allele frequency of the ith SNP in the barcode.
#' The final value shows the mean of heterozygosity across all the loci.
#' This measure is able to detect smaller changes in a population structure, and stable as it uses the mean of all SNPs.
#' The output ranges from 0-1, with a low value corresponding to low diversity (low probability of a different allele), and a high value corresponding to high diversity (high probability of a different allele).
#' The method "medianDistance" is a measure that gives the median genetic distance in a assessed group.
#' The is calculated using pairwise genetic distance as the proportion of SNPs differing between two samples. Then taking the median of these in the group of samples of interest
#' @param aggregate The administrative level at which we aggregate
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown
#' @param markerSize Allows adjustment of the size of markers on the map
#' @param markerFontSize Set the font size of the value labels inside the markers on the map, default=6
#' @param markerColours The colour to indicate the level of genetic diversity. Default: "red3", other examples: "forestgreen", "cornflowerblue","darkgoldenrod3"
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... any of plot parameters, including: width, height, units, dpi, legendPosition, legendWidth. 
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default:"NULL"
#'
#' @export
#'
#' @examples \dontrun{
#' #' #Given lists of time periods of interest, in this case calendar years
#' periods <- list(
#'                list(name="2018", type="year", start="1-Jan-2018"),
#'                list(name="2019", type="year", start="1-Jan-2019"),
#'                list(name="2020", type="year", start="1-Jan-2020")
#'            )
#'
#' # Produce genetic diversity maps for all available measures, for sampleSet "Laos", 
#' # at District level, and for calendar years 2018-2020.
#' mapDiversity (ctx, sampleSet="Laos", timePeriods=periods,
#'               measures="ALL", aggregate="District", markerColours="red3",
#'               minAggregateCount=10, showNames=TRUE, markerSize=16,
#'               width=20, height=20)
#'}
#
mapDiversity <- function (ctx, sampleSet, timePeriods=NULL,
                   measures="ALL",
                   aggregate="Province", minAggregateCount=10,
                   markerSize=16, markerFontSize=6, markerColours="red3",
                   showNames=TRUE, nameFontSize=5, 
                   ...) {

    task <- "map/diversity"
    params <- param.makeParameterList (ctx, task, 
                  timePeriods=timePeriods, measures=measures, aggregate=aggregate, minAggregateCount=minAggregateCount, 
                  markerSize=markerSize, markerFontSize=markerFontSize, markerColours=markerColours, showNames=showNames, nameFontSize=nameFontSize, 
                  ...)
    execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
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
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param measures default is "ALL", other options are "meanDistance", or "similarity".
#' "meanDistance" produces a map with pairwise comparisons between administrative divisions (e.g. sites) in which the thickness of the connection line corresponds to
#' the proportion of sample pairs with a mean genetic distance of at least the meanDistanceLevels
#' "similarity" produces a map with pairwise comparisons between administrative divisions (e.g. sites) in which the thickness of the connection line corresponds to
#' the proportion of sample pairs with at least the percentage of similarity given in minIdentity.
#' @param minIdentity value between 0 and 1, default is 1 
#' @param meanDistanceLevels meanDistanceLevels value between 0 and 1, default is 0.5
#' @param aggregate The administrative level at which we perform pairwise comparisons
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param markerSize Allows adjustment of the size of markers on the map.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District).
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... any of plot parameters, including: width, height, units, dpi, legendPosition, legendWidth. 
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default:"NULL"
#'
#' @export
#'
#' @examples \dontrun{
#' # Map connectedness of haplotypes between districts for a similarity of 1.0, and 0.95, 
#' # and a mean genetic distance of 0.5, 0.6, and 0.7
#' mapConnections (ctx, sampleSet="Laos",
#'                 measures="ALL", aggregate="District", minIdentity=c(1.0,0.95),
#'                 meanDistanceLevels=c(0.5,0.6,0.7),
#'                 minAggregateCount=10, showNames=TRUE)
#'}
#
mapConnections <- function (ctx, sampleSet,
                   measures="ALL",
                   minIdentity=1.0,
                   meanDistanceLevels=0.5,
                   aggregate="Province", minAggregateCount=10,
                   markerSize=6,
                   showNames=TRUE, nameFontSize=5,
                   ...) {

    task <- "map/connect"
    params <- param.makeParameterList (ctx, task, 
                  measures=measures, minIdentity=minIdentity, meanDistanceLevels=meanDistanceLevels, 
                  aggregate=aggregate, minAggregateCount=minAggregateCount, markerSize=markerSize, 
                  showNames=showNames, nameFontSize=nameFontSize, 
                  ...)
    execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
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
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param type Frequency of barcode groups can be visualized either as a "bar" and/or a "pie".
#' @param aggregate The administrative level at which we perform pairwise comparisons
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param markerScale Allows adjustment of the size of markers on the map, default: 0.8.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District).
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... any of plot parameters, including: width, height, units, dpi, legendPosition, legendWidth. 
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default:"NULL"
#'
#' @export
#'
#' @examples \dontrun{
#' # Map barcode group frequencies in both bar and pie form, for sampleSet "Laos", 
#' # at both Province and District level
#' mapBarcodeFrequencies (ctx, sampleSet="Laos",
#'                        type=c("bar","pie"),
#'                        aggregate=c("Province","District"),
#'                        minAggregateCount=10, showNames=TRUE, markerScale=0.8)
#'}
#
mapBarcodeFrequencies <- function (ctx, sampleSet,
                   type=c("bar","pie"),
                   aggregate="Province", minAggregateCount=10, 
                   markerScale=0.8,
                   showNames=TRUE, nameFontSize=5, 
                   ...) {

    task <- "map/barcodeFrequency"
    params <- param.makeParameterList (ctx, task, 
                  type=type, aggregate=aggregate, minAggregateCount=minAggregateCount, markerScale=markerScale, 
                  showNames=showNames, nameFontSize=nameFontSize, 
                  ...)
    execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Clustering
#'
#' This function partitions the Plasmodium samples into clusters with a similar genetic background, based on the genetic barcode similarity.
#' The genetic barcode consists of 101 SNPs determined for each sample. These SNPs are distributed across all nuclear chromosomes.
#' SNP positions were chosen based on their ability to differentiate populations and their power to recapitulate genetic distance.
#' Functions providing cluster analyses include: mapClusterSharing(), and mapClusterPrevalence(). Use commands ?mapClusterSharing, and ?mapClusterPrevalence for more information.
#' findClusters() also produces three .tab files that can be opened with for example Excel.
#' Output files are located in .../output/out/(sampleSetName)/cluster/data/(clusterSetname)/ge(minIdentity).
#' The output files include: i) clusterMembers.tab: A list of samples and the cluster they belong to. Samples that are missing from the list are the one that do not belong to a cluster.
#' ii) clusters.tab: A summary of cluster size and members. iii) clusterStats.tab: Frequencies of antimalarial drug resistance predictions and genetic mutation types for each cluster.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param clusterSet The name of the clustering set.
#' @param minIdentity The minimal similarity level set for a pair of samples to be in a cluster. For example, "0.95" corresponds to at least 95 percent genetic barcode similarity.
#' The default is 1.
#' @param impute To use imputed or filtered data. The default is TRUE.
#' @param clusteringMethod The clustering method. Two methods are available: "allNeighbours" and "louvain". The default is "louvain".
#' The "allNeighbours" method clusters samples together that are above the set "minIdentity" threshold.
#' This method is less informative at low similarity levels, because each sample will be assigned to a single cluster.
#' The "louvain" method is the preferred method. Also known as Louvain Community-based clustering,
#' this method uses an algorithm to identify clusters within a network that are strongly connected to each other, and more weakly connected to other clusters.
#' This method is superior to "allNeighbours", particularly in sample sets that have low genetic similarity.
#' @param minClusterSize To avoid creating very small clusters, one can set a minimum cluster size. The default is 10 samples.
#'
#' @export
#'
#' @examples \dontrun{
#' # Partition sampleset "Laos" into clusters of at least 5 samples, using the Louvain method, 
#' # with a genetic similarity threshold of 100 and 95 percent
#' # The clusterset will be called "Laos_clusters"
#' findClusters(ctx, sampleSet="Laos", clusterSet="Laos_clusters",
#'              minIdentity=c(1, 0.95),
#'              clusteringMethod="louvain", minClusterSize=5)
#'}
#
findClusters <- function (ctx, sampleSet, clusterSet,
                   minIdentity=1.0, impute=TRUE,
                   clusteringMethod="louvain", 
                   minClusterSize=10) {
    # Construct a list of arguments so we can use the params.getArgParameter utility function
    args <- list(clusterSet=clusterSet, minIdentity=minIdentity, impute=impute, clusteringMethod=clusteringMethod, minClusterSize=minClusterSize)
    p <- new.env ()
    p$cluster.clusterSet.name <- param.getArgParameter (args, "clusterSet")
    p$cluster.identity.min    <- param.getArgParameter (args, "minIdentity",      type="numeric", multiValue=TRUE, defaultValue=1.0)
    p$cluster.impute          <- param.getArgParameter (args, "impute",           type="logical", defaultValue=TRUE)
    p$cluster.method          <- param.getArgParameter (args, "clusteringMethod", defaultValue="louvain", validValues=c("louvain","allNeighbours"))
    p$cluster.minSize         <- param.getArgParameter (args, "minClusterSize",   type="integer", defaultValue=10)
    cluster.findClusters (ctx, sampleSetName=sampleSet, params=p)
}

#############################################################
#
#' Plot graphs of barcode identity, using cluster
#'
#' Creates graphs of cluster networks within a sampleset of Plasmodium samples.
#' The graphs are based on a cluster analysis, which needs to be produced with the findClusters function in advance.
#' The findClusters function partitions the Plasmodium samples into clusters with a similar genetic background, based on the 'genetic barcode' similarity.
#' Use command ?findClusters for more information.
#' The genetic barcode consists of 101 SNPs determined for each sample. These SNPs are distributed across all nuclear chromosomes.
#' SNP positions were chosen based on their ability to differentiate populations and their power to recapitulate genetic distance.
#' Additionally, data files are created. Created .tab files can be opened with for example Excel.
#' Maps and data are located in .../output/out/(sampleSetName)/graph/.
#'
#' @param ctx The analysis context, created by intializeContext()
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet()
#' @param clusterSet The name of the clustering set, defined in findClusters()
#' @param graphLayout TBD, the default="fr"
#' @param weightPower TBD, the default is 2
#' @param ... any of plot parameters, including: width, height, units, dpi, legendPosition, legendWidth
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default:"NULL"
#'
#' @return Produces network graphs of clusters assigned in findClusters function.
#' @export
#'
#' @examples \dontrun{
#' #Plot cluster network graph for sampleset "Laos", and clusterset "LAclust"
#' plotClusterGraph (ctx, sampleSet="Laos", clusterSet="LAclust",
#'                   graphLayout="fr", weightPower=2)
#'}
#
plotClusterGraph <- function (ctx, sampleSet, clusterSet,
                              graphLayout="fr", weightPower=2,
                              ...) {
                   
    task <- "graph"
    params <- param.makeParameterList (ctx, task,
                  clusterSet=clusterSet, graphLayout=graphLayout, weightPower=weightPower,
                  ...)
    execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Map cluster sharing
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
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param clusterSet The name of the clustering set, defined in findClusters().
#' @param type Frequency of clusters can be visualized either as a "bar" and/or a "pie". The default is c("bar", "pie").
#' @param aggregate The administrative level at which we perform pairwise comparisons, either "Province" and/or "District".
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples,
#' below which the marker is not shown.The default is 5.
#' @param markerScale Allows adjustment of the size of markers on the map, default: 0.8.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District).
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... any of plot parameters, including: width, height, units, dpi, legendPosition, legendWidth
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default:"NULL"
#'
#' @return Maps with cluster frequencies across administrative divisions
#' @export
#'
#' @examples \dontrun{
#' # Map cluster frequencies in both bar and pie form, for sampleset "Laos", 
#' # and clusterset "LAclust" at both Province and District level,
#' mapClusterSharing (ctx, sampleSet="Laos", clusterSet = "LAclust",
#'                    type=c("bar", "pie"),
#'                    aggregate=c("Province","District"),
#'                    minAggregateCount=5, showNames=TRUE, markerScale=0.8)
#'}
#
mapClusterSharing <- function (ctx, sampleSet, clusterSet,
                   type=c("bar","pie"),
                   aggregate="Province", minAggregateCount=5, 
                   markerScale=0.8, 
                   showNames=TRUE, nameFontSize=5, 
                   ...) {
                   
    task <- "map/clusterSharing"
    params <- param.makeParameterList (ctx, task, 
                  clusterSet=clusterSet, type=type, aggregate=aggregate, minAggregateCount=minAggregateCount, markerScale=markerScale, 
                  showNames=showNames, nameFontSize=nameFontSize, 
                  ...)
    execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Map cluster prevalence at different sites, connecting these sites
#'
#' For each assigned cluster a map is created, showing the prevalence and overlap in prevalence of that cluster across administrative divisions.
#' Each map also contains information on sample count, drug-resistance, and mutation frequencies.
#' The maps are based on the cluster analysis, which needs to be produced with the findClusters function in advance.
#' The findClusters function partitions the Plasmodium samples into clusters with a similar genetic background, based on the genetic barcode similarity.
#' Use command ?findClusters for more information.
#' The genetic barcode consists of 101 SNPs determined for each sample. These SNPs are distributed across all nuclear chromosomes.
#' SNP positions were chosen based on their ability to differentiate populations and their power to recapitulate genetic distance.
#' The prevalence of the assessed cluster is displayed in the nodes.
#' The overlap in prevalence is displayed in the edges (i.e. lines).
#' The definition of prevalence overlap is the proportion of sample pairs that have at least the assigned "minIdentity" of genetic barcode similarity in the findCluster function.
#' Additionally, data files are created. Created .tab files can be opened with for example Excel.
#' Maps and data are located in .../output/out/(sampleSetName)/map-clusterPrevalence/.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param clusterSet The name of the clustering set, defined in findClusters().
#' @param aggregate The administrative level at which we perform pairwise comparisons, either "Province" and/or "District".
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples,
#' below which the marker is not shown, default=5
#' @param markerSize Allows adjustment of the size of markers on the map.
#' @param markerFontSize Set the font size of the value labels inside the markers on the map, default=6
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... any of plot parameters, including: width, height, units, dpi, legendPosition, legendWidth
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default:"NULL"
#'
#' @return For each assigned cluster a map with cluster prevalence and overlap in cluster prevalence across administrative divisions
#' @export
#'
#' @examples \dontrun{
#' # Produce a network and prevalence map for each assigned cluster,
#' # for sampleSet "Laos", and clusterSet "LAclust" at both a Province and District level,
#' mapClusterPrevalence (ctx, sampleSet="Laos", clusterSet = "LAclust",
#'                       aggregate=c("Province","District"),
#'                       minAggregateCount=5, showNames=TRUE)
#'}
#
mapClusterPrevalence <- function (ctx, sampleSet, clusterSet,
                   aggregate="Province", minAggregateCount=5, 
                   markerSize=16, markerFontSize=6,
                   showNames=TRUE, nameFontSize=5, 
                   ...) {
                   
    task <- "map/clusterPrevalence"
    params <- param.makeParameterList (ctx, task, 
                  clusterSet=clusterSet, aggregate=aggregate, minAggregateCount=minAggregateCount, markerSize=markerSize, 
                  markerFontSize=markerFontSize, showNames=showNames, nameFontSize=nameFontSize, 
                  ...)
    execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Set to a user-defined colour palette
#'
#' Users can define their own colour palette instead of using the default colours. 
#' This function needs to be executed before sample selection, and after running initializeContext(), otherwise the defined colour palette will not be applied.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param palette The desired palette, specified as separate HTML colour codes
#'
#' @export
#'
#' @examples \dontrun{
#' # To view default 25-colour palette
#' ctx$config$defaultPalette 
#' 
#' # Set colours
#' setColourPalette(ctx, palette=c("red","lightblue","#93c4d2"))
#'
#' # then select sample set
#' selectSampleSet(ctx, sampleSet="SamplesetName",
#'                 select=list(
#'                            list(field="Species", values=c("Pf"))
#'                 ))
#'
#' # then proceed with the analysis
#'}
#
setColourPalette <- function (ctx, palette) {
    graphics.setColourPalette (ctx, palette)
}

#############################################################
#
#' Reset to the default colour palette, removing user-defined palettes
#'
#' Return to the default 25-colour palette of grcMalaria.
#'
#' @param ctx The analysis context, created by intializeContext().
#'
#' @return Return to the default colour palette of grcMalaria.
#' @export
#'
#' @examples \dontrun{
#' # return to the default 25-colour palette:
#' resetColourPalette(ctx)
#'
#' # To view default 25-colour palette
#' ctx$config$defaultPalette
#'}
#
resetColourPalette <- function (ctx) {
    graphics.resetColourPalette (ctx)
}

#############################################################
#
#' Graphical Attributes for Principal Component Analysis (PCA) and Neighbor-Joining (NJ)-tree
#'
#' This function allows for specification of the graphical display of samples based on their metadata,
#' that can subsequently be used to plot the population structure with a principal component analysis (plotPrincipalComponents()),
#' or creation of a neighbour-joining tree (plotTree()).
#' Give command ?plotPrincipalComponents, or ?plotTree for more information.
#' Specification based on metadata (i.e. attributes) is done by loading an Excel spreadsheet into R.
#' The spreadsheet and instructions can be found in the grcMalaria user-guide, see link in 'See also' section.
#'
#' @param ctx The analysis context, created by intializeContext()
#' @param name The name of the graphical attribute, this can be any metadata category of interest
#' Once defined, graphical attributes can be applied by referencing their name
#' Examples of names are: 'Countries', 'Provinces', 'Kelch13', etc
#' @param field The corresponding column-name in the original data file of the specified graphical attributes
#' @param file The file name of the Excel spreadsheet where all graphical attributes are specified
#' @param sheet The sheet name of the Excel spreadsheet where the to parameter 'name' corresponding graphical attribute is specified
#'
#' @return The analysis context, augmented with the new sample set
#' @export
#'
#' @seealso grcMalaria user-guide:
#' \itemize{
#'  \item {\url{https://cutt.ly/grcMalaria-R-UserGuide}}
#' }
#'
#' @examples \dontrun{
#' #Load Excel spreadsheet into R with specified graphical attributes
#' gaFile <- "/path/to/my/file/GenRe-GraphicAttributes.xlsx"
#'
#' #Load graphical attributes for province, kelch13, and plasmepsin2/3 amplification status
#' ctx <- loadGraphicAttributes (ctx, name="province", field="AdmDiv1",       
#'                               file=gaFile, sheet="Provinces")
#' ctx <- loadGraphicAttributes (ctx, name="k13", field="Kelch",
#'                               file=gaFile, sheet="Kelch13")
#' ctx <- loadGraphicAttributes (ctx, name="pm23", field="Plasmepsin2/3",
#'                               file=gaFile, sheet="Plasmepsin23")
#' ctx <- loadGraphicAttributes (ctx, name="pm23-noColour", field="Plasmepsin2/3",
#'                               file=gaFile, sheet="Plasmepsin23-noColour")
#'}
#
loadGraphicAttributes <- function (ctx, name, field, file, sheet) {
    graphics.loadAttributes (ctx, name, field, file, sheet)
}

#############################################################
#
#' Plot Population structure (PCA)
#'
#' This function produces principal component analysis (PCA) plots, in which the population (i.e. samples in the sampleSet),
#' is structured in principal components based on their genetic barcode.
#' The genetic barcode consists of 101 SNPs determined for each sample. These SNPs are distributed across all nuclear chromosomes.
#' SNP positions were chosen based on their ability to differentiate populations and their power to recapitulate genetic distance.
#' One plot is created for each specified attribute in the parameter 'plots', in which the samples are coloured according to user input.
#' Graphical display of attributes needs to be specified before using this function.
#' Give command ?loadGraphicAttributes for more information.
#' Additionally, data files are created. Created .tab files can be opened with for example Excel.
#' Plots and data are located in .../output/out/(sampleSetName)/(PCA-method, see parameter 'type')/.
#'
#' @param ctx The analysis context, created by intializeContext()
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet()
#' @param type The type of method used for the principal component analysis
#' "PCoA", performs multidimensional scaling on a distance matrix, in which each sample is a variable
#' "nipals", "bpca" are methods applied to the set of genetic barcode genotypes (each barcode position is a variable).
#' The default is "PCoA"
#' @param plots The list of attributes for which plots will be created. See example below
#' @param ... any of plot parameters, including: width, height, units, dpi, legendPosition, legendWidth
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default:"NULL"
#'
#' @return Principal component plots and .tab data files.
#' @export
#'
#' @examples \dontrun{
#' # Perform a principal component analysis for province, kelch13, plasmepsin2/3 amplification status, 
#' # and kelch13 in combination with plasmepsin2/3 amplification status,
#' # using the "PCoA" method on sampleSet "Laos"
#' plotPrincipalComponents (ctx, sampleSet="Laos", type="PCoA",
#'     plots=list(
#'               list(name="ByProvince",       attributes="province"),
#'               list(name="ByKelch13",        attributes="k13"),
#'               list(name="ByPm23",           attributes="pm23"),
#'               list(name="ByKelch13AndPm23", attributes=c("k13","pm23-noColour"))
#'            ), 
#'     width=15, height=10)
#'}
#
plotPrincipalComponents <- function (ctx, sampleSet, 
                                     type="PCoA", plots,
                                     ...) {
                                     
    # Construct a list of arguments so we can use the params.getArgParameter utility function
    args <- list(type=type, plots=plots)
    pcaType  <- param.getArgParameter (args, "type", defaultValue="PCoA", validValues=c("PCoA", "nipals", "bpca"))
    task <- paste("pca", pcaType, sep="/")
    params <- param.makeParameterList (ctx, task,
                  plots=plots, 
                  ...)
    execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
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
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param type The method used to create the tree.
#' Currently the only option "njt" which creates a neighbour-joining tree.
#' @param plots The list of attributes for which plots will be created. See example below.
#' @param ... any of plot parameters, including: width, height, units, dpi, legendPosition, legendWidth
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default:"NULL"
#'
#' @return Phylogenetic tree plots, .tab data files, and a .newick tree file.
#' @export
#'
#' @examples \dontrun{
#' # Produce a neighbour-joining tree for sampleSet "Laos"
#' # with plots in which the samples (i.e. nodes) are coloured according to i) province, ii) kelch13, 
#' # iii) plasmepsin2/3 amplification status, iv) kelch13 and plasmepsin2/3 amplification status
#' plotTree (ctx, sampleSet="Laos", type="njt",
#'           plots=list(
#'                     list(name="ByProvince",       attributes="province"),
#'                     list(name="ByKelch13",        attributes="k13"),
#'                     list(name="ByPm23",           attributes="pm23"),
#'                     list(name="ByKelch13AndPm23", attributes=c("k13","pm23-noColour"))
#'           ), 
#'           width=15, height=10)
#'}
#
plotTree <- function (ctx, sampleSet, 
                      type="njt", plots,
                      ...) {

    # Construct a list of arguments so we can use the params.getArgParameter utility function
    args <- list(type=type, plots=plots)
    treeType  <- param.getArgParameter (args, "type", defaultValue="njt", validValues=c("njt"))
    task <- paste("tree", treeType, sep="/")
    params <- param.makeParameterList (ctx, task,
                  plots=plots, 
                  ...)
    execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}
#
#
#
#############################################################
