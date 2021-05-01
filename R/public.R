
#############################################################
# Data file loading
#
#' Load a Genetic Report Cards data file
#' Load a GRC data file, ready for analysis.
#' The data frame returned constaind a data table whose values and structure may have been manipulated after being read from the file,
#' and may have been updated to the latest version fo the format.
#'
#' @param file The path to the GRC data Microsoft Excel file (use forward slashes in the path)
#' @param sheet The name of the datasheet within the Excel file containing the GRC data
#' @param species The species being analyzed ("Pf" or "Pv")
#' @param version The version number of the GRC data file format. This is in the documentation when you receive or download the file.
#'
#' @return A list containing a data frame with the data ready to be analyzed, plus some metadata
#' @export
#'
#' @examples
#' pfData <- loadGrc ("D:/MyData/GRC/20210705-GenRe.xlsx", sheet="GenRe-Mekong", species="Pf", version="1.0")
#
loadGrc <- function (file, sheet="GenRe-Mekong", species="Pf", version="1.0") {
    grcData <- grcData.load (file, sheet, species, version)
    grcData
}

#############################################################
# Initialization
#
#' Initialize a Genetics Report Cards dataset
#' Initialize a GRC dataset (obtained from function loadGrc()) so that it is ready for analysis.
#' It performs barcode imputation, computes genetic distances, and initializes data that will subsequently be used in analyses.
#'
#' @param grcData The data obtained from reading the GRC Excel data.
#' @param dir The folder where the outputs from this and subsequent analyses will be stored.
#'
#' @return An analysis context object, which is a list that contains all the data for analysis, which will be passed to subsequent analysis tasks.
#' @export
#'
#' @examples
#' pfData <- loadGrc ("D:/MyData/GRC/20210705-GenRe.xlsx", sheet="GenRe-Mekong", species="Pf", version="1.0")
#' ctx <- initializeContext (pfData, dir="D:/MyAnalysis/GRC")
#
initializeContext <- function (grcData, dir=".") {
    options(scipen=10)
    options(stringsAsFactors=FALSE)
    
    config <- setup.getConfig (grcData, dir)
    ctx <- analysis.createContext (grcData$data, config)
    ctx
}

#############################################################
# Sample Selection
#
#' Select a set of samples using metadata
#' Selects a set of samples for analysis, based on their metadata values.
#' A given analysis context can contain multiple sampe sets, with different names.
#' The selection criteria are specified as a list containing a sequence of lists wth two elements each:
#' "field" which is the column name to be checked, and "values" which is an array of possible values that can be matched for a sample to be selected.
#' A selected sample must match all the criteria.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param name The name of the sample set, to be used to identifty it when calling analysis tasks.
#' @param select The criteria for selecting the samples. The value must be a list.
#'
#' @return The analysis context, augmented with the new sample set
#' @export
#'
#' @examples
#' # Select all samples from Vietnam and Laos which were collected upon hospital admission (day 0, hour 0)
#' ctx <- selectSampleSet (ctx, sampleSet="Day0-VN-LA",
#'              select=list(
#'                list(field="TimePoint", values="D00H00"),
#'                list(field="Country", values=c("VN","LA"))
#'            )
#
selectSampleSet <- function (ctx, sampleSet, select) {
    ctx <- analysis.selectDataset (ctx, name, select)
    ctx
}

#############################################################
# Drug resistance prevalence maps
#
#' Map prevalence of Drug resistance
#' Creates a map showing the levels of resistance to a particular drug for different administrative divisions.
#' If multiple drugs are specified, then different maps will be created for different drugs.
#' The predictions of resistance for any given drugs are aggregated at the desired administrative level: Province (level 1), or District (level2) and separate maps are created.
#' For each aggregation unit, we place a marker on the map, coloured according to the level of resistance to the drug, with a label indicating the prevalence.
#' To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used; must have been previusly created by selectSampleSet().
#' @param drugs An array of drugs for which prevalence of resistance will be estimated; "ALL" creats maps for all the drugs for which we have resistance predictions.
#' @param aggregate The administrative level at which we aggregate
#' @param minAggregateCount The minimum count of aggregated samples, below which a marker is not shown.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param markerSize Allows adjustment of the size of markers on the map.
#' @param width The width (in inches) of the map image.
#' @param height The heigt (in inches) of the map image.
#'
#' @export
#'
#' @examples
#' # Get a context for a P. falciparum GRC dataset
#' pfData <- loadGrc ("D:/MyData/GRC/20210705-GenRe.xlsx", sheet="GenRe-Mekong", species="Pf", version="1.0")
#' ctx <- initializeContext (pfData, dir="D:/MyAnalysis/GRC")
#
#' # Select all samples from Vietnam and Laos which were collected upon hospital admission (day 0, hour 0)
#' ctx <- selectSampleSet (ctx, sampleSet="Day0-VN-LA",
#'              select=list(
#'                list(field="TimePoint", values="D00H00"),
#'                list(field="Country", values=c("VN","LA"))
#'            )
#'
#' # Map prevalence of artemisinin and piperaquine resistance at provincial and district level in Veitnam and Laos
#' mapDrugResistancePrevalence (ctx, sampleSet="Day0-VN-LA",
#'                    mutation=c("Artemisinin", "Piperaquine"), aggregate=c(1,2),
#'                    minAggregateCount=10)
#'
#
mapDrugResistancePrevalence <- function (ctx, sampleSet,
                   drugs="ALL",
                   aggregate=1, minAggregateCount=10,
                   showNames=TRUE, markerSize=16,
                   width=15, height=15) {

    mapParams <- list(
        map.aggregateCountMin=minAggregateCount,
        map.markerSize=markerSize,			# Use this for constant marker size
        #map.markerSize=c(4,40),			# Use this for count-proportional marker size
        map.markerNames=showNames,
        map.size=c(width=width,height=height)
    )
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/drug", 
                                 aggregation=aggregate, measures=drugs, params=params)
)
#
#############################################################
# Mutation prevalence maps
#
#' Title
#'
#' @param ctx
#' @param sampleSet
#' @param drugs
#' @param aggregate
#' @param minAggregateCount
#' @param showNames
#' @param markerSize
#' @param width
#' @param height
#'
#' @return
#' @export
#'
#' @examples
mapMutationPrevalence <- function (ctx, sampleSet,
                   mutation="ALL",
                   aggregate=1,
                   minAggregateCount=10,
                   showNames=TRUE,
                   markerSize=16,
                   width=15,
                   height=15) {}

#############################################################
# Diversity maps
#
#' Title
#'
#' @param ctx
#' @param sampleSet
#' @param measures
#' @param aggregate
#' @param minAggregateCount
#' @param showNames
#' @param markerColours
#' @param markerSize
#' @param width
#' @param height
#'
#' @return
#' @export
#'
#' @examples
mapDiversity <- function (ctx, sampleSet,
                   measures="ALL",
                   aggregate=1,
                   minAggregateCount=10,
                   showNames=TRUE,
                   markerColours=c("white","red3"),
                   markerSize=16,
                   width=15,
                   height=15) {}

#############################################################
# Connectedness maps
#
#' Title
#'
#' @param ctx
#' @param sampleSet
#' @param measures
#' @param similarityLevels
#' @param aggregate
#' @param minAggregateCount
#' @param showNames
#' @param markerColours
#' @param markerSize
#' @param width
#' @param height
#'
#' @return
#' @export
#'
#' @examples
mapConnections <- function (ctx, sampleSet,
                   measures="ALL",
                   similarityLevels=1.0,
                   aggregate=1,
                   minAggregateCount=10,
                   showNames=TRUE,
                   markerColours=c("white","red3"),
                   markerSize=16,
                   width=15,
                   height=15) {}

# Q: Check this:     map.connect.meanDistance.min=c(0.55,0.6)
#                    measures=allConnectednessMeasures,

#############################################################
# Haplotype frequency maps
#
#' Title
#'
#' @param ctx
#' @param sampleSet
#' @param type
#' @param similarityLevels
#' @param aggregate
#' @param minGroupSize
#' @param showNames
#' @param markerScale
#' @param width
#' @param height
#'
#' @return
#' @export
#'
#' @examples
mapGroupFrequencies <- function (ctx, sampleSet,
                   type=c("bar","pie"),
                   similarityLevels=1.0,
                   aggregate=1,
                   minGroupSize=10,
                   showNames=TRUE,
                   markerScale=0.8,
                   width=15,
                   height=15) {}

# Q: Check this:     map.haplo.markerSampleCount="mean"  #"none", or a count

#############################################################
# Population structure analysis - PCA and NJT
#
#' Title
#'
#' @param ctx
#' @param sampleSet
#' @param type
#' @param attributes
#' @param width
#' @param height
#'
#' @return
#' @export
#'
#' @examples
plotPopulationStructure <- function (ctx, sampleSet,
                   type=c("njt","PCoA"),
                   attributes=NULL,
                   width=15,
                   height=15) {}

# Q: DEFINE A DEFAULT SET OF ATTRIBUTES

#############################################################
# Population structure analysis - haplotype networks
#
#' Title
#'
#' @param ctx
#' @param sampleSet
#' @param attributes
#' @param aggregateCountMin
#' @param width
#' @param height
#'
#' @return
#' @export
#'
#' @examples
plotBarcodeNetwork <- function (ctx, sampleSet,
                   attributes=NULL,
                   aggregateCountMin=10,
                   width=15,
                   height=15) {}

#############################################################
# Population structure analysis - haplotype similarity graphs
#
#' Title
#'
#' @param ctx
#' @param sampleSet
#' @param attributes
#' @param similarityLevels
#' @param minGroupSize
#' @param minConnectionSimilarity
#' @param graphLayout
#' @param weighting
#' @param width
#' @param height
#'
#' @return
#' @export
#'
#' @examples
mapGroupGraph <- function (ctx, sampleSet,
                   attributes=NULL,
                   similarityLevels=1.0,
                   minGroupSize=5,
                   minConnectionSimilarity=0.7,
                   graphLayout="fr",
                   weighting="identitySquared",
                   width=15,
                   height=15) {}

# Q: No min aggregation count?

#############################################################
# Haplotype sharing maps
#
#' Title
#'
#' @param ctx
#' @param sampleSet
#' @param type
#' @param similarityLevels
#' @param aggregate
#' @param minAggregateCount
#' @param showNames
#' @param markerScale
#' @param minGroupSize
#' @param width
#' @param height
#'
#' @return
#' @export
#'
#' @examples
mapGroupSharing <- function (ctx, sampleSet,
                   type=c("bar","pie"),
                   similarityLevels=1.0,
                   aggregate=1,
                   minAggregateCount=5,
                   showNames=TRUE,
                   markerScale=0.8,
                   minGroupSize=10,
                   width=15,
                   height=15) {}

# Q: Check this:     map.haplo.markerSampleCount="mean"  #"none", or a count

#############################################################
# Haplotype sharing maps - subgroup mapping
#
#' Title
#'
#' @param ctx
#' @param sampleSet
#' @param similarityLevels
#' @param aggregate
#' @param minAggregateCount
#' @param showNames
#' @param minGroupSize
#' @param markerScale
#' @param width
#' @param height
#'
#' @return
#' @export
#'
#' @examples
mapGroupPrevalence <- function (ctx, sampleSet,
                   similarityLevels=1.0,
                   aggregate=1,
                   minAggregateCount=5,
                   showNames=TRUE,
                   minGroupSize=10,
                   markerScale=0.8,
                   width=15,
                   height=15) {}

# Q: Check this:     map.haplo.markerSampleCount="mean"  #"none", or a count



