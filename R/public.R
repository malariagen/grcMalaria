
#############################################################
#
#' Load a Genetic Report Cards data file
#'
#' Load a GRC data file, ready for analysis.
#' The data frame returned contains a data table whose values and structure may have been manipulated after being read from the file,
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
#' #TBD
#
loadGrc <- function (file, sheet="GenRe-Mekong", species="Pf", version="1.0") {
    grcData <- grcData.load (file, sheet, species, version)
    grcData
}

#############################################################
#
#' Initialize a Genetics Report Cards dataset
#'
#' Initialize a GRC dataset (obtained from function loadGrc()) so that it is ready for analysis.
#' It performs barcode imputation, computes genetic distances, and initializes data that will subsequently be used in analyses.
#'
#' @param grcData The data obtained from reading the GRC Excel data.
#' @param dir The folder where the outputs from this and subsequent analyses will be stored.
#' @param minSnpTypability The minimum proportion of non-missing samples for a barcode position to be retained in analysis.
#' @param minSampleTypability The minimum proportion of non-missing positions for a sample to be retained in analysis.
#'
#' @return An analysis context object, which is a list that contains all the data for analysis, which will be passed to subsequent analysis tasks.
#' @export
#'
#' @examples
#' #TBD
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
#' Selects a set of samples for analysis, based on their metadata values.
#' A given analysis context can contain multiple sampe sets, with different names.
#' The selection criteria are specified as a list containing a sequence of lists wth two elements each:
#' "field" which is the column name to be checked, and "values" which is an array of possible values that can be matched for a sample to be selected.
#' A selected sample must match all the criteria.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set, to be used to identifty it when calling analysis tasks.
#' @param select The criteria for selecting the samples. The value must be a list.
#'
#' @return The analysis context, augmented with the new sample set
#' @export
#'
#' @examples
#' #TBD
#
selectSampleSet <- function (ctx, sampleSetName, select) {
    ctx <- analysis.selectSampleSet (ctx, sampleSetName, select)
    ctx
}

#############################################################
#
#' Map prevalence of Drug resistance
#'
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
#' #TBD
#
mapDrugResistancePrevalence <- function (ctx, sampleSet,
                   drugs="ALL",
                   aggregate="Province", minAggregateCount=10,
                   showNames=TRUE, markerSize=16,
                   width=15, height=15) {

    params <- list(
        map.aggregateCountMin=minAggregateCount,
        map.markerSize=markerSize,			# Use this for constant marker size
        #map.markerSize=c(4,40),			# Use this for count-proportional marker size
        map.markerNames=showNames,
        map.size=c(width=width,height=height)
    )
    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/drug", plotList=NULL,
                                 aggregation=aggLevels, measures=drugs, params=params)
}
#
#############################################################
#
#' Map the prevalence of genetic mutations.
#'
#' @param ctx TBD
#' @param sampleSet TBD
#' @param mutations TBD
#' @param aggregate TBD
#' @param minAggregateCount TBD
#' @param showNames TBD
#' @param markerSize TBD
#' @param width TBD
#' @param height TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
mapMutationPrevalence <- function (ctx, sampleSet,
                   mutations="ALL",
                   aggregate="Province", minAggregateCount=10,
                   showNames=TRUE, markerSize=16,
                   width=15, height=15) {

    params <- list(
        map.aggregateCountMin=minAggregateCount,
        map.markerSize=markerSize,			# Use this for constant marker size
        #map.markerSize=c(4,40),			# Use this for count-proportional marker size
        map.markerNames=showNames,
        map.size=c(width=width,height=height)
    )
    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/mutation", plotList=NULL,
                                 aggregation=aggLevels, measures=mutations, params=params)
}
                   
#############################################################
#
#' Map genetic Diversity
#'
#' @param ctx TBD
#' @param sampleSet TBD
#' @param measures TBD
#' @param aggregate TBD
#' @param minAggregateCount TBD
#' @param showNames TBD
#' @param markerColours TBD
#' @param markerSize TBD
#' @param width TBD
#' @param height TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
mapDiversity <- function (ctx, sampleSet,
                   measures="ALL",
                   aggregate="Province", minAggregateCount=10,
                   showNames=TRUE, markerSize=16,
                   markerColours="red3",
                   width=15, height=15) {
    
    params <- list(
        map.aggregateCountMin=minAggregateCount,
        map.markerSize=markerSize,			# Use this for constant marker size
        #map.markerSize=c(4,40),			# Use this for count-proportional marker size
        map.markerNames=showNames,
        map.diversity.markerColours=markerColours,
        map.size=c(width=width,height=height)
    )
    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/diversity", plotList=NULL,
                                 aggregation=aggLevels, measures=measures, params=params)
}

#############################################################
#
#' Map Connectedness between sites
#'
#' @param ctx TBD
#' @param sampleSet TBD
#' @param measures TBD
#' @param similarityLevels TBD TBD
#' @param aggregate TBD
#' @param minAggregateCount TBD
#' @param showNames TBD
#' @param markerColours TBD
#' @param markerSize TBD
#' @param width TBD
#' @param height TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
mapConnections <- function (ctx, sampleSet,
                   measures="ALL",
                   similarityLevels=1.0,
                   meanDistanceLevels=0.5,
                   aggregate="Province", minAggregateCount=10,
                   showNames=TRUE,
                   markerSize=16,
                   width=15, height=15) {
                   
    params <- list(
        map.aggregateCountMin=minAggregateCount,
        map.markerSize=markerSize,			# Use this for constant marker size
        #map.markerSize=c(4,40),			# Use this for count-proportional marker size
        map.markerNames=showNames,
	map.connect.similarity.min=similarityLevels,
	map.connect.meanDistance.min=meanDistanceLevels,
        map.size=c(width=width,height=height)
    )
    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/connect", plotList=NULL,
                                 aggregation=aggLevels, measures=measures, params=params)
}

#############################################################
#
#' Map barcode group frequencies
#'
#' @param ctx TBD
#' @param sampleSet TBD
#' @param type TBD
#' @param similarityLevels TBD
#' @param aggregate TBD
#' @param minGroupSize TBD
#' @param showNames TBD
#' @param markerScale TBD
#' @param width TBD
#' @param height TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
mapGroupFrequencies <- function (ctx, sampleSet,
                   type=c("bar","pie"),
                   similarityLevels=1.0,
                   aggregate="Province",
                   minGroupSize=10,
                   showNames=TRUE,
                   markerScale=0.8,
                   width=15,
                   height=15) {}

# Q: Check this:     map.haplo.markerSampleCount="mean"  #"none", or a count

#############################################################
#
#' Plot a tree (e.g. neighbour joining tree, NJT)
#'
#' @param ctx TBD
#' @param sampleSet TBD TBD
#' @param type TBD
#' @param attributes TBD
#' @param width TBD
#' @param height TBD
#' @param colouredBranches TBD
#' @param showLeaves TBD
#' @param nonLeafColour TBD
#' @param branchThickness TBD
#' @param popLabels TBD
#' @param sampleLabels TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
plotTree <- function (ctx, sampleSet, method="njt",
                   attributes=NULL,
                   width=15,
                   height=15,
                   colouredBranches=TRUE,
                   showLeaves=TRUE,
                   nonLeafColour="darkgrey",
                   branchThickness=3,
                   popLabels=FALSE,
                   sampleLabels=FALSE) {}

# Q: DEFINE A DEFAULT SET OF ATTRIBUTES

#############################################################
#
#' Plot Population structure (PCA)
#'
#' @param ctx TBD
#' @param sampleSet TBD TBD
#' @param method TBD
#' @param attributes TBD
#' @param width TBD
#' @param height TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
plotPCA <- function (ctx, sampleSet, method="PCoA",
                   attributes=NULL,
                   width=15,
                   height=15) {}

# Q: DEFINE A DEFAULT SET OF ATTRIBUTES

#############################################################
#
#' Plot Population structure (barcode networks)
#'
#' @param ctx TBD
#' @param sampleSet TBD
#' @param attributes TBD
#' @param aggregateCountMin TBD
#' @param width TBD
#' @param height TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
plotBarcodeNetwork <- function (ctx, sampleSet,
                   attributes=NULL,
                   aggregateCountMin=10,
                   width=15,
                   height=15) {}

#############################################################
#
#' Plot barcode similarity graphs
#'
#' @param ctx TBD
#' @param sampleSet TBD
#' @param attributes TBD
#' @param similarityLevels TBD
#' @param minGroupSize TBD
#' @param minConnectionSimilarity TBD
#' @param graphLayout TBD
#' @param weighting TBD
#' @param width TBD
#' @param height TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
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
#
#' Map barcode sharing
#'
#' @param ctx TBD
#' @param sampleSet TBD
#' @param type TBD
#' @param similarityLevels TBD
#' @param aggregate TBD
#' @param minAggregateCount TBD
#' @param showNames TBD TBD
#' @param markerScale TBD
#' @param minGroupSize TBD
#' @param width TBD
#' @param height TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
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
#
#' Map Barcode Group prevalence
#'
#' @param ctx TBD
#' @param sampleSet TBD
#' @param similarityLevels TBD
#' @param aggregate TBD
#' @param minAggregateCount TBD
#' @param showNames TBD
#' @param minGroupSize TBD
#' @param markerScale TBD
#' @param width TBD
#' @param height TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
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
