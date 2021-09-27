
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
#' @return A list containing a data frame with the data ready to be analyzed, plus some configuration metadata
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
#' #TBD
#
mergeGrc <- function (srcData, newData, overwrite=FALSE, extendColumns=FALSE) {
    grcData <- grcData.merge (srcData, newData, overwrite, extendColumns)
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
#' Map of Sample Counts
#'
#' Creates a map showing the numbers of samples collected at different aggregation levels.
#' For each aggregation unit, we place a marker on the map, coloured according to the administrative
#' division it is in
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used; must have been previusly created by selectSampleSet().
#' @param aggregate The administrative level at which we aggregate (Province or District)
#' @param minAggregateCount The minimum count of aggregated samples, below which a marker is not shown.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param colourBy Shows the aggregation level to be used to colour the markers (Country or Province)
#' @param markerSize Allows adjustment of the size of markers on the map. If only one value is passed, 
#'                   all markers will be that size; if two values are passed, they will be used as the min 
#'                   and max sizeof the marker, whose size will reflect the number of samples.
#' @param width The width (in inches) of the map image.
#' @param height The heigt (in inches) of the map image.
#'
#' @export
#'
#' @examples
#' #TBD
#
mapSampleCounts <- function (ctx, sampleSet,
                   aggregate="Province", minAggregateCount=1,
                   markerSize=c(4,40), colourBy="Province", showNames=TRUE,
                   width=15, height=15) {

    if (length(colourBy) > 1) {
        die ("colourBy parameter can only accet a single value (\"Country\" or \"Province\")")
    }
    colourAggLevel <- map.getAggregationLevelsFromLabels (colourBy)
    params <- list(
        map.aggregateCountMin=minAggregateCount,
        map.markerColourAggLevel=colourAggLevel,
        map.markerSize=markerSize,
        map.markerNames=showNames,
        map.size=c(width=width,height=height)
    )
    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/sampleCount", plotList=NULL,
                                 aggregation=aggLevels, measures=NULL, params=params)
}
#
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
#' @param aggregate TBD
#' @param markerScale TBD
#' @param showNames TBD
#' @param width TBD
#' @param height TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
mapBarcodeFrequencies <- function (ctx, sampleSet,
                   type=c("bar","pie"),
                   aggregate="Province", minAggregateCount=10, showNames=TRUE,
                   markerScale=0.8,
                   width=15, height=15) {

    mapParams <- list(
        map.aggregateCountMin=minAggregateCount,
        map.markerNames=showNames,
        map.haplo.visualizations=type,        
        map.haplo.markerScale=markerScale,
        map.size=c(width=width,height=height)
    )
    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/haploFreq", plotList=NULL,
                                 aggregation=aggLevels, measures=measures, params=mapParams)
}
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
#' @param clusteringMethod TBD
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
                   aggregate="Province", 
                   minAggregateCount=5, showNames=TRUE, markerScale=0.8,
                   similarityLevels=1.0, clusteringMethod="allNeighbours", minGroupSize=10,
                   width=15, height=15) {
                   
    mapParams <- list(
        cluster.identity.thresholds=similarityLevels,
        cluster.method=clusteringMethod,
        cluster.identity.minCount=minGroupSize,
        map.aggregateCountMin=minAggregateCount,
        map.markerNames=showNames,
        map.haplo.visualizations=type,        
        map.haplo.markerScale=markerScale,
        map.size=c(width=width,height=height)
    )
    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/haploShare", plotList=NULL,
                                 aggregation=aggLevels, measures=measures, params=mapParams)
                   
}

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
#' @param clusteringMethod TBD
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
                   aggregate="Province", 
                   minAggregateCount=5, showNames=TRUE, markerScale=0.8,
                   similarityLevels=1.0, clusteringMethod="allNeighbours", minGroupSize=10,
                   width=15, height=15) {

    mapParams <- list(
        cluster.identity.thresholds=similarityLevels,
        cluster.method=clusteringMethod,
        cluster.identity.minCount=minGroupSize,
        map.aggregateCountMin=minAggregateCount,
        map.markerNames=showNames,
        map.haplo.visualizations="group",        
        map.haplo.markerScale=markerScale,
        map.size=c(width=width,height=height)
    )
    aggLevels <- map.getAggregationLevelsFromLabels (aggregate)
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/haploShare", plotList=NULL,
                                 aggregation=aggLevels, measures=measures, params=mapParams)
}





##############################################################
#
#                 T O   B E   D O N E
#
##############################################################
# 
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

#         ############## PREV ##########################
#analysis.execute(ctx=ctx, datasetList=list(dataset.Laos), tasks=c("njt"), 
#                plotList=list(plot.def$byProvince, plot.def$byProvinceK1P1, plot.def$byK13, plot.def$byPm23))
#
#analysis.execute(ctx=ctx, datasetList=list(dataset.Laos), tasks=c("pca/PCoA"), 
#                plotList=list(plot.def$byProvince, plot.def$byProvinceK1P1, plot.def$byK13, plot.def$byPm23))
#         ############## PREV ##########################



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


#         ############## PREV ##########################
#haploNetParams <- list(
#    haploNet.minHaploCount=3
#)
#
#analysis.execute(ctx=ctx, datasetList=list(dataset.Laos), tasks="haploNet", params=haploNetParams, 
#		plotList=list(plot.def$byProvince, plot.def$byK13, plot.def$byPm23))
#         ############## PREV ##########################

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

#         ############## PREV ##########################
#graphParams <- list(
#    cluster.identity.thresholds=c(1.0,0.95,0.9,0.85,0.8,0.75),
#    cluster.identity.minCount=5,
#    graph.connectIdentityMin=0.7,
#    graph.layoutAlgorithm="fr",			#"drl", "fr","kk","lgl","mds","graphopt"
#    graph.weightFunction="identitySquared"
#)
#
#analysis.execute(ctx=ctx, datasetList=list(dataset.Laos), tasks="graph", params=graphParams, 
#                plotList=list(plot.def$byProvince, plot.def$byK13, plot.def$byPm23))
#         ############## PREV ##########################

