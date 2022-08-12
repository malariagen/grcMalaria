
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
#' @param sampleSetName The name of the sample set, to be used to identifty it when calling analysis tasks.
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
#' Map of Locations
#'
#' Creates a map showing markers for all sites/districts/provinces where samples are collected.
#' For each location, we place a marker on the map, coloured according to the administrative
#' division it is in
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used; must have been previusly created by selectSampleSet().
#' @param aggregate The administrative level at which we aggregate (Province/District/Site)
#' @param markerSize Allows adjustment of the size of markers on the map. If only one value is passed, 
#'                   all markers will be that size; if two values are passed, they will be used as the min 
#'                   and max sizeof the marker, whose size will reflect the number of samples.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param colourBy Shows the aggregation level to be used to colour the markers (Country or Province)
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @export
#'
#' @examples
#' #TBD
#
mapLocations <- function (ctx, sampleSet,
                   aggregate="Province",
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
    params <- list(
        aggegation.levels=aggLevels,
        map.markerColourAggLevel=colourAggLevel,
        map.markerSize=markerSize,
        map.markerNames=showNames,
        plot.size=list(width=width,height=height)
    )
    analysis.executeOnSampleSet (ctx=ctx, sampleSetName=sampleSet, tasks="map/location", params=params)
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
#' @param timePeriods The list of time period object for partitioning samples into time-interval plots
#' @param aggregate The administrative level at which we aggregate (Province or District)
#' @param minAggregateCount The minimum count of aggregated samples, below which a marker is not shown.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param colourBy Shows the aggregation level to be used to colour the markers (Country or Province)
#' @param markerSize Allows adjustment of the size of markers on the map. If only one value is passed, 
#'                   all markers will be that size; if two values are passed, they will be used as the min 
#'                   and max sizeof the marker, whose size will reflect the number of samples.
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @export
#'
#' @examples
#' #TBD
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
#' @param timePeriods The list of time period object for partitioning samples into time-interval plots
#' @param drugs An array of drugs for which prevalence of resistance will be estimated; "ALL" creats maps for all the drugs for which we have resistance predictions.
#' @param aggregate The administrative level at which we aggregate
#' @param minAggregateCount The minimum count of aggregated samples, below which a marker is not shown.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param markerSize Allows adjustment of the size of markers on the map.
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @export
#'
#' @examples
#' #TBD
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
#' Map the prevalence of genetic mutations.
#'
#' @param ctx TBD
#' @param sampleSet TBD
#' @param timePeriods The list of time period object for partitioning samples into time-interval plots
#' @param mutations TBD
#' @param aggregate TBD
#' @param minAggregateCount TBD
#' @param showNames TBD
#' @param markerSize TBD
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
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
#' @param ctx TBD
#' @param sampleSet TBD
#' @param timePeriods The list of time period object for partitioning samples into time-interval plots
#' @param measures can be "ALL", or any vector containing one or more of ("maxHaploFreq","haploHet", "meanSnpHet","medianDistance")
#' @param aggregate TBD
#' @param minAggregateCount TBD
#' @param showNames TBD
#' @param markerColours TBD
#' @param markerSize TBD
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
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
#' @param ctx TBD
#' @param sampleSet TBD
#' @param measures TBD
#' @param minIdentity TBD 
#' @param meanDistanceLevels TBD
#' @param aggregate TBD
#' @param minAggregateCount TBD
#' @param showNames TBD
#' @param markerSize TBD
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
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
#' @param ctx TBD
#' @param sampleSet TBD
#' @param type TBD
#' @param aggregate TBD
#' @param minAggregateCount TBD
#' @param showNames TBD
#' @param markerScale TBD
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
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
#' Clustering
#'
#' @param ctx TBD
#' @param sampleSet TBD
#' @param clusterSet TBD
#' @param minIdentity TBD
#' @param clusteringMethod TBD
#' @param minClusterSize TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
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
#' Plot graphs of bracode identity, using cluster
#'
#' @param ctx TBD
#' @param sampleSet TBD
#' @param clusterSet TBD
#' @param graphLayout TBD
#' @param weightPower TBD
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
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
#' Map cluster sharing
#'
#' @param ctx TBD
#' @param sampleSet TBD
#' @param clusterSet TBD
#' @param type TBD
#' @param aggregate TBD
#' @param minAggregateCount TBD
#' @param showNames TBD
#' @param markerScale TBD
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
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
#' Map Cluster prevalence at different sites, connecting these sites
#'
#' @param ctx TBD
#' @param sampleSet TBD
#' @param clusterSet TBD
#' @param aggregate TBD
#' @param minAggregateCount TBD
#' @param showNames TBD
#' @param markerScale TBD
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
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
#' Graphical Attributes
#'
#' @param ctx TBD
#' @param name TBD
#' @param field TBD
#' @param file TBD
#' @param sheet TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
loadGraphicAttributes <- function (ctx, name, field, file, sheet) {
    ctx <- graphics.loadAttributes (ctx, name, field, file, sheet)
    ctx
}

#############################################################
#
#' Plot Population structure (PCA)
#'
#' @param ctx TBD
#' @param sampleSet TBD TBD
#' @param type TBD
#' @param plots TBD
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
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
#' @param ctx TBD
#' @param sampleSet TBD TBD
#' @param type TBD
#' @param plots TBD
#' @param width The width (in inches) of the map image.
#' @param height The height (in inches) of the map image.
#'
#' @return TBD
#' @export
#'
#' @examples
#'  #TBD
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
