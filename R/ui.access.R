
UI.SUPPORTED_MAP_TYPES <- c("location", "sampleCount", "drug", "mutation", "alleleProp", "diversity")

ui.makeSingleMap <- function (ctx, sampleSet, type,
                   #timePeriods=NULL,	TBD
                   measure=NULL,
                   aggregate="Province", minAggregateCount=1,
                   markerSize=c(4,40), markerColours="red3", showNames=TRUE, colourBy="Province",
                   ...) {
                   
    if (!(type %in% UI.SUPPORTED_MAP_TYPES)) {
        stop (paste("Unsupported map type for UI:", type))
    }
                   
    # aggregate param 
    ui.checkSingleValue (aggregate, "aggregate")
    aggLevel <- map.getAggregationLevelsFromLabels (aggregate)

    # colourBy param (locations, sampleCount only)
    if (length(colourBy) > 1) {
        stop ("colourBy can only accet a single value (\"Country\" or \"Province\")")
    }
    colourAggLevel <- map.getAggregationLevelsFromLabels (colourBy)
    if (colourAggLevel > 1) {
        stop ("colourBy parameter can only accet values \"Country\" or \"Province\"")
    }
    
    # measure param
    if (type=="location") {
        measure <- "Location"
    } else if (type=="sampleCount") {
        measure <- "NumberOfSamples"
    }
    ui.checkSingleValue (measure, "measure")

    params <- list(
        analysis.timeIntervals=NULL,	# TBD
        analysis.measures=measure,
        aggregation.levels=aggLevel,
        map.aggregateCountMin=minAggregateCount,
        map.markerColourAggLevel=colourAggLevel,
        map.markerSize=markerSize,
        map.markerNames=showNames,
        map.diversity.markerColours=markerColours
    )
    params <- c(params, parsePlotParams(...))

    plotInfo <- execute.makeMapOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=paste0("map/", type), params=params)
    plotInfo
}
#
#
#
ui.checkSingleValue <- function (value, paramName) {
    if (length(value) > 1) {
        stop (paste(paramName, "can only accet a single value"))
    }
}
#
#
execute.makeMapOnSampleSet <- function (userCtx, sampleSetName, task, params) {
    #
    # Execute computations and plots 
    #
    execute.checkSampleSet (userCtx, sampleSetName)
    #
    # Execute the task
    #
    t <- execute.parseTask (task)
    if (t$name != "map") {
        stop (paste ("Invalid analysis task:", t$taskName))
    }
    plotInfo <- map.makeSingleMap (userCtx, sampleSetName, t$method, params)
    plotInfo
}
#
#
#
map.makeSingleMap <- function(userCtx, sampleSetName, mapType, params) {		#;print(mapType)	;print(params)
    #
    # Get the set of specs for the maps to be produced
    #
    mapSpecs <- map.createMapSpecs (userCtx, sampleSetName, mapType, params)
    mapSpec <- mapSpecs[[1]]
    #
    # Generate the ggplot2 plot object
    #
    mapPlot <- map.generateMapPlot (mapSpec, params)
    #
    # Add the "look and feel" elements to the map plot
    #
    mapPlot <- mapPlot + mapSpec$master$theme
    #
    # Deal with the legend - TBD
    #
    #mapPlot <- map.processLegend (mapPlot, mapSpec)
    plotInfo <- list(mainPlot=mapPlot, legendPlot=NULL, data=NULL)
    plotInfo
}
