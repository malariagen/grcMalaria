###################################################################
# Task Parameters retrieval, with defaults
###################################################################
#
param.getParam <- function (paramName, paramList) {
    value <- param.defaultParams[[paramName]]
    if (!is.null(paramList)) {
        if (!is.null(paramList[[paramName]])) {
            value <- paramList[[paramName]]
        }
    }
    value
}
#
param.getParamIfDefined <- function (paramName, paramList) {
    if (paramName %in% names(paramList)) {
        return (param.getParam(paramName, paramList))
    }
    return (NULL)
}
#
param.defaultParams <- list (
    graph.connectIdentityMin=0.4,
    graph.layoutAlgorithm="fr",
    graph.weightPower=2,
    
    cluster.method="allNeighbours",
    cluster.minSize=5,
    cluster.identity.min=1.0,
    
    #haploNet.minHaploCount=1,
    
    map.connect.identity.min=1.0,
    map.connect.meanDistance.min=0.5,
    map.diversity.markerColours="red3",
    map.prevalence.markerColours="red3",
    
    # The number of samples that correspond to the "standard" size marker. Can be "mean", or an integer count
    map.cluster.markerSampleCount="mean",
    # The multiplier for sizing the marker wrt the "standard" size
    map.cluster.markerScale=1.0,
    
    map.aggregateCountMin=5,
    
    map.markerNames=FALSE,
    map.markerSize=16,
    plot.size=list(width=12,height=12)
)



