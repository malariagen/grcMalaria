###################################################################
# Task Parameters retrieval, with defaults
###################################################################
#
analysis.getParam <- function (paramName, paramList) {
    value <- analysis.defaultParams[[paramName]]
    if (!is.null(paramList)) {
        if (!is.null(paramList[[paramName]])) {
            value <- paramList[[paramName]]
        }
    }
    value
}

analysis.defaultParams <- list (
    graph.connectIdentityMin=0.4,
    graph.layoutAlgorithm="fr",
    graph.weightFunction="identitySquared",
    
    cluster.method="allNeighbours",
    cluster.identity.minCount=5,
    cluster.identity.thresholds=c(1.0),
    
    haploNet.minHaploCount=1,
    
    map.connect.similarity.min=1.0,
    map.connect.meanDistance.min=0.5,
    map.diversity.markerColours="red3",
    map.prevalence.markerColours="red3",
    
    # The number of samples that correspond to the "standard" size marker. Can be "mean", or an integer count
    map.haplo.markerSampleCount="mean",
    # The multiplier for sizing the marker wrt the "standard" size
    map.haplo.markerScale=1.0,
    
    map.aggregateCountMin=5,
    
    map.markerNames=FALSE,
    map.markerSize=16,
    map.size=c(12,12)
)



