% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public.R
\name{plotClusterGraph}
\alias{plotClusterGraph}
\title{Plot graphs of barcode identity, using cluster}
\usage{
plotClusterGraph(
  ctx,
  sampleSet,
  clusterSet,
  graphLayout = "fr",
  weightPower = 2,
  ...
)
}
\arguments{
\item{ctx}{The analysis context, created by intializeContext()}

\item{sampleSet}{The name of the sample set being used, which must have been previously created by selectSampleSet()}

\item{clusterSet}{The name of the clustering set, defined in findClusters()}

\item{graphLayout}{TBD, the default="fr"}

\item{weightPower}{TBD, the default is 2}

\item{...}{Aesthetics: any of the following plot parameters,
width: the width of the plot (numeric), default=15
height: the height of the plot (numeric), default=15
units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
dpi: the resolution of the plot output, expressed as dots per inch, default=300
format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
legendWidth: specifies how wide a fixed width legend space should be, default="NULL"
legendDirection: specifies location of legend. Supported values: "vertical" (default), "horizontal"
legendFontSize: specifies font size of the legend (numeric), default=4
axisTitleSize: specifies axis label font size (numeric), default=1}
}
\value{
Produces network graphs of clusters assigned in findClusters function.
}
\description{
Creates graphs of cluster networks within a sampleset of Plasmodium samples.
The graphs are based on a cluster analysis, which needs to be produced with the findClusters function in advance.
The findClusters function partitions the Plasmodium samples into clusters with a similar genetic background, based on the 'genetic barcode' similarity.
}
\examples{
\dontrun{
#Plot cluster network graph for sampleset "Laos", and clusterset "LAclust"
plotClusterGraph (ctx, sampleSet="Laos", clusterSet="LAclust",
                  graphLayout="fr", weightPower=2)
}
}
