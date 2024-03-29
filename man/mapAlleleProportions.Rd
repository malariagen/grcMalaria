% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public.R
\name{mapAlleleProportions}
\alias{mapAlleleProportions}
\title{Map the proportion of alleles for a given mutation (including Pfkelch13 and gene amplifications).}
\usage{
mapAlleleProportions(
  ctx,
  sampleSet,
  timePeriods = NULL,
  mutations = "ALL",
  aggregate = "Province",
  minAggregateCount = 10,
  markerSize = 16,
  alleleColours = NULL,
  showNames = TRUE,
  nameFontSize = 5,
  ...
)
}
\arguments{
\item{ctx}{The analysis context, created by intializeContext()}

\item{sampleSet}{The name of the sample set being used; must have been previously created by selectSampleSet()}

\item{timePeriods}{Time-sequence maps can be implemented using a parameter called timePeriods parameter. When this is passed to a function, it will partition samples into time-interval plots. Time intervals parameters are available in 5 analyses: mapSampleCounts(), mapDrugResistancePrevalence(), mapMutationPrevalence(), mapAlleleProportions() and mapDiversity().
Why use time-intervals parameter:
It will produce maps with consistent geographical boundaries
Users can slice up the dataset in any specified time period}

\item{mutations}{The genetic mutation types for which prevalence will be calculated; "ALL" creates maps for all 11 available genetic mutations
Available genetic mutations include: "Pfkelch13", "pm23-Amp", "mdr1-Amp", "crt", "dhfr", "dhps", "mdr1", "mdr2", "arps10", "fd", "exo"
To specify a mutation put it in between quotation marks e.g. "Pfkelch13", or  c("Pfkelch13","dhfr","fd") to select several genetic mutation types}

\item{aggregate}{The administrative level at which we aggregate. Separate maps are created for each administrative level.}

\item{minAggregateCount}{The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.}

\item{markerSize}{Allows adjustment of the size of markers on the map. If only one value is passed, all markers will be of that size; if two values are passed, they will be used as the min and max size of the marker, whose size will reflect the number of samples.}

\item{alleleColours}{A vector of R colours (specified either as R colour names, or HTML colour codes) which is used to colour the slices in the pie symbols of the plot.
All elements of the vector must be named with the name of the value (allele) they will represent.
The sequence of names will determine the order in which the values are displayed.
The name "Other" will indicate the colour used to group together samples that carry alleles other than those names. Usually, "Other" will be listed at the end of the sequence.
It no colour is specified for "Other" in the sequence, then an additional "Other"="white" element will be automatically added to the end of the sequence.}

\item{showNames}{If TRUE, labels are shown with the name of the aggregation unit (Province or District).}

\item{nameFontSize}{Allows adjustment of the Province or District label name font size (numeric), default=5.}

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
\description{
This function creates maps with pie-charts showing the proportion of selected alleles in each administrative division.
Markers are displayed as pie-chart showing proportion of all the homozygous alleles for the specified gene.
Missing genotypes and heterozygous genotypes are excluded in the calculations.
}
\examples{
\dontrun{
# Produce maps for all available genetic mutations for sample set called Cambodia,
# at province level, with a minimum of 10 samples
mapAlleleProportions (ctx, sampleSet="Cambodia",
                      mutations="ALL",
                      aggregate="Province", minAggregateCount=10,
                      showNames=TRUE, markerSize=c(10,25))

# Produce maps of some of the most common Pfkelch13 variants (WT, C580Y, R539T and Y493H),
# and group together the other alleles
kelch13Col <- c(WT="green",C580Y="red",R539T="orange",Y493H="yellow",Other="gray")
mapAlleleProportions (ctx, sampleSet="Cambodia",
                      mutations="Pfkelch13", alleleColours=kelch13Col,
                      aggregate="District", minAggregateCount=10,
                      showNames=TRUE, markerSize=c(10,25))
}
}
\seealso{
Useful links:
\itemize{
 \item user guide for information on functionality of specific genetic mutation types {\url{https://genremekong.org/tools/grcmalaria-r-package-user-guide}}
}
}
