% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public.R
\name{mapSampleCounts}
\alias{mapSampleCounts}
\title{Map of Sample Counts}
\usage{
mapSampleCounts(
  ctx,
  sampleSet,
  timePeriods = NULL,
  aggregate = "Province",
  minAggregateCount = 1,
  markerSize = c(4, 40),
  markerFontSize = 6,
  colourBy = "Province",
  showNames = TRUE,
  nameFontSize = 5,
  ...
)
}
\arguments{
\item{ctx}{The analysis context, created by intializeContext().}

\item{sampleSet}{The name of the sample set being used, which must have been previously created by selectSampleSet().}

\item{timePeriods}{Time-sequence maps can be implemented using a parameter called timePeriods parameter. When this is passed to a function, it will partition samples into time-interval plots. Time intervals parameters are available in 5 analyses: mapSampleCounts(), mapDrugResistancePrevalence(), mapMutationPrevalence(), mapAlleleProportions() and mapDiversity().
Why use time-intervals parameter:
It will produce maps with consistent geographical boundaries
Users can slice up the dataset in any specified time period}

\item{aggregate}{The administrative level at which we aggregate (Province or District). Separate maps are created for each administrative level.}

\item{minAggregateCount}{The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.}

\item{markerSize}{Allows adjustment of the size of markers on the map. If only one value is passed, all markers will be of that size; if two values are passed, they will be used as the min and max size of the marker, whose size will reflect the number of samples.}

\item{markerFontSize}{Allows adjustment of the font size shown on the markers (numeric), default=6}

\item{colourBy}{Shows the aggregation level to be used to colour the markers (”Province” or “Country”). For each aggregation unit, we place a marker on the map, coloured according to the level of resistance to the drug or mutation, with a label indicating the prevalence.}

\item{showNames}{If TRUE, labels are shown with the name of the aggregation unit (Province or District)}

\item{nameFontSize}{Allows adjustment of the Province or District label name font size (numeric), default=5}

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
Creates maps showing the numbers of samples collected at different aggregation levels.
Number of samples are shown on the marker.
If number of samples is less than minimum count of aggregated samples (minAggregateCount), the marker will not appear.
File name with 'filtered’ is showing samples that passed the quality filtering threshold (default is minSampleTypability=0.75, meaning samples that have more than 25% missingness in their barcode are filtered out and not shown on the map).
File name with ‘unfiltered’ is showing all the samples without applying quality filtering.
}
\examples{
\dontrun{
## To implement time-sequence maps ##
# 1. Specify time periods
# name: the name parameter is required and is added to the end of the file name of the produced 
# files.
# There are two types of time interval: type="year" uses start date and by default show until 
# the end of the defined year; type="period" allows user to specify any time period using both 
# the start and end date.
# start and end - start and end date. The parameter start must be a date in the 
# format dd-MMM-yyyy.
# NOTE: Public GRC data files contains only year of collection, therefore only the year in the 
# provided “start” date is used to produce time-sequence maps.
periods <- list(
               list(name="2021", type="year", start="1-Jan-2021"),
               list(name="2020", type="year", start="1-Jan-2020"),
               list(name="2017-19", type="period", start="1-Jan-2017", end="31-Dec-2019")
            )


# 2. Apply timePeriods parameter to mapSampleCounts
mapSampleCounts(ctx, sampleSet="Laos", timePeriods = periods,
                aggregate=c("Province","District"), minAggregateCount=1,
                markerSize=c(10,40), colourBy="Province", showNames=TRUE,
                width=15, height=15)
}
}
