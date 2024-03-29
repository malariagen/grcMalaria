% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public.R
\name{setColourPalette}
\alias{setColourPalette}
\title{Set to a user-defined colour palette}
\usage{
setColourPalette(ctx, palette)
}
\arguments{
\item{ctx}{The analysis context, created by intializeContext().}

\item{palette}{The desired palette, specified as separate HTML colour codes}
}
\description{
Users can define their own colour palette instead of using the default colours.
This function needs to be executed before sample selection, and after running initializeContext(), otherwise the defined colour palette will not be applied.
}
\examples{
\dontrun{
# To view default 25-colour palette
ctx$config$defaultPalette

# Set colours
setColourPalette(ctx, palette=c("red","lightblue","#93c4d2"))

# then select sample set
selectSampleSet(ctx, sampleSet="SamplesetName",
                select=list(
                           list(field="Species", values=c("Pf"))
                ))

# then proceed with the analysis
}
}
