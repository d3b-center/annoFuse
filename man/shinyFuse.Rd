% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shinyFuse.R
\name{shinyFuse}
\alias{shinyFuse}
\title{shinyFuse}
\usage{
shinyFuse(out_annofuse = NULL)
}
\arguments{
\item{out_annofuse}{The character string specifying the location of the file
output by the annoFuse pipeline. This file needs to be structured with the set
of columns required for the later exploration steps in the interactive app.}
}
\value{
A Shiny app object
}
\description{
Exploring interactively the results of the annoFuse pipeline
}
\details{
The application can also be started without specifying the location of the data,
which can be provided (via upload) at runtime.
}
\examples{
out_annofuse <- system.file("extdata", "PutativeDriverAnnoFuse.tsv",
  package = "annoFuseData"
)
if (interactive()) {
  shinyFuse(out_annofuse)
}
}
