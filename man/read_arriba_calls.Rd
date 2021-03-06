% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annoFuse-utils.R
\name{read_arriba_calls}
\alias{read_arriba_calls}
\title{Read in fusion calls from  Arriba v1.1.0}
\usage{
read_arriba_calls(arriba_calls)
}
\arguments{
\item{arriba_calls}{Please refer to software documenation [fusions.tsv](https://arriba.readthedocs.io/en/latest/output-files/)}
}
\value{
A data.frame object with correct column specifications
}
\description{
Read in fusion calls from  Arriba v1.1.0
}
\examples{
fusionfileArriba <- read_arriba_calls(
  system.file("extdata", "arriba_example.tsv", package = "annoFuseData")
)
}
