% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadGEO.R
\name{loadCounts}
\alias{loadCounts}
\title{Loads expression data from .h5 count files.
Only samples with counted expression are kept.
If es already containts expression data it is returned as is.}
\usage{
loadCounts(es, counts_dir)
}
\arguments{
\item{es}{ExpressionSet from GEO to check for expression in ARCHS4/dee2 or other h5 files}

\item{counts_dir}{directory with  .h5 files  collections. There must be meta.rda file
in counts_dir and each collection's sub directory must have meta.txt file with description.
Also \code{counts_dir} must contain \code{counts_priority.txt} file.}
}
\value{
either original es or an ExpressionSet with loaded count data from ARCHS4
}
\description{
Loads expression data from .h5 count files.
Only samples with counted expression are kept.
If es already containts expression data it is returned as is.
}
\keyword{internal}
