% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biscale.R
\name{biscale.control}
\alias{biscale.control}
\title{Control for standard procedure}
\usage{
biscale.control(
  row.mean = FALSE,
  row.std = FALSE,
  col.mean = FALSE,
  col.std = FALSE
)
}
\arguments{
\item{row.mean}{if \code{row.mean = TRUE} (the default), row centering will be performed resulting in a matrix with row means zero.
If \code{row.mean} is a vector, it will be used in the iterative process.
If \code{row.mean = FALSE} nothing is done.}

\item{row.std}{if \code{row.std = TRUE} , row scaling will be performed resulting in a matrix with row variance one.
If \code{row.std} is a vector, it will be used in the iterative process.
If \code{row.std = FALSE} (the default) nothing is done.}

\item{col.mean}{similar to \code{row.mean}.}

\item{col.std}{similar to \code{row.std}.}
}
\value{
A list with components named as the arguments.
}
\description{
Various parameters that control aspects of the standard procedure.
}
