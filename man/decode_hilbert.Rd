% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hilbert.R
\name{decode_hilbert}
\alias{decode_hilbert}
\title{Given a Hilbert z-value, translate it into an x and y value.}
\usage{
decode_hilbert(z)
}
\arguments{
\item{z}{An integer to index a 2D matrix element.}
}
\value{
An array \code{c(x, y)} of array coordinates.
}
\description{
Given a Hilbert z-value, translate it into an x and y value.
}
