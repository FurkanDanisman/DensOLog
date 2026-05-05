\name{L2-functions}
\alias{L2_DensOLog}
\alias{L2_BK2002}
\alias{L2_BinnedNP}
\alias{L2_KernSmooth}
\title{L2 error calculations}
\usage{
L2_DensOLog(counts, grid, pdf, smooth = FALSE, length_grid = 1001, ...)

L2_BK2002(counts, grid, pdf, length_grid = 1001, ...)

L2_BinnedNP(counts, grid, pdf, length_grid = 1001, ...)

L2_KernSmooth(counts, grid, pdf, length_grid = 1001, ...)
}
\arguments{
\item{counts}{Integer observed frequencies for each bin.}
\item{grid}{Numeric bin breaks.}
\item{pdf}{True density function used for comparison.}
\item{smooth}{Logical; use smoothed DensOLog when \code{TRUE}.}
\item{length_grid}{Number of evaluation grid points.}
\item{...}{Additional arguments passed to the fitting function.}
}
\description{
Compute numerical L2 errors on a dense evaluation grid, matching the simulation
study convention.
}
\value{A numeric L2 error.}
