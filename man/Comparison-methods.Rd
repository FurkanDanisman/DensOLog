\name{Comparison-methods}
\alias{BK2002}
\alias{BinnedNP}
\alias{KernSmooth}
\alias{evaluate_BK2002}
\alias{evaluate_BinnedNP}
\alias{evaluate_KernSmooth}
\title{Comparison density estimators}
\usage{
BK2002(counts, grid, h = NULL, bw_B = 50, bw_seed = 1,
  bw_method = c("ste", "dpi"), max_iter = 200, tol_L1 = 1e-6,
  verbose = FALSE)

BinnedNP(counts, grid, seed_start = 1L, max_tries = 200L)

KernSmooth(counts, grid, seed = 1, scalest = "minim", level = 2L,
  kernel = "normal", canonical = FALSE, truncate = TRUE)

evaluate_BK2002(object, x)

evaluate_BinnedNP(object, x, bandwidth = c("boot", "plugin"))

evaluate_KernSmooth(object, x)
}
\arguments{
\item{counts}{Integer observed frequencies for each bin.}
\item{grid}{Numeric bin breaks. Must have length \code{length(counts) + 1}.}
\item{h}{Optional bandwidth for BK2002.}
\item{bw_B, bw_seed, bw_method}{Bandwidth-selection controls for BK2002.}
\item{max_iter, tol_L1, verbose}{Fixed-point iteration controls for BK2002.}
\item{seed_start, max_tries}{Bandwidth retry controls for BinnedNP.}
\item{seed, scalest, level, kernel, canonical, truncate}{Controls passed to KernSmooth bandwidth and density routines.}
\item{object}{A fitted comparison-method object.}
\item{x}{Evaluation points.}
\item{bandwidth}{BinnedNP bandwidth choice.}
}
\description{
Fit and evaluate the three baseline estimators used in the DensOLog simulation
study: BK2002, BinnedNP, and KernSmooth. All methods use observed frequencies
and grid breaks as inputs.
}
\value{
Fitting functions return method-specific fitted objects. Evaluation functions
return numeric density values at \code{x}.
}
