# Method Questions

## Multivariate Support and Non-Uniform Grids

This note records what can be directly supported from the original papers or the package documentation used for the comparison methods.

| Method | Multivariate? | Non-uniform observed grid? | Direct support |
|---|---:|---:|---|
| `KernSmooth::bkde` | No, univariate | No for the output grid used by `bkde` | The `bkde` documentation defines `x` as a numeric vector and `gridsize` as the number of equally spaced points where the estimate is computed. |
| `KernSmooth::bkde2D` | Yes, but only two-dimensional | No for the output mesh | The `bkde2D` documentation defines `x` as a two-column numeric matrix and says `gridsize` gives the number of equally spaced points in each direction. |
| `binnednp` | No, essentially univariate interval-grouped data | Yes | Barreiro-Ures et al. describe the target as the distribution and density functions of "the random variable" relating weed emergence to hydrothermal time. The `bw.dens.binned` interface defines `y` as a vector of observed values giving the extremes of the interval sequence. |
| Blower and Kelsall (2002) / NKDEBD | Yes, in the paper | Yes, in the paper | Blower and Kelsall state that the same method works for binned data in `R^d`, for any finite collection of bounded and abutting bins. In their spatial example, they further state that bins need not be rectangular and may be arbitrarily defined adjoining geographic areas. |

## Suggested Wording

`KernSmooth::bkde()` is a univariate binned kernel density estimator whose output grid is equally spaced. The `KernSmooth` package also includes `bkde2D()` for two-dimensional data, again on an equally spaced output mesh.

`binnednp` is designed for one-dimensional interval-grouped data and supports arbitrary interval endpoints through its endpoint vector.

For Blower and Kelsall's nonlinear KDE for binned data, the paper supports a stronger statement: the method is theoretically presented beyond one dimension and beyond uniform rectangular grids. The current implementation in this repository is more restrictive than the paper because `nonlinear_kde_binned_BK2002()` uses a one-dimensional numeric `breaks` vector and validates an approximately uniform computational grid.

## Supporting Details

### `binnednp`

Evidence that `binnednp` is univariate:

- The PubMed abstract for Barreiro-Ures et al. (2019) describes the target as the distribution and density functions of "the random variable" that relates weed emergence with hydrothermal time.
- The `bw.dens.binned` documentation defines `y` as a vector of observed values that define the extremes of the sequence of intervals in which data are binned.
- The examples in the `bw.dens.binned` documentation generate a one-dimensional complete sample with `x <- rnorm(n, 6, 1)` and define one-dimensional interval endpoints with `y <- seq(...)`.

Evidence that `binnednp` supports non-uniform interval endpoints:

- The argument `y` is documented as a vector defining interval extremes, not as a fixed-width grid.
- The method uses proportions or counts within each interval (`w` or `ni`) paired with the endpoint vector `y`, so unequal interval widths are representable.

### Blower and Kelsall (2002)

Evidence for multivariate support:

- In the list of properties of the estimator, the paper states that the same method works for binned data in `R^d`, given any finite collection of bounded and abutting bins.
- Later, the paper says the results of Sections 3, 4, and 5 extend to bounded and abutting rectangular regions in higher dimensions.
- The implementation section includes a two-dimensional spatial-data example.

Evidence for non-uniform or arbitrary bins:

- The one-dimensional introduction uses consecutive, abutting, bounded intervals in `R`; it does not require equal interval widths.
- The property list uses "any finite collection" of bounded and abutting bins.
- In the spatial example, the paper says there is no need for the bins to be rectangular and that the method has potential applicability to counts over arbitrarily defined adjoining geographical areas.

## References

- Wand, M. P. and Jones, M. C. (1995). *Kernel Smoothing*. Chapman and Hall.
- Ripley, B. D. `KernSmooth`: functions for kernel smoothing supporting Wand and Jones (1995). R package documentation for `bkde` and `bkde2D`.
- Barreiro-Ures, D., Francisco-Fernandez, M., Cao, R., Fraguela, B. B., Doallo, R., Gonzalez-Andujar, J. L., and Reyes, M. (2019). Analysis of interval-grouped data in weed science: The `binnednp` Rcpp package.
- Blower, G. and Kelsall, J. E. (2002). Nonlinear kernel density estimation for binned data: convergence in entropy. *Bernoulli*, 8(4), 423-449.
