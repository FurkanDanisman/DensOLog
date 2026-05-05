L2_DensOLog <- function(counts, grid, pdf, smooth = FALSE, length_grid = 1001, ...) {
  fit <- DensOLog(counts, grid, smooth = smooth, ...)
  eval_grid <- .eval_grid_from_breaks(grid, length_grid = length_grid)
  dx <- eval_grid[2L] - eval_grid[1L]
  sqrt(sum((pdf(eval_grid) - ddensolog(fit, eval_grid, smooth = smooth))^2) * dx)
}

L2_BK2002 <- function(counts, grid, pdf, length_grid = 1001, ...) {
  fit <- BK2002(counts, grid, ...)
  eval_grid <- .eval_grid_from_breaks(grid, length_grid = length_grid)
  dx <- eval_grid[2L] - eval_grid[1L]
  sqrt(sum((pdf(eval_grid) - evaluate_BK2002(fit, eval_grid))^2) * dx)
}

L2_BinnedNP <- function(counts, grid, pdf, length_grid = 1001, ...) {
  fit <- BinnedNP(counts, grid, ...)
  eval_grid <- .eval_grid_from_breaks(grid, length_grid = length_grid)
  dx <- eval_grid[2L] - eval_grid[1L]
  sqrt(sum((pdf(eval_grid) - evaluate_BinnedNP(fit, eval_grid))^2) * dx)
}

L2_KernSmooth <- function(counts, grid, pdf, length_grid = 1001, ...) {
  fit <- KernSmooth(counts, grid, ...)
  eval_grid <- .eval_grid_from_breaks(grid, length_grid = length_grid)
  dx <- eval_grid[2L] - eval_grid[1L]
  sqrt(sum((pdf(eval_grid) - evaluate_KernSmooth(fit, eval_grid))^2) * dx)
}
