.check_counts_grid <- function(counts, grid) {
  if (!is.numeric(counts) || !is.numeric(grid)) {
    stop("counts and grid must be numeric.", call. = FALSE)
  }
  if (length(grid) != length(counts) + 1L) {
    stop("grid must contain bin breaks, so length(grid) must equal length(counts) + 1.", call. = FALSE)
  }
  if (any(!is.finite(grid)) || any(diff(grid) <= 0)) {
    stop("grid must be finite and strictly increasing.", call. = FALSE)
  }
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("counts must be finite and non-negative.", call. = FALSE)
  }
  if (sum(counts) <= 0) {
    stop("counts must contain at least one observation.", call. = FALSE)
  }
  if (any(abs(counts - round(counts)) > sqrt(.Machine$double.eps))) {
    stop("counts must be observed frequencies, so they must be integer-valued.", call. = FALSE)
  }
  invisible(TRUE)
}

.check_uniform_grid <- function(grid, rtol = 1e-8, atol = 1e-12) {
  d <- diff(grid)
  dx <- stats::median(d)
  dev <- max(abs(d - dx))
  if (dev > atol + rtol * max(1, abs(dx))) {
    stop("grid must be uniformly spaced for this implementation.", call. = FALSE)
  }
  invisible(TRUE)
}

.eval_grid_from_breaks <- function(grid, length_grid = 1001) {
  seq(min(grid), max(grid), length.out = length_grid)
}

.left_centers <- function(grid) {
  (grid[-length(grid)] + grid[-1L]) / 2
}

.restore_seed <- function(seed) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  old_seed
}

.on_exit_restore_seed <- function(old_seed) {
  if (is.null(old_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }
}
