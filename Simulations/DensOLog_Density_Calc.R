L2_Distance_calc = function(x,grid,pdf,range){
  
  # Running the LCD Algorithm 
  
  fhat2 <- get_fhatn(x, grid, alpha=2)
  x_obs <- fhat2$fhatn$x  
  delta = max(diff(grid))
  rng <- range(x_obs, na.rm = TRUE)
  grid_narrow <- grid[ grid >= rng[1] - delta & grid <= rng[2] + delta]
  # grid_narrow = grid
  true_density_values <- pdf(grid_narrow)
  estimated_grid_density = evaluateLogConDens(grid_narrow, fhat2$fhatn, which=4)[,5]

  # Compute the absolute difference for L2 distance
  sqr_differance <- (true_density_values - estimated_grid_density)^2
  L2 = sum(sqr_differance*delta)
  
  L2 = L2 + integral(function(u) pdf(u)^2, range[1], grid_narrow[1])
  L2 = L2 + integral(function(u) pdf(u)^2, grid_narrow[length(grid_narrow)], range[2])
  L2 = sqrt(L2)
  
  return(L2)
  
}


L2_binnednp = function(x,grid,pdf,range){
  
  # Running the LCD Algorithm 
  
  estimated_grid_density <- binned_np.func(x, grid)  
  true_density_values <- pdf(grid)
  delta = max(diff(grid))
  
  # Compute the absolute difference for L2 distance
  sqr_differance <- (true_density_values - estimated_grid_density)^2
  L2 = sum(sqr_differance*delta)
  
  #L2 = integral(pdf,range[1],grid[1])^2 + L2 
  #L2 = integral(pdf,grid[length(grid)],range[2])^2 + L2
  #L2 = sqrt(L2)
  
  L2 = L2 + integral(function(u) pdf(u)^2, range[1], grid[1])
  L2 = L2 + integral(function(u) pdf(u)^2, grid[length(grid)], range[2])
  L2 = sqrt(L2)
  
  return(L2)
  
}

L2_binnednp2 = function(x,grid,pdf,range){
  
  # Running the LCD Algorithm 
  
  estimated_grid_density <- binned_np.func2(x, grid)  
  true_density_values <- pdf(grid)
  delta = max(diff(grid))
  
  # Compute the absolute difference for L2 distance
  sqr_differance <- (true_density_values - estimated_grid_density)^2
  L2 = sum(sqr_differance*delta)
  
  L2 = L2 + integral(function(u) pdf(u)^2, range[1], grid[1])
  L2 = L2 + integral(function(u) pdf(u)^2, grid[length(grid)], range[2])
  L2 = sqrt(L2)
  
  return(L2)
  
}

L2_from_res <- function(x,grid, pdf, range) {
  
  hobj <- hist(x, breaks = grid, plot = FALSE)
  
  res <- nonlinear_kde_binned_BK2002(
    counts = hobj$counts,
    breaks = hobj$breaks,
    grid = grid,
    h = NULL,        # <- forces SJ-from-binned as suggested in the paper
    bw_B = 50,
    bw_method = "ste"
  )
  
  x <- res$x
  fhat <- res$fhat
  dx <- x[2] - x[1]   # assumes equally spaced grid (as in our implementation)
  
  true <- pdf(x)
  L2_mid <- sum((true - fhat)^2) * dx
  
  # tails (L2 uses integral of pdf^2 outside grid)
  L2_left  <- integral(function(u) pdf(u)^2, range[1], x[1])
  L2_right <- integral(function(u) pdf(u)^2, x[length(x)], range[2])
  
  sqrt(L2_left + L2_mid + L2_right)
}

L2_from_kernsmooth <- function(x,grid, pdf, range) {
  
  hobj <- hist(x, breaks=grid, plot=FALSE)
  counts <- hobj$counts
  centers <- (grid[-1] + grid[-length(grid)]) / 2
  
  y <- rep(centers, times=counts)  # expand at centers (still discrete, but correct locations)

  # 1) Bandwidth via dpik (your settings)
  
  h <- dpik_with_fallback_jitter(
    y = y,
    counts = hobj$counts,
    breaks = hobj$breaks,
    grid = grid,
    seed = 1
  )
  
  #h <- dpik(y, scalest="minim", level=2L, kernel="normal",
  #          canonical=FALSE, gridsize=length(grid), range.x=range(grid), truncate=TRUE)
  
  
  # 2) KDE via bkde (your settings)
  fit <- bkde(y, bandwidth = h,
              kernel = "normal", canonical = FALSE,
              gridsize = length(grid), range.x = range(grid),
              truncate = TRUE)
  
  x   <- fit$x
  fhat <- fit$y
  dx  <- x[2] - x[1]   # bkde grid is equally spaced with gridsize
  
  # 3) L2 on bkde grid
  true <- pdf(x)
  L2_mid <- sum((true - fhat)^2) * dx
  
  # 4) tails (truncate=TRUE => estimate is 0 outside [x[1], x[end]])
  L2_left  <- integral(function(u) pdf(u)^2, range[1], x[1])
  L2_right <- integral(function(u) pdf(u)^2, x[length(x)], range[2])
  
  sqrt(L2_left + L2_mid + L2_right)
  
}

L1_Distance_calc = function(x,grid,pdf,range){
  
  # Running the LCD Algorithm 
  
  fhat2 <- get_fhatn(x, grid, alpha=2)  
  
  true_density_values <- pdf(grid)
  estimated_grid_density = evaluateLogConDens(grid, fhat2$fhatn, which=4)[,5]
  delta = max(diff(grid))
  
  
  # Compute the absolute difference for L1 distance
  abs_difference <- abs(true_density_values - estimated_grid_density)
  L1 <- sum(abs_difference*delta)
  
  L1 = integral(pdf,range[1],grid[1]) + L1
  L1 = integral(pdf,grid[length(grid)],range[2]) + L1
  
  return(L1)
  
}


