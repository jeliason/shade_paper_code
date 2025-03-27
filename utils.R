#' Perform 2D FFT-Based Density Smoothing (Simple Version)
#'
#' @param X A point pattern (`ppp` object) or density grid (`im` object)
#' @param kernel_func A function that computes kernel values dynamically
#' @param resolution Resolution of the density grid if `X` is a `ppp` object
#' @return A smoothed density estimate as an `im` object
smooth_density_fft <- function(X, kernel_func, resolution = 128) {
  # Convert point pattern to a density grid if needed
  if (spatstat.geom::is.ppp(X)) {
    X <- spatstat.geom::pixellate(X, dimyx = resolution, padzero = TRUE)
  }
  
  # Ensure input is an `im` object
  if (!spatstat.geom::is.im(X)) stop("X must be a `ppp` or `im` object")
  
  # Extract matrix and dimensions
  Y <- X$v
  nr <- nrow(Y)
  nc <- ncol(Y)
  
  # Pad the image to prevent wrap-around effects
  Ypad <- matrix(0, nrow = 2 * nr, ncol = 2 * nc)
  Ypad[1:nr, 1:nc] <- Y
  
  # Generate kernel values dynamically
  xcol.ker <- X$xstep * c(0:(nc - 1), -(nc:1))
  yrow.ker <- X$ystep * c(0:(nr - 1), -(nr:1))
  
  Kern <- outer(yrow.ker, xcol.ker, kernel_func)  # Compute dynamically
  
  # Compute FFT of image and kernel
  fft_Y <- stats::fft(Ypad)
  fft_Kern <- stats::fft(Kern)
  
  # Multiply in Fourier space and take inverse FFT
  smooth_Y <- Re(stats::fft(fft_Y * fft_Kern, inverse = TRUE)) / (4 * nc * nr)
  
  # Extract valid region
  smooth_Y <- smooth_Y[1:nr, 1:nc]
  
  # Convert back to `im` object
  smoothed_image <- spatstat.geom::im(smooth_Y, xcol = X$xcol, yrow = X$yrow, unitname = spatstat.geom::unitname(X))
  
  return(smoothed_image)
}

make_rbfs <- function(max_dist,
                      n_basis_functions=6,
                      basis_function_sigma=8) {
  gaussian_rbf <- function(x, mu, sigma) {
    exp(-(x - mu)^2 / (2 * sigma^2))
  }
  
  basis_function_centers <- seq(0, max_dist, length.out = n_basis_functions)  # Equally spaced centers
  
  rbfs <- lapply(basis_function_centers,\(mu) {
    function(x) {
      gaussian_rbf(x,mu,basis_function_sigma)
    }
  })
  
  names(rbfs) <- paste0("rbf",1:n_basis_functions)
  
  rbfs
}