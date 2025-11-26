# Global constants for spatial interaction analysis
# Minimum interaction radius (microns) - excludes very short distances where interpretation
# is ambiguous due to potential cell overlap or segmentation artifacts.
# Set to 25 microns (~1.5-2 cell diameters) to ensure clear intercellular spacing.
MIN_INTERACTION_RADIUS <- 30

#' Get data path for simulation/analysis
#'
#' Returns the appropriate data path based on environment (HPC vs local)
#' and creates the directory if it doesn't exist.
#'
#' @param sim_dir Simulation directory name (e.g., "sim_timing", "crc_analysis")
#' @param local_path Local path to use when not on HPC (default: "./data/")
#' @return Character string with full path to data directory
get_data_path <- function(sim_dir, local_path = "./data/") {
  system_env <- Sys.getenv("SYSTEM_ENV", "")

  if (system_env == "HPC") {
    # Use HPC_DATA_PATH from environment
    base_path <- Sys.getenv("HPC_DATA_PATH", "./data")
    path <- file.path(base_path, sim_dir, "")  # Empty string forces trailing slash
  } else {
    # Local execution - use simulation-specific directory
    path <- file.path(".", sim_dir, "data", "")  # Empty string forces trailing slash
  }

  # Create directory if it doesn't exist
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  return(path)
}

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

#' Create truncated radial basis functions
#'
#' Creates Gaussian RBFs that are zeroed out below a minimum distance.
#' This encodes prior knowledge that very short distances are unreliable.
#'
#' @param n_basis_functions Number of basis functions
#' @param max_dist Maximum distance for RBF support
#' @param basis_function_sigma Spread of Gaussian RBFs
#' @param min_dist Minimum distance - basis functions are zero below this
#' @return List of basis functions
make_truncated_rbfs <- function(n_basis_functions,
                                max_dist,
                                basis_function_sigma,
                                min_dist = MIN_INTERACTION_RADIUS) {

  # Create standard Gaussian RBF
  gaussian_rbf <- function(x, mu, sigma) {
    exp(-(x - mu)^2 / (2 * sigma^2))
  }

  # Equally spaced centers from min_dist to max_dist
  basis_centers <- seq(min_dist, max_dist, length.out = n_basis_functions)

  rbfs <- lapply(basis_centers, function(mu) {
    function(x) {
      # Zero out below minimum distance
      val <- ifelse(x < min_dist, 0, gaussian_rbf(x, mu, basis_function_sigma))
      return(val)
    }
  })

  names(rbfs) <- paste0("rbf", 1:n_basis_functions)

  return(rbfs)
}

#' Compute Simultaneous Credible Bands for Spatial Interaction Curves
#'
#' Computes simultaneous credible bands that control the family-wise error rate
#' across all distances, using the supremum of standardized deviations approach.
#' This ensures the entire curve falls within the band with (1-alpha)% probability,
#' rather than just pointwise coverage.
#'
#' @param lp_data Data frame containing rvar columns (from posterior package) to compute bands for.
#'   Should NOT include the x coordinate or any "true" value columns.
#' @param x_seq Numeric vector of x-coordinates (e.g., distances) corresponding to rows of lp_data.
#' @param alpha Significance level (default 0.05 for 95% simultaneous bands).
#' @param estimate_cols Character vector of column names to process. If NULL (default),
#'   processes all rvar columns in lp_data.
#'
#' @return A data frame with columns:
#'   \item{x}{Distance coordinates}
#'   \item{variable}{Name of the estimate column}
#'   \item{mean}{Posterior mean at each distance}
#'   \item{lower}{Simultaneous lower bound}
#'   \item{upper}{Simultaneous upper bound}
#'
#' @details
#' The algorithm:
#' 1. Extracts mean and SD at each distance for each rvar column
#' 2. Standardizes posterior deviations: (draw - mean) / SD
#' 3. Finds maximum absolute deviation across all distances for each draw
#' 4. Uses the (1-alpha) quantile of these maxima as the critical value
#' 5. Constructs bands: mean ± critical_value × SD
#'
#' @examples
#' # Single estimate
#' bands <- compute_simultaneous_bands(
#'   lp_data = tibble(estimate = lp_e),
#'   x_seq = seq(0, 100, 1),
#'   alpha = 0.05
#' )
#'
#' # Multiple estimates
#' bands <- compute_simultaneous_bands(
#'   lp_data = tibble(estimate_hier = lp_h, estimate_flat = lp_f),
#'   x_seq = seq(0, 100, 1),
#'   alpha = 0.05
#' )
compute_simultaneous_bands <- function(lp_data, x_seq, alpha = 0.05, estimate_cols = NULL) {

  lp_data <- as.data.frame(lp_data)

  # Determine which columns to process
  if (is.null(estimate_cols)) {
    # Find all rvar columns
    estimate_cols <- names(lp_data)[sapply(lp_data, posterior::is_rvar)]
  }

  # Validate inputs
  if (length(estimate_cols) == 0) {
    stop("No rvar columns found in lp_data")
  }
  if (nrow(lp_data) != length(x_seq)) {
    stop("Number of rows in lp_data must match length of x_seq")
  }

  # Step 1: Extract mean and SD per distance for all estimates
  df_summary <- lp_data %>%
    as.data.frame(check.names = FALSE) %>%
    mutate(x = x_seq) %>%
    mutate(
      across(
        all_of(estimate_cols),
        list(
          mn = ~as.vector(posterior::E(.)),
          sd = ~as.vector(sd(.))
        )
      ),
      .keep = "unused"
    )

  # Step 2: Standardize residuals at each distance separately
  # For each column, standardize each rvar element (distance point) by its own mean and SD
  lp_std_list <- list()

  for (col in estimate_cols) {
    rvar_col <- lp_data[[col]]
    mn_col <- df_summary[[paste0(col, "_mn")]]
    sd_col <- df_summary[[paste0(col, "_sd")]]

    # Extract draws matrix: rows = draws, cols = distances
    draws_matrix <- draws_of(rvar_col)

    # Standardize each distance point: (draws - mean) / sd
    # Handle SD = 0 case
    sd_col_safe <- ifelse(sd_col < 1e-10, 1, sd_col)

    # Standardize each column (distance) separately
    std_draws <- draws_matrix
    for (i in 1:ncol(draws_matrix)) {
      std_draws[, i] <- (draws_matrix[, i] - mn_col[i]) / sd_col_safe[i]
    }

    lp_std_list[[col]] <- rvar(std_draws)
  }

  lp_std <- as.data.frame(lp_std_list)
  
  # Step 3: Find maximum deviation across distances for each posterior draw
  max_dev <- sapply(lp_std, function(col) posterior::rvar_max(abs(col)))

  # Step 4: Find critical value for simultaneous coverage
  z_score_band <- sapply(max_dev, quantile, probs = 1 - alpha, na.rm = TRUE)

  # Handle cases where z_score is NA, NaN, or Inf
  z_score_band[is.na(z_score_band) | is.nan(z_score_band) | is.infinite(z_score_band)] <- 0

  # Step 5: Construct simultaneous lower and upper bands
  result <- df_summary %>%
    mutate(x = x_seq)

  # Add lower and upper bounds for each estimate
  for (col in estimate_cols) {
    mn_col <- paste0(col, "_mn")
    sd_col <- paste0(col, "_sd")

    result <- result %>%
      mutate(
        !!paste0(col, "_lower") := !!sym(mn_col) - z_score_band[[col]] * !!sym(sd_col),
        !!paste0(col, "_upper") := !!sym(mn_col) + z_score_band[[col]] * !!sym(sd_col)
      )
  }

  # Reshape to long format for easier plotting
  result_long <- result %>%
    select(x, ends_with("_mn"), ends_with("_lower"), ends_with("_upper")) %>%
    pivot_longer(
      cols = -x,
      names_to = c("variable", ".value"),
      names_pattern = "(.+)_(mn|lower|upper)$"
    ) %>%
    rename(mean = mn)

  return(result_long)
}

#' Compute Pointwise Credible Bands for Spatial Interaction Curves
#'
#' Computes pointwise credible bands using simple quantiles at each distance.
#' Unlike simultaneous bands, these only guarantee pointwise coverage - each
#' individual point has (1-alpha) probability of being in the band, but the
#' entire curve may not be contained.
#'
#' @param lp_data Data frame containing rvar columns (from posterior package) to compute bands for.
#'   Should NOT include the x coordinate or any "true" value columns.
#' @param x_seq Numeric vector of x-coordinates (e.g., distances) corresponding to rows of lp_data.
#' @param alpha Significance level (default 0.05 for 95% pointwise bands).
#' @param estimate_cols Character vector of column names to process. If NULL (default),
#'   processes all rvar columns in lp_data.
#'
#' @return A data frame with columns:
#'   \item{x}{Distance coordinates}
#'   \item{variable}{Name of the estimate column}
#'   \item{mean}{Posterior mean at each distance}
#'   \item{lower}{Pointwise lower bound}
#'   \item{upper}{Pointwise upper bound}
#'
#' @details
#' The algorithm:
#' 1. Extracts mean at each distance for each rvar column
#' 2. Computes (alpha/2) and (1-alpha/2) quantiles at each distance independently
#' 3. Returns mean ± quantile-based bounds
#'
#' @examples
#' # Single estimate
#' bands <- compute_pointwise_bands(
#'   lp_data = tibble(estimate = lp_e),
#'   x_seq = seq(0, 100, 1),
#'   alpha = 0.05
#' )
compute_pointwise_bands <- function(lp_data, x_seq, alpha = 0.05, estimate_cols = NULL) {

  lp_data <- as.data.frame(lp_data)

  # Determine which columns to process
  if (is.null(estimate_cols)) {
    # Find all rvar columns
    estimate_cols <- names(lp_data)[sapply(lp_data, posterior::is_rvar)]
  }

  # Validate inputs
  if (length(estimate_cols) == 0) {
    stop("No rvar columns found in lp_data")
  }
  if (nrow(lp_data) != length(x_seq)) {
    stop("Number of rows in lp_data must match length of x_seq")
  }

  # Compute mean and quantiles at each distance
  result_list <- list()

  for (col in estimate_cols) {
    rvar_col <- lp_data[[col]]

    result_list[[col]] <- tibble(
      x = x_seq,
      variable = col,
      mean = as.vector(posterior::E(rvar_col)),
      lower = sapply(rvar_col, quantile, probs = alpha / 2),
      upper = sapply(rvar_col, quantile, probs = 1 - alpha / 2)
    )
  }

  # Combine all estimates
  result_long <- bind_rows(result_list)

  return(result_long)
}
