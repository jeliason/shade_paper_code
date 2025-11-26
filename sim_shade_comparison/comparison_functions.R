# Shared functions for SHADE comparison simulations
# Used by both 01_simulation_runner.R and 02_null_calibration.R

# ============================================================================
# CONSTANTS
# ============================================================================

# Distance range for all analyses (0-75 μm)
# Covers all three RBF basis functions
# Note: envelope() in spatstat requires r to start at 0
DISTANCE_RANGE_MIN <- 0
DISTANCE_RANGE_MAX <- 75

# ============================================================================
# FUNCTIONS
# ============================================================================

#' Plot SIC with simultaneous credible bands
#'
#' @param df_summ Data frame with columns: x, mn, lo_simul, hi_simul, lp_true, lp_true_patient, lp_true_cohort
#' @param image_id Image identifier for plot title
#' @param method Method name for plot title (e.g., "Simple SHADE", "SHADE Hierarchical")
#' @param verbose If TRUE, create and display plot
#' @param plot_detections If TRUE, plot only detections; if FALSE, plot only non-detections; if NULL, plot all
#' @param is_single_image For hierarchical models: if TRUE, main estimate is patient-level; if FALSE, it's image-level
plot_sic_detection <- function(df_summ, image_id, method = "SHADE", verbose = FALSE, plot_detections = TRUE, is_single_image = FALSE) {
  if(!verbose) return(invisible(NULL))

  # Check if detection occurred
  excludes_zero <- any(df_summ$lo_simul > 0 | df_summ$hi_simul < 0, na.rm = TRUE)

  # Filter based on plot_detections parameter
  if(!is.null(plot_detections)) {
    if(plot_detections && !excludes_zero) return(invisible(NULL))
    if(!plot_detections && excludes_zero) return(invisible(NULL))
  }

  # Log detection or non-detection
  if(excludes_zero) {
    detection_indices <- which(df_summ$lo_simul > 0 | df_summ$hi_simul < 0)
    detection_r <- df_summ$x[detection_indices]

    cat(sprintf("  *** Detection in %s image %d ***\n", method, image_id))
    cat(sprintf("      r values with detection: %s\n",
                paste(round(detection_r, 2), collapse = ", ")))
  } else {
    cat(sprintf("  --- No detection in %s image %d ---\n", method, image_id))
  }

  # Create ggplot
  status <- ifelse(excludes_zero, "Detection", "No Detection")
  p <- ggplot(df_summ, aes(x = x)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray50", size = 0.5) +
    geom_ribbon(aes(ymin = lo_simul, ymax = hi_simul, fill = "95% Credible Band"),
                alpha = 0.3)

  # Determine label for main estimate based on single-image status
  # For single-image patients: mn and lp_true are patient-level
  # For multi-image patients: mn and lp_true are image-level
  main_level <- if(is_single_image) "Patient" else "Image"

  # Add estimated SICs (solid lines)
  if("mn" %in% names(df_summ)) {
    p <- p + geom_line(aes(y = mn, color = !!main_level, linetype = "Estimated"), size = 1)
  }
  if("lp_patient" %in% names(df_summ) && !is_single_image) {
    # Only show patient line separately if we're looking at image-level
    p <- p + geom_line(aes(y = lp_patient, color = "Patient", linetype = "Estimated"), size = 0.8)
  }
  if("lp_global" %in% names(df_summ)) {
    p <- p + geom_line(aes(y = lp_global, color = "Global", linetype = "Estimated"), size = 0.8)
  }

  # Add true SICs (dashed lines)
  if("lp_true" %in% names(df_summ)) {
    p <- p + geom_line(aes(y = lp_true, color = !!main_level, linetype = "True"), size = 0.8)
  }
  if("lp_true_patient" %in% names(df_summ) && !is_single_image) {
    # Only show patient line separately if we're looking at image-level
    p <- p + geom_line(aes(y = lp_true_patient, color = "Patient", linetype = "True"), size = 0.8)
  }
  if("lp_true_cohort" %in% names(df_summ)) {
    p <- p + geom_line(aes(y = lp_true_cohort, color = "Global", linetype = "True"), size = 0.8)
  }

  # Define colors (by hierarchical level) and linetypes (estimated vs true)
  p <- p +
    scale_color_manual(
      name = "Level",
      values = c(
        "Image" = "blue",
        "Patient" = "orange",
        "Global" = "darkgreen"
      ),
      breaks = c("Image", "Patient", "Global")
    ) +
    scale_linetype_manual(
      name = "Type",
      values = c(
        "Estimated" = "solid",
        "True" = "dashed"
      ),
      breaks = c("Estimated", "True")
    ) +
    scale_fill_manual(
      name = "",
      values = c("95% Credible Band" = "lightblue")
    ) +
    labs(
      title = sprintf("%s: Image %d (%s)%s", method, image_id, status,
                      ifelse(is_single_image, " [Patient-level]", "")),
      x = "Distance (μm)",
      y = "Spatial Interaction Curve"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  print(p)

  return(invisible(NULL))
}

#' Create patient metadata with proper indexing
create_patient_structure <- function(num_patients, num_images, seed = 2024) {
  set.seed(seed)

  # Generate images per patient
  images_per_patient <- rep(num_images, num_patients)

  # Create patient metadata (single group)
  patients_df <- tibble(
    patient_id = 1:num_patients,
    images_count = images_per_patient
  )

  # Create image metadata
  total_images <- sum(images_per_patient)
  images_df <- tibble(
    image_id = 1:total_images,
    patient_id = rep(patients_df$patient_id, patients_df$images_count),
    image_within_patient = sequence(patients_df$images_count)
  )

  return(list(
    patients_df = patients_df,
    images_df = images_df,
    total_patients = num_patients,
    total_images = total_images
  ))
}

#' Generate hierarchical spatial coefficients
generate_spatial_coefficients <- function(structure, num_potentials, cohort_mean = c(1.5, 1.0, 0.5),
                                          sigma_individual = 0.4, sigma_image = 0.25, seed = 2024) {
  set.seed(seed)

  # Cohort-level mean (fixed effect) - single group with clustering
  cohort_mean_vec <- cohort_mean

  # Individual-level effects (random effects around cohort mean)
  individual_effects <- array(0, dim = c(structure$total_patients, num_potentials))
  rownames(individual_effects) <- paste0("patient_", 1:structure$total_patients)
  colnames(individual_effects) <- paste0("potential_", 1:num_potentials)

  for(i in 1:structure$total_patients) {
    individual_effects[i, ] <- rnorm(num_potentials, mean = cohort_mean_vec, sd = sigma_individual)
  }

  # Image-level effects (random effects around individual means)
  image_effects <- array(0, dim = c(structure$total_images, num_potentials))
  rownames(image_effects) <- paste0("image_", 1:structure$total_images)
  colnames(image_effects) <- paste0("potential_", 1:num_potentials)

  for(i in 1:structure$total_images) {
    patient_id <- structure$images_df$patient_id[i]
    individual_mean_vec <- individual_effects[patient_id, ]
    image_effects[i, ] <- rnorm(num_potentials, mean = individual_mean_vec, sd = sigma_image)
  }

  return(list(
    cohort_mean = cohort_mean_vec,
    individual_effects = individual_effects,
    image_effects = image_effects,
    sigma_individual = sigma_individual,
    sigma_image = sigma_image
  ))
}

#' Generate spatial patterns for all images
generate_spatial_patterns <- function(structure, coefficients, density_params, size_im = 1500, seed = 2024) {
  set.seed(seed)

  W <- owin(c(0, size_im), c(0, size_im))
  area <- size_im^2

  # Create RBF potentials
  potentials <- make_rbfs(n_basis_functions = 3, max_dist = 75, basis_function_sigma = 15)

  patterns <- vector("list", structure$total_images)

  for(i in 1:structure$total_images) {
    if(i %% 5 == 0) cat("Generating pattern", i, "of", structure$total_images, "\n")

    # Get image-specific parameters
    image_coeffs <- coefficients$image_effects[i, ]

    # Generate source cells (T cells and B cells)
    pat_sources <- rmpoispp(lambda = c(density_params$np_t/area, density_params$np_b/area), win = W)

    # Generate tumor cell density based on T cells only
    if(sum(pat_sources$marks == 1) > 0) {
      t_cells <- unmark(subset(pat_sources, marks == 1))

      # Create custom kernel using image-specific coefficients
      custom_kernel <- Vectorize(function(x, y) {
        d <- sqrt(x^2 + y^2)
        kernel_value <- sum(sapply(seq_along(image_coeffs), function(k) {
          image_coeffs[k] * potentials[[k]](d)
        }))
        return(kernel_value)
      })

      # Generate density surface
      dens <- smooth_density_fft(t_cells, custom_kernel, resolution = 128)

      # Generate tumor cells
      lambda_integral <- sum(exp(dens$v)) * (dens$xstep * dens$ystep)

      if(lambda_integral > 0) {
        beta0 <- log(density_params$np_tumor / lambda_integral)
        tumor_cells <- rpoispp(lambda = exp(dens + beta0))
      } else {
        tumor_cells <- rpoispp(lambda = density_params$np_tumor/area, win = W)
        beta0 <- log(density_params$np_tumor/area)
      }
    } else {
      # No T cells - uniform tumor distribution
      tumor_cells <- rpoispp(lambda = density_params$np_tumor/area, win = W)
      beta0 <- log(density_params$np_tumor/area)
    }

    # Set tumor cell marks
    marks(tumor_cells) <- factor(3)
    # Combine all cell types
    combined_pattern <- superimpose(pat_sources, tumor_cells)

    patterns[[i]] <- list(
      pattern = combined_pattern,
      beta0 = beta0,
      image_coeffs = image_coeffs,
      patient_id = structure$images_df$patient_id[i]
    )
  }

  return(patterns)
}

#' Scenario 2: Generate spatial patterns with spatial gradient confounder
#'
#' Target cells driven by source cells + unmeasured spatial gradient
#' @param gradient_type Type of gradient: "linear", "radial", "sinusoidal"
#' @param gradient_strength Strength of gradient effect (log-scale)
generate_spatial_patterns_gradient <- function(structure, coefficients, density_params,
                                                size_im = 1500, seed = 2024,
                                                gradient_type = "radial",
                                                gradient_strength = 1.0) {
  set.seed(seed)

  W <- owin(c(0, size_im), c(0, size_im))
  area <- size_im^2

  # Create RBF potentials
  potentials <- make_rbfs(n_basis_functions = 3, max_dist = 75, basis_function_sigma = 15)

  patterns <- vector("list", structure$total_images)

  for(i in 1:structure$total_images) {
    if(i %% 5 == 0) cat("Generating gradient pattern", i, "of", structure$total_images, "\n")

    image_coeffs <- coefficients$image_effects[i, ]

    # Generate source cells (Poisson)
    pat_sources <- rmpoispp(lambda = c(density_params$np_t/area, density_params$np_b/area), win = W)

    # Create gradient function
    gradient_fn <- function(x, y) {
      if(gradient_type == "linear") {
        # Linear gradient from left to right
        return(gradient_strength * (x - size_im/2) / (size_im/2))
      } else if(gradient_type == "radial") {
        # Radial gradient from center
        cx <- size_im / 2
        cy <- size_im / 2
        dist_from_center <- sqrt((x - cx)^2 + (y - cy)^2)
        max_dist <- sqrt(2) * size_im / 2
        return(gradient_strength * (dist_from_center / max_dist - 0.5))
      } else if(gradient_type == "sinusoidal") {
        # Sinusoidal wave pattern
        return(gradient_strength * sin(2 * pi * x / size_im))
      }
    }

    # Generate tumor cell density based on T cells + gradient
    if(sum(pat_sources$marks == 1) > 0) {
      t_cells <- unmark(subset(pat_sources, marks == 1))

      # Create custom kernel for source effect
      custom_kernel <- Vectorize(function(x, y) {
        d <- sqrt(x^2 + y^2)
        kernel_value <- sum(sapply(seq_along(image_coeffs), function(k) {
          image_coeffs[k] * potentials[[k]](d)
        }))
        return(kernel_value)
      })

      # Source-driven density
      dens <- smooth_density_fft(t_cells, custom_kernel, resolution = 128)

      # Add gradient effect
      gradient_surface <- as.im(function(x, y) gradient_fn(x, y), W = W, eps = dens$xstep)

      # Combine: log(lambda) = source effect + gradient effect
      combined_log_lambda <- dens + gradient_surface

      # Normalize to get desired number of cells
      lambda_integral <- sum(exp(combined_log_lambda$v)) * (dens$xstep * dens$ystep)

      if(lambda_integral > 0) {
        beta0 <- log(density_params$np_tumor / lambda_integral)
        tumor_cells <- rpoispp(lambda = exp(combined_log_lambda + beta0))
      } else {
        tumor_cells <- rpoispp(lambda = density_params$np_tumor/area, win = W)
        beta0 <- log(density_params$np_tumor/area)
      }
    } else {
      # No T cells - gradient only
      gradient_surface <- as.im(function(x, y) gradient_fn(x, y), W = W)
      lambda_integral <- sum(exp(gradient_surface$v)) * prod(gradient_surface$xstep, gradient_surface$ystep)
      beta0 <- log(density_params$np_tumor / lambda_integral)
      tumor_cells <- rpoispp(lambda = exp(gradient_surface + beta0))
    }

    marks(tumor_cells) <- factor(3)
    combined_pattern <- superimpose(pat_sources, tumor_cells)

    patterns[[i]] <- list(
      pattern = combined_pattern,
      beta0 = beta0,
      image_coeffs = image_coeffs,
      patient_id = structure$images_df$patient_id[i],
      gradient_type = gradient_type,
      gradient_strength = gradient_strength
    )
  }

  return(patterns)
}

#' Scenario 3: Generate spatial patterns with tissue compartment confounder
#'
#' Discrete spatial regions with different baseline target densities
#' @param n_compartments Number of spatial compartments (2-4)
#' @param compartment_effect Strength of compartment effect on baseline density
generate_spatial_patterns_compartments <- function(structure, coefficients, density_params,
                                                     size_im = 1500, seed = 2024,
                                                     n_compartments = 3,
                                                     compartment_effect = 0.5) {
  set.seed(seed)

  W <- owin(c(0, size_im), c(0, size_im))
  area <- size_im^2

  # Create RBF potentials
  potentials <- make_rbfs(n_basis_functions = 3, max_dist = 75, basis_function_sigma = 15)

  patterns <- vector("list", structure$total_images)

  for(i in 1:structure$total_images) {
    if(i %% 5 == 0) cat("Generating compartment pattern", i, "of", structure$total_images, "\n")

    image_coeffs <- coefficients$image_effects[i, ]

    # Create random compartments using Voronoi tessellation
    n_seeds <- n_compartments
    seeds_x <- runif(n_seeds, 0, size_im)
    seeds_y <- runif(n_seeds, 0, size_im)

    # Assign compartment labels (vectorized for grid evaluation)
    assign_compartment <- function(x, y) {
      # x and y are vectors of coordinates (grid points)
      # Compute distance from each point to each seed
      n_points <- length(x)

      # Create distance matrix: each row is a point, each column is a seed
      x_diff <- outer(x, seeds_x, "-")  # n_points x n_seeds
      y_diff <- outer(y, seeds_y, "-")  # n_points x n_seeds
      dist_mat <- sqrt(x_diff^2 + y_diff^2)

      # Find which seed is closest for each point
      apply(dist_mat, 1, which.min)
    }

    # Random compartment effects (multiplicative on baseline density)
    compartment_multipliers <- exp(rnorm(n_compartments, 0, compartment_effect))

    # Generate source cells (Poisson, independent of compartments)
    pat_sources <- rmpoispp(lambda = c(density_params$np_t/area, density_params$np_b/area), win = W)

    # Generate tumor cell density based on T cells + compartment
    if(sum(pat_sources$marks == 1) > 0) {
      t_cells <- unmark(subset(pat_sources, marks == 1))

      # Create custom kernel for source effect
      custom_kernel <- Vectorize(function(x, y) {
        d <- sqrt(x^2 + y^2)
        kernel_value <- sum(sapply(seq_along(image_coeffs), function(k) {
          image_coeffs[k] * potentials[[k]](d)
        }))
        return(kernel_value)
      })

      # Source-driven density
      dens <- smooth_density_fft(t_cells, custom_kernel, resolution = 128)

      # Create compartment image explicitly
      x_grid <- dens$xcol
      y_grid <- dens$yrow
      xy_expand <- expand.grid(y = y_grid, x = x_grid)  # Note: y first for matrix ordering
      comp_ids <- assign_compartment(xy_expand$x, xy_expand$y)
      comp_mat <- matrix(comp_ids, nrow = length(y_grid), ncol = length(x_grid))
      compartment_image <- im(comp_mat, xcol = x_grid, yrow = y_grid)

      # Create compartment effect surface
      effect_mat <- matrix(log(compartment_multipliers[comp_mat]),
                          nrow = length(y_grid), ncol = length(x_grid))
      compartment_surface <- im(effect_mat, xcol = x_grid, yrow = y_grid)

      # Combine: log(lambda) = source effect + compartment effect
      combined_log_lambda <- dens + compartment_surface

      lambda_integral <- sum(exp(combined_log_lambda$v)) * (dens$xstep * dens$ystep)

      if(lambda_integral > 0) {
        beta0 <- log(density_params$np_tumor / lambda_integral)
        tumor_cells <- rpoispp(lambda = exp(combined_log_lambda + beta0))
      } else {
        tumor_cells <- rpoispp(lambda = density_params$np_tumor/area, win = W)
        beta0 <- log(density_params$np_tumor/area)
      }
    } else {
      # No T cells - compartment effect only
      # Create grid
      resolution <- 128
      x_grid <- seq(W$xrange[1], W$xrange[2], length.out = resolution)
      y_grid <- seq(W$yrange[1], W$yrange[2], length.out = resolution)
      xy_expand <- expand.grid(y = y_grid, x = x_grid)

      # Assign compartments
      comp_ids <- assign_compartment(xy_expand$x, xy_expand$y)
      comp_mat <- matrix(comp_ids, nrow = length(y_grid), ncol = length(x_grid))
      compartment_image <- im(comp_mat, xcol = x_grid, yrow = y_grid)

      # Create compartment effect surface
      effect_mat <- matrix(log(compartment_multipliers[comp_mat]),
                          nrow = length(y_grid), ncol = length(x_grid))
      compartment_surface <- im(effect_mat, xcol = x_grid, yrow = y_grid)

      lambda_integral <- sum(exp(compartment_surface$v)) * (compartment_surface$xstep * compartment_surface$ystep)
      beta0 <- log(density_params$np_tumor / lambda_integral)
      tumor_cells <- rpoispp(lambda = exp(compartment_surface + beta0))
    }

    marks(tumor_cells) <- factor(3)
    combined_pattern <- superimpose(pat_sources, tumor_cells)

    patterns[[i]] <- list(
      pattern = combined_pattern,
      beta0 = beta0,
      image_coeffs = image_coeffs,
      patient_id = structure$images_df$patient_id[i],
      n_compartments = n_compartments,
      compartment_effect = compartment_effect,
      compartment_image = compartment_image,
      compartment_effects = compartment_multipliers
    )
  }

  return(patterns)
}

#' Scenario 4: Generate spatial patterns with unmeasured competing source
#'
#' Target cells driven by measured source + unmeasured competing source
#' @param unmeasured_strength Strength of unmeasured source effect relative to measured
generate_spatial_patterns_unmeasured_source <- function(structure, coefficients, density_params,
                                                          size_im = 1500, seed = 2024,
                                                          unmeasured_strength = 1.0) {
  set.seed(seed)

  W <- owin(c(0, size_im), c(0, size_im))
  area <- size_im^2

  # Create RBF potentials
  potentials <- make_rbfs(n_basis_functions = 3, max_dist = 75, basis_function_sigma = 15)

  patterns <- vector("list", structure$total_images)

  for(i in 1:structure$total_images) {
    if(i %% 5 == 0) cat("Generating unmeasured source pattern", i, "of", structure$total_images, "\n")

    image_coeffs <- coefficients$image_effects[i, ]

    # Generate measured source cells (T cells and B cells)
    pat_sources <- rmpoispp(lambda = c(density_params$np_t/area, density_params$np_b/area), win = W)

    # Generate UNMEASURED competing source cells (e.g., macrophages)
    # Use similar density as one of the measured sources
    unmeasured_density <- density_params$np_t  # Same as T cells
    unmeasured_source <- rpoispp(lambda = unmeasured_density/area, win = W)
    marks(unmeasured_source) <- factor(99)  # Different mark to distinguish

    # Generate tumor cell density based on BOTH measured T cells AND unmeasured source
    has_t_cells <- sum(pat_sources$marks == 1) > 0
    has_unmeasured <- npoints(unmeasured_source) > 0

    if(has_t_cells || has_unmeasured) {
      # Create custom kernel for measured source effect
      custom_kernel_measured <- Vectorize(function(x, y) {
        d <- sqrt(x^2 + y^2)
        kernel_value <- sum(sapply(seq_along(image_coeffs), function(k) {
          image_coeffs[k] * potentials[[k]](d)
        }))
        return(kernel_value)
      })

      # Create kernel for unmeasured source effect (same shape, scaled strength)
      custom_kernel_unmeasured <- Vectorize(function(x, y) {
        d <- sqrt(x^2 + y^2)
        kernel_value <- unmeasured_strength * sum(sapply(seq_along(image_coeffs), function(k) {
          image_coeffs[k] * potentials[[k]](d)
        }))
        return(kernel_value)
      })

      # Density from measured source
      if(has_t_cells) {
        t_cells <- unmark(subset(pat_sources, marks == 1))
        dens_measured <- smooth_density_fft(t_cells, custom_kernel_measured, resolution = 128)
      } else {
        dens_measured <- as.im(0, W = W)
      }

      # Density from unmeasured source
      if(has_unmeasured) {
        dens_unmeasured <- smooth_density_fft(unmeasured_source, custom_kernel_unmeasured, resolution = 128)
      } else {
        dens_unmeasured <- as.im(0, W = W)
      }

      # Combine: log(lambda) = measured effect + unmeasured effect
      combined_log_lambda <- dens_measured + dens_unmeasured

      lambda_integral <- sum(exp(combined_log_lambda$v)) * (dens_measured$xstep * dens_measured$ystep)

      if(lambda_integral > 0) {
        beta0 <- log(density_params$np_tumor / lambda_integral)
        tumor_cells <- rpoispp(lambda = exp(combined_log_lambda + beta0))
      } else {
        tumor_cells <- rpoispp(lambda = density_params$np_tumor/area, win = W)
        beta0 <- log(density_params$np_tumor/area)
      }
    } else {
      tumor_cells <- rpoispp(lambda = density_params$np_tumor/area, win = W)
      beta0 <- log(density_params$np_tumor/area)
    }

    marks(tumor_cells) <- factor(3)
    # Only include MEASURED sources in the pattern (unmeasured source is hidden)
    combined_pattern <- superimpose(pat_sources, tumor_cells)

    patterns[[i]] <- list(
      pattern = combined_pattern,
      beta0 = beta0,
      image_coeffs = image_coeffs,
      patient_id = structure$images_df$patient_id[i],
      unmeasured_strength = unmeasured_strength,
      unmeasured_pattern = unmeasured_source  # Store for diagnostics
    )
  }

  return(patterns)
}

#' Run SHADE analysis
run_shade_analysis <- function(patterns, structure, num_potentials,
                               method = "variational",
                               parameterization = "centered",
                               scale_sigma = 5) {
  cat("Preparing data for SHADE...\n")

  # Prepare coordinates data
  coords <- do.call(rbind, lapply(seq_along(patterns), function(i) {
    pat <- patterns[[i]]$pattern

    tibble(
      x = pat$x,
      y = pat$y,
      type = marks(pat),
      image_id = paste0("img_", i)
    )
  }))

  # Create patient metadata for SHADE
  image_ids <- unique(coords$image_id)
  coords$image_id <- factor(coords$image_id, level = image_ids)
  patient_metadata <- tibble(
    Spot = image_ids,
    Patient = paste0("pt_", structure$images_df$patient_id),
    Group = 1  # Single group
  )

  # Prepare model data
  prep <- prepare_spatial_model_data(
    x = coords$x,
    y = coords$y,
    cell_type = coords$type,
    image_id = coords$image_id,
    patient_metadata = patient_metadata,
    type_idx = 3,  # Tumor cells are type 3
    n_dummy = 1e3,
    n_basis_functions = num_potentials
  )
  stan_data <- prep$stan_data

  stan_data$scale_sigma_betas <- rep(scale_sigma, num_potentials)
  stan_data$scale_sigmas <- scale_sigma

  cat("Fitting SHADE model with method =", method, ", parameterization =", parameterization, ", scale =", scale_sigma, "\n")

  if(method == "sampling") {
    shade_fit <- run_SHADE_model(stan_data,
                                 method = "sample",  # cmdstanr uses "sample" not "sampling"
                                 parameterization = parameterization,
                                 chains = 1,
                                 iter_warmup = 500,
                                 iter_sampling = 500,
                                 refresh = 20,
                                 threads = 2)
  } else {
    shade_fit <- run_SHADE_model(stan_data,
                                 method = "variational",
                                 parameterization = parameterization,
                                 draws = 1e3,
                                 refresh = 20,
                                 threads = 2)
  }

  try(shade_fit$draws(),silent=TRUE)

  return(list(fit = shade_fit, prep = prep))
}

#' Fit simple (flat) SHADE model to a single pattern
simple_shade <- function(pattern, potentials, mod) {
  tryCatch({
    # SHADE model parameters
    n_dummy <- 1000
    type <- "3"

    Q <- make_quadrature(pattern, n_dummy = n_dummy)

    offset <- log(spatstat.geom::intensity(Q$dummy)) |>
      tibble::enframe() |>
      dplyr::right_join(tibble::tibble(name = spatstat.geom::marks(Q)), by = "name") |>
      dplyr::filter(name == type) |>
      dplyr::pull(value)

    data_list <- make_data(Q, potentials, type, verbose = FALSE)

    x_cells <- as.matrix(data_list$data)[, -1]
    is_cell <- data_list$response

    data_stan <- list(
      N = nrow(x_cells),
      K = ncol(x_cells),
      y = is_cell,
      X = x_cells,
      oset = -offset
    )

    fit <- mod$variational(data_stan, draws = 1e3, refresh = 0, show_messages=FALSE)

    try(fit$draws(),silent=TRUE)

    results <- list(fit = fit, data_stan = data_stan)
  }, error = function(e) {
    print(e)
    return(NULL)
  })
}

#' Run simple (flat) SHADE analysis on each image independently
run_simple_shade_analysis <- function(patterns, potentials, mod) {
  cat("Running simple SHADE analysis...\n")

  # Fit simple SHADE model to each image independently
  image_fits <- lapply(seq_along(patterns), function(i) {
    if(i %% 10 == 0) {
      cat("Fitting image", i, "of", length(patterns), "\n")
    }

    tryCatch({
      pattern <- patterns[[i]]$pattern
      # Fit simple SHADE model - pass potentials explicitly
      results <- simple_shade(pattern, potentials, mod)
      return(results)
    }, error = function(e) {
      cat("Error in simple SHADE for image", i, ":", conditionMessage(e), "\n")
      return(NULL)
    })
  })

  return(image_fits)
}

#' Run G-cross analysis using MAD test
run_gcross_analysis <- function(patterns, structure, verbose = FALSE) {
  cat("Running G-cross analysis...\n")

  # Track detections
  n_detected <- 0
  n_processed <- 0

  # Image-level G-cross results
  image_results <- lapply(seq_along(patterns), function(i) {
    out <- tryCatch({

      pat <- patterns[[i]]$pattern

      # Check cell counts
      n_t <- sum(pat$marks == "1")
      n_tumor <- sum(pat$marks == "3")

      if(n_t < 3 || n_tumor < 3) {
        stop("not enough cells!")
      }

      # Calculate envelope
      r_seq <- seq(DISTANCE_RANGE_MIN, DISTANCE_RANGE_MAX, by = 1)

      # Two-sided test: detect any deviation from CSR
      env <- envelope(pat, Gcross, i = "1", j = "3",
                     nsim = 99,
                     nrank = 10,
                     correction = "km",
                     global = TRUE,
                     r = r_seq,
                     verbose = FALSE,
                     savefuns = TRUE)

      # Detect any deviation from envelope (two-sided)
      detect <- as.numeric(any(env$obs > env$hi | env$obs < env$lo, na.rm = TRUE))
      detection_indices <- which(env$obs > env$hi | env$obs < env$lo)

      # If detection and verbose, log and plot
      if(detect == 1 && verbose) {
        detection_r <- env$r[detection_indices]
        cat(sprintf("  *** Detection in image %d ***\n", i))
        cat(sprintf("      r values with detection: %s\n",
                    paste(round(detection_r, 2), collapse = ", ")))
        plot(env, main = sprintf("G-cross: Image %d", i))
      }

      # Update running total
      n_processed <<- n_processed + 1
      n_detected <<- n_detected + detect

      # Log progress every 10 images
      if(i %% 10 == 0) {
        cat(sprintf("  Image %d/%d: %d detections (%.1f%%)\n",
                    i, length(patterns), n_detected, 100 * n_detected / n_processed))
      }

      out <- list(mean_detect = detect)
    }, error = function(e) {
      cat("Error in G-cross for image", i, ":", conditionMessage(e), "\n")
      out <- list(mean_detect = NA)
    })
    return(out)
  })

  return(image_results)
}

#' Run K-cross analysis using MAD test
run_kcross_analysis <- function(patterns, structure, verbose = FALSE) {
  cat("Running K-cross analysis...\n")

  # Track detections
  n_detected <- 0
  n_processed <- 0

  # Image-level K-cross results
  image_results <- lapply(seq_along(patterns), function(i) {
    out <- tryCatch({

      pat <- patterns[[i]]$pattern

      # Check cell counts
      n_t <- sum(pat$marks == "1")
      n_tumor <- sum(pat$marks == "3")

      if(n_t < 3 || n_tumor < 3) {
        stop("not enough cells!")
      }

      # Calculate envelope
      r_seq <- seq(DISTANCE_RANGE_MIN, DISTANCE_RANGE_MAX, by = 1)

      # Two-sided test: detect any deviation from CSR
      env <- envelope(pat, Lcross, i = "1", j = "3",
                     nsim = 99,
                     nrank = 8,
                     correction = "iso",
                     global = TRUE,
                     r = r_seq,
                     verbose = FALSE,
                     savefuns = TRUE)

      # Detect any deviation from envelope (two-sided)
      detect <- as.numeric(any(env$obs > env$hi | env$obs < env$lo, na.rm = TRUE))
      detection_indices <- which(env$obs > env$hi | env$obs < env$lo)

      # If detection and verbose, log and plot
      if(detect == 1 && verbose) {
        detection_r <- env$r[detection_indices]
        cat(sprintf("  *** Detection in image %d ***\n", i))
        cat(sprintf("      r values with detection: %s\n",
                    paste(round(detection_r, 2), collapse = ", ")))
        plot(env, main = sprintf("L-cross: Image %d", i))
      }

      # Update running total
      n_processed <<- n_processed + 1
      n_detected <<- n_detected + detect

      # Log progress every 10 images
      if(i %% 10 == 0) {
        cat(sprintf("  Image %d/%d: %d detections (%.1f%%)\n",
                    i, length(patterns), n_detected, 100 * n_detected / n_processed))
      }

      out <- list(mean_detect = detect)
    }, error = function(e) {
      cat("Error in K-cross for image", i, ":", conditionMessage(e), "\n")
      out <- list(mean_detect = NA)
    })
    return(out)
  })

  return(image_results)
}

#' Apply hard-core thinning to a point pattern
#'
#' Removes points to enforce minimum distance r_hc between any two points.
#' Uses sequential thinning: keep first point, then keep each subsequent point
#' only if it is >= r_hc from all previously kept points.
#'
#' @param pattern A marked point pattern (ppp object)
#' @param r_hc Hard-core radius in same units as pattern (default: 10 microns)
#' @return A thinned marked point pattern
apply_hardcore_thinning <- function(pattern, r_hc = 10) {
  # Convert to data frame
  df <- data.frame(
    x = pattern$x,
    y = pattern$y,
    mark = marks(pattern)
  )

  n <- nrow(df)
  if(n == 0) return(pattern)

  # Sequential thinning: keep point i if distance to all kept points >= r_hc
  keep <- rep(FALSE, n)
  keep[1] <- TRUE  # Always keep first point

  for(i in 2:n) {
    # Calculate distances to all previously kept points
    kept_indices <- which(keep[1:(i-1)])
    if(length(kept_indices) == 0) {
      keep[i] <- TRUE
      next
    }

    dists <- sqrt((df$x[i] - df$x[kept_indices])^2 +
                  (df$y[i] - df$y[kept_indices])^2)

    # Keep if minimum distance >= r_hc
    if(min(dists) >= r_hc) {
      keep[i] <- TRUE
    }
  }

  # Create thinned pattern
  if(sum(keep) == 0) {
    # Return empty pattern with same window
    return(ppp(numeric(0), numeric(0), window = pattern$window,
               marks = factor(levels = levels(marks(pattern)))))
  }

  df_thinned <- df[keep, ]

  pattern_thinned <- ppp(
    x = df_thinned$x,
    y = df_thinned$y,
    window = pattern$window,
    marks = factor(df_thinned$mark, levels = levels(marks(pattern)))
  )

  return(pattern_thinned)
}
