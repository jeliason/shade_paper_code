# Shared functions for SHADE comparison simulations
# Used by both simulation_runner.R and null_calibration.R

#' Create patient metadata with proper indexing
create_patient_structure <- function(num_pts_per_group, num_images, seed = 2024) {
  set.seed(seed)

  # Generate images per patient
  total_patients <- num_pts_per_group * 2
  images_per_patient <- rep(num_images, total_patients)

  # Create patient metadata
  patients_df <- tibble(
    patient_id = 1:total_patients,
    group = c(rep("responder", num_pts_per_group), rep("nonresponder", num_pts_per_group)),
    group_idx = c(1:num_pts_per_group, 1:num_pts_per_group),
    images_count = images_per_patient
  )

  # Create image metadata
  total_images <- sum(images_per_patient)
  images_df <- tibble(
    image_id = 1:total_images,
    patient_id = rep(patients_df$patient_id, patients_df$images_count),
    group = rep(patients_df$group, patients_df$images_count),
    image_within_patient = sequence(patients_df$images_count)
  )

  return(list(
    patients_df = patients_df,
    images_df = images_df,
    total_patients = total_patients,
    total_images = total_images
  ))
}

#' Generate hierarchical spatial coefficients
generate_spatial_coefficients <- function(structure, num_potentials, seed = 2024) {
  set.seed(seed)

  num_responders <- sum(structure$patients_df$group == "responder")
  num_nonresponders <- sum(structure$patients_df$group == "nonresponder")

  # Group-level means (fixed effects)
  group_means <- matrix(0, nrow = 2, ncol = num_potentials)
  rownames(group_means) <- c("responder", "nonresponder")
  colnames(group_means) <- paste0("potential_", 1:num_potentials)

  # Responders: strong clustering (positive coefficients, decreasing with distance)
  group_means["responder", ] <- c(1.5, 1.0, 0.5)  # Short, medium, long range

  # Non-responders: weak repulsion (negative coefficients, decreasing with distance)
  group_means["nonresponder", ] <- c(-1.5, -1.0, -0.5)  # Short, medium, long range

  # Individual-level effects (random effects around group means)
  sigma_individual <- 0.1

  individual_effects <- array(0, dim = c(structure$total_patients, num_potentials))
  rownames(individual_effects) <- paste0("patient_", 1:structure$total_patients)
  colnames(individual_effects) <- paste0("potential_", 1:num_potentials)

  for(i in 1:structure$total_patients) {
    patient_group <- structure$patients_df$group[i]
    group_mean_vec <- group_means[patient_group, ]
    individual_effects[i, ] <- rnorm(num_potentials, mean = group_mean_vec, sd = sigma_individual)
  }

  # Image-level effects (random effects around individual means)
  sigma_image <- 0.1

  image_effects <- array(0, dim = c(structure$total_images, num_potentials))
  rownames(image_effects) <- paste0("image_", 1:structure$total_images)
  colnames(image_effects) <- paste0("potential_", 1:num_potentials)

  for(i in 1:structure$total_images) {
    patient_id <- structure$images_df$patient_id[i]
    individual_mean_vec <- individual_effects[patient_id, ]
    image_effects[i, ] <- rnorm(num_potentials, mean = individual_mean_vec, sd = sigma_image)
  }

  return(list(
    group_means = group_means,
    individual_effects = individual_effects,
    image_effects = image_effects
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
      patient_id = structure$images_df$patient_id[i],
      group = structure$images_df$group[i]
    )
  }

  return(patterns)
}

#' Run SHADE analysis
run_shade_analysis <- function(patterns, structure, num_potentials) {
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
    Group = structure$images_df$group
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

  stan_data$scale_sigma_betas <- c(10, 10, 10)
  stan_data$scale_sigmas <- 10

  cat("Fitting SHADE model...\n")
  shade_fit <- run_SHADE_model(stan_data,
                               method = "variational",
                               draws = 1e3,
                               refresh = 20,
                               threads = 2)

  return(list(fit = shade_fit, prep = prep))
}

#' Fit simple (flat) SHADE model to a single pattern
simple_shade <- function(pattern, potentials) {
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

    Q_grid <- make_quadrature(pattern, n_dummy = 50, dist = "grid")
    data_list_grid <- make_data(Q_grid, potentials, type, verbose = FALSE)

    X_new <- as.matrix(data_list_grid$data)[, -1]
    resp_new <- data_list_grid$response
    X_new <- X_new[resp_new == 0, ]
    N_new <- nrow(X_new)
    offset_new <- rep(unique(offset), N_new)

    data_stan <- list(
      N = nrow(x_cells),
      K = ncol(x_cells),
      y = is_cell,
      X = x_cells,
      oset = -offset,
      N_new = N_new,
      X_new = X_new,
      offset_new = offset_new
    )

    mod <- cmdstan_model("sim_shade_comparison/simple_shade.stan")
    fit <- mod$variational(data_stan, draws = 1e3, refresh = 100)

    list(fit = fit, data_stan = data_stan)
  }, error = function(e) {
    print(e)
    return(NULL)
  })
}

#' Run simple (flat) SHADE analysis on each image independently
run_simple_shade_analysis <- function(patterns, potentials) {
  cat("Running simple SHADE analysis...\n")

  # Fit simple SHADE model to each image independently
  image_fits <- lapply(seq_along(patterns), function(i) {
    if(i %% 10 == 0) {
      cat("Fitting image", i, "of", length(patterns), "\n")
    }

    tryCatch({
      pattern <- patterns[[i]]$pattern
      # Fit simple SHADE model - pass potentials explicitly
      results <- simple_shade(pattern, potentials)
      return(results)
    }, error = function(e) {
      cat("Error in simple SHADE for image", i, ":", conditionMessage(e), "\n")
      return(NULL)
    })
  })

  return(image_fits)
}

#' Run G-cross analysis
run_gcross_analysis <- function(patterns, structure) {
  cat("Running G-cross analysis...\n")

  # Image-level G-cross results
  image_results <- lapply(seq_along(patterns), function(i) {
    if(i %% 10 == 0) {
      print(i)
    }

    out <- tryCatch({

      pat <- patterns[[i]]$pattern
      true_group <- patterns[[i]]$group

      # Check cell counts
      n_t <- sum(pat$marks == "1")
      n_tumor <- sum(pat$marks == "3")

      if(n_t < 3 || n_tumor < 3) {
        stop("not enough cells!")
      }

      # Calculate G-cross with global (simultaneous) envelope
      # For global envelope at α=0.05, use nsim=199 minimum
      # Restrict to relevant distance range (10-75 μm) for proper calibration
      r_seq <- seq(10, 75, by = 1)
      env <- envelope(pat, Gcross, i = "1", j = "3", nsim = 199,
                      correction = "km", verbose = FALSE,
                      global = TRUE, r = r_seq)

      # Check for clustering or repulsion
      # Global envelope: check if observed function exits envelope ANYWHERE in range
      clustering_detected <- any(env$obs > env$hi, na.rm = TRUE)
      repulsion_detected  <- any(env$obs < env$lo, na.rm = TRUE)

      # Return detection result
      if(true_group == "responder") {
        detect <- as.numeric(clustering_detected)  # Should detect clustering
      } else {
        detect <- as.numeric(repulsion_detected)   # Should detect repulsion
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

#' Run K-cross analysis
run_kcross_analysis <- function(patterns, structure) {
  cat("Running K-cross analysis...\n")

  # Image-level K-cross results
  image_results <- lapply(seq_along(patterns), function(i) {
    if(i %% 10 == 0) {
      print(i)
    }

    out <- tryCatch({

      pat <- patterns[[i]]$pattern
      true_group <- patterns[[i]]$group

      # Check cell counts
      n_t <- sum(pat$marks == "1")
      n_tumor <- sum(pat$marks == "3")

      if(n_t < 3 || n_tumor < 3) {
        stop("not enough cells!")
      }

      # Calculate K-cross with global (simultaneous) envelope
      # For global envelope at α=0.05, use nsim=199 minimum
      # Restrict to relevant distance range (10-75 μm) for proper calibration
      r_seq <- seq(10, 75, by = 1)
      env <- envelope(pat, Kcross, i = "1", j = "3", nsim = 199,
                      correction = "iso", verbose = FALSE,
                      global = TRUE, r = r_seq)

      # Check for clustering or repulsion
      # Global envelope: check if observed function exits envelope ANYWHERE in range
      clustering_detected <- any(env$obs > env$hi, na.rm = TRUE)
      repulsion_detected  <- any(env$obs < env$lo, na.rm = TRUE)

      # Return detection result
      if(true_group == "responder") {
        detect <- as.numeric(clustering_detected)  # Should detect clustering
      } else {
        detect <- as.numeric(repulsion_detected)   # Should detect repulsion
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
