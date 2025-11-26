library(spatstat)
library(tidyverse)
library(Matrix)
library(SHADE)
library(cmdstanr)

# Load utility functions and constants
source("./utils.R")

# Load shared comparison functions

source("sim_shade_comparison/comparison_functions.R")
start <- Sys.time()

# ============================================================================
# MAIN SIMULATION SCRIPT
# ============================================================================

# Set environment and get data path
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  sim_idx <- 48
  verbose <- TRUE  # Set to TRUE to enable plotting

} else {
  args <- commandArgs(trailingOnly=TRUE)
  sim_idx <- as.numeric(args[1])
  verbose <- FALSE
  cat("Running simulation", sim_idx, "\n")
}

# Get data path (handles HPC vs local automatically)
path <- get_data_path("sim_shade_comparison")
print(paste("Data path:", path))

file_out <- paste0(path,"sim_",sim_idx,".rds")

# Check for force rerun flag
force_rerun <- Sys.getenv("FORCE_RERUN", unset = "0") == "1"

if (file.exists(file_out) & SYSTEM_ENV == "HPC" & !force_rerun) {
  cat("Simulation", sim_idx, "already completed. Skipping.\n")
  quit(save="no")
}

if (force_rerun & file.exists(file_out)) {
  cat("FORCE MODE: Rerunning simulation", sim_idx, "even though output exists.\n")
}

# EXPERIMENTAL GRID (testing MCMC/VI, centered/non-centered, prior scales)
grid <- expand.grid(
  method = c("variational"),# "sampling"),
  parameterization = c("centered"),
  scale_sigma = c(5),
  t_density = c("high", "low"),        # Test both to see if data quantity matters
  tumor_density = c("high","low"),              # Keep tumor density high
  num_images = c(1,2,3),                      # Best case: 3 images per patient
  sim_rep = 1:50
)

# Extract current condition parameters
condition <- list(
  method = as.character(grid$method[sim_idx]),
  parameterization = as.character(grid$parameterization[sim_idx]),
  scale_sigma = grid$scale_sigma[sim_idx],
  num_images = as.character(grid$num_images[sim_idx]),
  t_density = as.character(grid$t_density[sim_idx]),
  tumor_density = as.character(grid$tumor_density[sim_idx]),
  sim_rep = grid$sim_rep[sim_idx]
)

cat("Condition:", condition$method, condition$parameterization,
    "scale=", condition$scale_sigma,
    condition$t_density, condition$tumor_density, condition$num_images,
    "Rep:", condition$sim_rep, "\n")

# Fixed parameters
num_patients <- 40  # Single group, all with non-zero coefficients
num_potentials <- 3
size_im <- 1500

# Density parameters
density_params <- list(
  np_t = ifelse(condition$t_density == "high", 150, 15),
  np_tumor = ifelse(condition$tumor_density == "high", 150, 15),
  np_b = ifelse(condition$t_density == "high", 150, 15)  # B cells same as T cells
)

# ============================================================================
# RUN SIMULATION
# ============================================================================

cat("Creating patient structure...\n")
structure <- create_patient_structure(
  num_patients = num_patients,
  num_images = grid$num_images[sim_idx],
  seed = 2024 + sim_idx
)

cat("Generating spatial coefficients...\n")
coefficients <- generate_spatial_coefficients(
  structure = structure,
  num_potentials = num_potentials,
  cohort_mean = c(1.5, 1.0, 0.5),  # All patients have clustering
  seed = 2024 + condition$sim_rep
)


cat("Generating spatial patterns...\n")
patterns <- generate_spatial_patterns(
  structure = structure,
  coefficients = coefficients,
  density_params = density_params,
  size_im = size_im,
  seed = 2024 + sim_idx + 1000
)

cat("Running SHADE analysis...\n")
shade_results <- run_shade_analysis(
  patterns,
  structure,
  num_potentials,
  method = condition$method,
  parameterization = condition$parameterization,
  scale_sigma = condition$scale_sigma
)

# Define potentials for simple SHADE

potentials <- make_rbfs(max_dist = 75, n_basis_functions = 3, basis_function_sigma = 15)
mod <- cmdstan_model("sim_shade_comparison/super_simple_shade.stan", quiet=TRUE)
simple_shade_results <- run_simple_shade_analysis(patterns, potentials, mod)

cat("Running G-cross analysis...\n")

gcross_results <- run_gcross_analysis(patterns, structure)

cat("Running K-cross analysis...\n")
kcross_results <- run_kcross_analysis(patterns, structure)

# ============================================================================
# EXTRACT RESULTS AND CALCULATE METRICS
# ============================================================================

# Extract SHADE results
library(posterior)
shade_draws <- as_draws_rvars(shade_results$fit$draws())

# Extract simple SHADE results and calculate image-level power
simple_shade_image_results <- sapply(1:structure$total_images, function(i) {
  tryCatch({
    results <- simple_shade_results[[i]]
    if(is.null(results)) return(NULL)

    shade_draws <- as_draws_rvars(results$fit$draws())

    # Get image-level coefficients for all 3 potentials
    beta <- as.vector(shade_draws$beta[1:3])
    beta_true <- coefficients$image_effects[i,]

    # Get patient-level and cohort-level true coefficients
    patient_id <- structure$images_df$patient_id[i]
    beta_true_patient <- coefficients$individual_effects[patient_id,]
    beta_true_cohort <- coefficients$cohort_mean

    # Compute SIC using same distance range as global envelopes
    x_seq <- seq(DISTANCE_RANGE_MIN, DISTANCE_RANGE_MAX, 1)
    x_des <- lapply(potentials, \(pot) pot(x_seq)) %>% do.call(cbind,.)

    lp <- as.vector(x_des %*% beta)
    lp_true <- as.vector(x_des %*% beta_true)
    lp_true_patient <- as.vector(x_des %*% beta_true_patient)
    lp_true_cohort <- as.vector(x_des %*% beta_true_cohort)

    # Compute simultaneous 95% credible bands
    lp_data <- lp %>% as.matrix() %>% as.data.frame()
    colnames(lp_data) <- "SIC"

    df_summ <- compute_simultaneous_bands(
      lp_data = lp_data,
      x_seq = x_seq,
      alpha = 0.05
    ) %>%
      rename(mn = mean, lo_simul = lower, hi_simul = upper) %>%
      mutate(
        lo_pw = as.vector(quantile(lp, probs = 0.025)),
        hi_pw = as.vector(quantile(lp, probs = 0.975)),
        lp_true = lp_true,
        lp_true_patient = lp_true_patient,
        lp_true_cohort = lp_true_cohort
      )

    # Check if simultaneous credible band excludes zero anywhere in the range
    excludes_zero <- any(df_summ$lo_simul > 0 | df_summ$hi_simul < 0, na.rm = TRUE)

    # Check if simultaneous credible band contains true SIC everywhere (coverage)
    coverage <- all(df_summ$lo_simul <= df_summ$lp_true &
                    df_summ$hi_simul >= df_summ$lp_true, na.rm = TRUE)

    # Plot if detection and verbose mode
    plot_sic_detection(df_summ, image_id = i, method = "Simple SHADE", verbose = verbose)

    return(list(detection = as.numeric(excludes_zero),
                coverage = as.numeric(coverage)))

  }, error = function(e) {
    cat("Error extracting simple SHADE results for image", i, ":", conditionMessage(e), "\n")
    return(NULL)
  })
}, simplify = FALSE)

# Calculate SHADE image-level power
shade_image_results <- sapply(1:structure$total_images, function(i) {
  # Get patient info
  patient_id <- structure$images_df$patient_id[i]
  num_images_for_patient <- sum(structure$images_df$patient_id == patient_id)

  # For single-image patients: use patient-level estimates and truth
  # For multi-image patients: use image-level estimates and truth
  if(num_images_for_patient == 1) {
    beta <- as.vector(shade_draws$beta_indiv[2:4,patient_id])
    beta_true <- coefficients$individual_effects[patient_id,]
  } else {
    beta <- as.vector(shade_draws$beta_local[i,2:4])
    beta_true <- coefficients$image_effects[i,]
  }

  # Get patient-level and cohort-level true coefficients
  beta_true_patient <- coefficients$individual_effects[patient_id,]
  beta_true_cohort <- coefficients$cohort_mean

  # Get estimated patient-level and global coefficients
  beta_patient <- as.vector(shade_draws$beta_indiv[2:4,patient_id])
  beta_global <- as.vector(shade_draws$beta_global[2:4])

  # Compute SIC using same distance range as global envelopes
  x_seq <- seq(DISTANCE_RANGE_MIN, DISTANCE_RANGE_MAX, 1)
  x_des <- lapply(potentials,\(pot) pot(x_seq)) %>% do.call(cbind,.)

  # Compute estimated SICs at all levels
  lp <- as.vector(x_des %*% beta)
  lp_patient_rvar <- as.vector(x_des %*% beta_patient)
  lp_global_rvar <- as.vector(x_des %*% beta_global)

  # Extract mean values for plotting (rvars need to be converted to numeric)
  lp_patient <- as.vector(mean(lp_patient_rvar))
  lp_global <- as.vector(mean(lp_global_rvar))

  # Compute true SICs at all levels
  lp_true <- as.vector(x_des %*% beta_true)
  lp_true_patient <- as.vector(x_des %*% beta_true_patient)
  lp_true_cohort <- as.vector(x_des %*% beta_true_cohort)

  # Compute simultaneous 95% credible bands using utility function
  lp_data <- lp %>% as.matrix() %>% as.data.frame()
  colnames(lp_data) <- "SIC"

  df_summ <- compute_simultaneous_bands(
    lp_data = lp_data,
    x_seq = x_seq,
    alpha = 0.05
  ) %>%
    rename(mn = mean, lo_simul = lower, hi_simul = upper) %>%
    mutate(
      lo_pw = as.vector(quantile(lp, probs = 0.025)),
      hi_pw = as.vector(quantile(lp, probs = 0.975)),
      lp_patient = lp_patient,
      lp_global = lp_global,
      lp_true = lp_true,
      lp_true_patient = lp_true_patient,
      lp_true_cohort = lp_true_cohort
    )

  # Check if simultaneous credible band excludes zero anywhere in the range
  excludes_zero <- any(df_summ$lo_simul > 0 | df_summ$hi_simul < 0, na.rm = TRUE)

  # Check if simultaneous credible band contains true SIC everywhere (coverage)
  # NOTE: Hierarchical models shrink image-level estimates toward patient means,
  # so coverage of image-level CIs will be lower than nominal when there is
  # substantial image-level variation (sigma_image = 0.25). This is expected
  # behavior and not directly comparable to G-cross/K-cross (which don't report coverage).
  coverage <- all(df_summ$lo_simul <= df_summ$lp_true &
                  df_summ$hi_simul >= df_summ$lp_true, na.rm = TRUE)

  # Plot if detection and verbose mode
  plot_sic_detection(df_summ, image_id = i, method = "SHADE Hierarchical", verbose = verbose,
                     plot_detections = NULL, is_single_image = (num_images_for_patient == 1))

  return(list(detection = as.numeric(excludes_zero),
              coverage = as.numeric(coverage)))
}, simplify = FALSE)

# Calculate performance metrics
# Extract detection and coverage for SHADE methods
shade_detections <- sapply(shade_image_results, function(x) if(is.list(x)) x$detection else NA)
shade_coverages <- sapply(shade_image_results, function(x) if(is.list(x)) x$coverage else NA)

simple_shade_detections <- sapply(compact(simple_shade_image_results), function(x) if(is.list(x)) x$detection else NA)
simple_shade_coverages <- sapply(compact(simple_shade_image_results), function(x) if(is.list(x)) x$coverage else NA)

# Power = proportion of images with detection (simultaneous band excludes zero)
shade_image_power <- mean(shade_detections, na.rm = TRUE)
simple_shade_image_power <- mean(simple_shade_detections, na.rm = TRUE)

# Coverage = proportion of images where credible band contains true SIC everywhere
shade_image_coverage <- mean(shade_coverages, na.rm = TRUE)
simple_shade_image_coverage <- mean(simple_shade_coverages, na.rm = TRUE)

# G-cross and K-cross (detection only, no coverage metric)
gx_detect <- unlist(lapply(gcross_results,\(o) o$mean_detect))
gcross_image_power <- mean(gx_detect, na.rm = TRUE)

kx_detect <- unlist(lapply(kcross_results,\(o) o$mean_detect))
kcross_image_power <- mean(kx_detect, na.rm = TRUE)

# ============================================================================
# NULL CALIBRATION (Type I Error)
# ============================================================================

cat("\n=== NULL CALIBRATION ===\n")
cat("Generating NULL patterns (no spatial association)...\n")

# Generate NULL coefficients (all zeros with no hierarchical variation)
coefficients_null <- generate_spatial_coefficients(
  structure = structure,
  num_potentials = num_potentials,
  cohort_mean = c(0, 0, 0),  # NULL: no spatial interaction
  sigma_individual = 0,      # NULL: no patient-level variation
  sigma_image = 0,           # NULL: no image-level variation
  seed = 2024 + condition$sim_rep + 5000
)

# Generate NULL spatial patterns (CSR for tumor cells)
W <- owin(c(0, size_im), c(0, size_im))
area <- size_im^2

patterns_null <- lapply(1:structure$total_images, function(i) {
  # Generate source cells (T cells and B cells) - homogeneous Poisson
  pat_sources <- rmpoispp(lambda = c(density_params$np_t/area, density_params$np_b/area), win = W)

  # Generate tumor cells UNIFORMLY (NULL: no dependence on source cells)
  tumor_cells <- rpoispp(lambda = density_params$np_tumor/area, win = W)
  marks(tumor_cells) <- factor(3)

  # Combine all cell types
  combined_pattern <- superimpose(pat_sources, tumor_cells)

  list(
    pattern = combined_pattern,
    beta0 = log(density_params$np_tumor/area),
    image_coeffs = coefficients_null$image_effects[i, ],  # All zeros
    patient_id = structure$images_df$patient_id[i]
  )
})

cat("Running NULL analyses...\n")

# Run SHADE on null patterns
shade_results_null <- run_shade_analysis(
  patterns_null,
  structure,
  num_potentials,
  method = condition$method,
  parameterization = condition$parameterization,
  scale_sigma = condition$scale_sigma
)
shade_draws_null <- as_draws_rvars(shade_results_null$fit$draws())

# Run simple SHADE on null patterns
simple_shade_results_null <- run_simple_shade_analysis(patterns_null, potentials, mod)

# Run G-cross on null patterns
gcross_results_null <- run_gcross_analysis(patterns_null, structure, verbose = FALSE)

# Run K-cross on null patterns
kcross_results_null <- run_kcross_analysis(patterns_null, structure, verbose = FALSE)

cat("Calculating Type I error rates...\n")

# SHADE hierarchical: false positives
shade_fp <- sapply(1:structure$total_images, function(i) {
  # Get patient info
  patient_id <- structure$images_df$patient_id[i]
  num_images_for_patient <- sum(structure$images_df$patient_id == patient_id)

  # For single-image patients: use patient-level estimates
  # For multi-image patients: use image-level estimates
  if(num_images_for_patient == 1) {
    beta <- as.vector(shade_draws_null$beta_indiv[2:4,patient_id])
  } else {
    beta <- as.vector(shade_draws_null$beta_local[i,2:4])
  }

  # Get patient-level and global coefficients for plotting
  beta_patient <- as.vector(shade_draws_null$beta_indiv[2:4,patient_id])
  beta_global <- as.vector(shade_draws_null$beta_global[2:4])

  x_seq <- seq(DISTANCE_RANGE_MIN, DISTANCE_RANGE_MAX, 1)
  x_des <- lapply(potentials, \(pot) pot(x_seq)) %>% do.call(cbind,.)

  # Compute estimated SICs at all levels
  lp <- as.vector(x_des %*% beta)
  lp_patient_rvar <- as.vector(x_des %*% beta_patient)
  lp_global_rvar <- as.vector(x_des %*% beta_global)

  # Extract mean values for plotting
  lp_patient <- as.vector(mean(lp_patient_rvar))
  lp_global <- as.vector(mean(lp_global_rvar))

  lp_data <- lp %>% as.matrix() %>% as.data.frame()
  colnames(lp_data) <- "SIC"

  df_summ <- compute_simultaneous_bands(lp_data = lp_data, x_seq = x_seq, alpha = 0.05) %>%
    rename(mn = mean, lo_simul = lower, hi_simul = upper) %>%
    mutate(
      lo_pw = as.vector(quantile(lp, probs = 0.025)),
      hi_pw = as.vector(quantile(lp, probs = 0.975)),
      lp_patient = lp_patient,
      lp_global = lp_global,
      lp_true = 0,  # TRUE SIC is zero under null
      lp_true_patient = 0,
      lp_true_cohort = 0
    )

  excludes_zero <- any(df_summ$lo_simul > 0 | df_summ$hi_simul < 0, na.rm = TRUE)

  # Plot false positives if detected and verbose mode
  if(excludes_zero) {
    plot_sic_detection(df_summ, image_id = i, method = "SHADE Hierarchical (NULL)", verbose = verbose,
                       plot_detections = NULL, is_single_image = (num_images_for_patient == 1))
  }

  return(as.numeric(excludes_zero))
})

# Simple SHADE: false positives
simple_shade_fp <- sapply(1:structure$total_images, function(i) {
  results <- simple_shade_results_null[[i]]
  if(is.null(results)) return(NA)

  beta_draws <- as_draws_rvars(results$fit$draws(variables = "beta"))$beta
  beta <- beta_draws[1:3]

  x_seq <- seq(DISTANCE_RANGE_MIN, DISTANCE_RANGE_MAX, 1)
  x_des <- lapply(potentials, \(pot) pot(x_seq)) %>% do.call(cbind,.)
  lp <- as.vector(x_des %*% beta)

  lp_data <- lp %>% as.matrix() %>% as.data.frame()
  colnames(lp_data) <- "SIC"

  df_summ <- compute_simultaneous_bands(lp_data = lp_data, x_seq = x_seq, alpha = 0.05) %>%
    rename(mn = mean, lo_simul = lower, hi_simul = upper) %>%
    mutate(
      lo_pw = as.vector(quantile(lp, probs = 0.025)),
      hi_pw = as.vector(quantile(lp, probs = 0.975)),
      lp_true = 0,  # TRUE SIC is zero under null
      lp_true_patient = 0,
      lp_true_cohort = 0
    )

  excludes_zero <- any(df_summ$lo_simul > 0 | df_summ$hi_simul < 0, na.rm = TRUE)

  # Plot false positives if detected and verbose mode
  if(excludes_zero) {
    plot_sic_detection(df_summ, image_id = i, method = "Simple SHADE (NULL)", verbose = verbose)
  }

  return(as.numeric(excludes_zero))
})

# G-cross and K-cross: false positives
gcross_fp <- sapply(gcross_results_null, function(res) {
  if(is.na(res$mean_detect)) return(NA)
  return(res$mean_detect)
})

kcross_fp <- sapply(kcross_results_null, function(res) {
  if(is.na(res$mean_detect)) return(NA)
  return(res$mean_detect)
})

# Type I error rates
type_i_error_shade <- mean(shade_fp, na.rm = TRUE)
type_i_error_simple <- mean(simple_shade_fp, na.rm = TRUE)
type_i_error_gcross <- mean(gcross_fp, na.rm = TRUE)
type_i_error_kcross <- mean(kcross_fp, na.rm = TRUE)

cat("\nType I Error Rates:\n")
cat(sprintf("  SHADE Hierarchical: %.1f%%\n", type_i_error_shade * 100))
cat(sprintf("  SHADE Flat:         %.1f%%\n", type_i_error_simple * 100))
cat(sprintf("  G-cross:            %.1f%%\n", type_i_error_gcross * 100))
cat(sprintf("  K-cross:            %.1f%%\n", type_i_error_kcross * 100))

# ============================================================================
# SAVE RESULTS
# ============================================================================

out <- list(
  # Power metrics
  shade_image_power = shade_image_power,
  simple_shade_image_power = simple_shade_image_power,
  gcross_image_power = gcross_image_power,
  kcross_image_power = kcross_image_power,
  # Coverage metrics (SHADE only)
  shade_image_coverage = shade_image_coverage,
  simple_shade_image_coverage = simple_shade_image_coverage,
  # Type I error metrics
  type_i_error_shade = type_i_error_shade,
  type_i_error_simple = type_i_error_simple,
  type_i_error_gcross = type_i_error_gcross,
  type_i_error_kcross = type_i_error_kcross,
  # Condition parameters
  condition = condition
)
print(out)

saveRDS(out, file_out)

end <- Sys.time()
cat("\nTotal runtime:", round(difftime(end, start, units = "mins"), 1), "minutes\n")

end - start
