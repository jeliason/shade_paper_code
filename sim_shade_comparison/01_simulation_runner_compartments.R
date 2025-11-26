library(spatstat)
library(tidyverse)
library(Matrix)
library(SHADE)
library(cmdstanr)

# Load utility functions and constants
source("./utils.R")
source("sim_shade_comparison/comparison_functions.R")
source("sim_shade_comparison/plot_confounder_diagnostics.R")
start <- Sys.time()

# ============================================================================
# SCENARIO 3: TISSUE COMPARTMENT CONFOUNDER
# Tests SHADE performance when discrete spatial regions have different
# baseline target densities (e.g., tumor islands, stroma, crypts)
# ============================================================================

SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  sim_idx <- 5
  verbose <- TRUE
} else {
  args <- commandArgs(trailingOnly=TRUE)
  sim_idx <- as.numeric(args[1])
  verbose <- FALSE
  cat("Running compartment confounder simulation", sim_idx, "\n")
}

path <- get_data_path("sim_shade_comparison")
compartment_path <- paste0(path, "compartments/")
dir.create(compartment_path, showWarnings = FALSE, recursive = TRUE)
file_out <- paste0(compartment_path, "sim_", sim_idx, ".rds")

force_rerun <- Sys.getenv("FORCE_RERUN", unset = "0") == "1"
if (file.exists(file_out) & SYSTEM_ENV == "HPC" & !force_rerun) {
  cat("Simulation", sim_idx, "already completed. Skipping.\n")
  quit(save="no")
}

# SIMULATION GRID
# Test all four density regimes to see how compartment confounding affects each
# Vary compartment effect strength to test confounding severity
grid <- expand.grid(
  method = c("variational"),
  parameterization = c("centered"),
  scale_sigma = c(5),
  t_density = c("high", "low"),
  tumor_density = c("high", "low"),
  num_images = 3,  # Stable regime from main results
  n_compartments = 3,  # Moderate complexity
  compartment_effect = c(0.8, 1.2, 1.5),  # Weak, moderate, strong confounding
  sim_rep = 1:100
)

condition <- list(
  method = as.character(grid$method[sim_idx]),
  parameterization = as.character(grid$parameterization[sim_idx]),
  scale_sigma = grid$scale_sigma[sim_idx],
  num_images = as.numeric(grid$num_images[sim_idx]),
  t_density = as.character(grid$t_density[sim_idx]),
  tumor_density = as.character(grid$tumor_density[sim_idx]),
  n_compartments = grid$n_compartments[sim_idx],
  compartment_effect = grid$compartment_effect[sim_idx],
  sim_rep = grid$sim_rep[sim_idx]
)

cat("Compartment Confounder Simulation:\n")
cat("  N compartments:", condition$n_compartments, "Effect strength:", condition$compartment_effect, "\n")
cat("  Densities:", condition$t_density, "/", condition$tumor_density,
    "Images:", condition$num_images, "Rep:", condition$sim_rep, "\n")

# Fixed parameters
num_patients <- 40
num_potentials <- 3
size_im <- 1500

density_params <- list(
  np_t = ifelse(condition$t_density == "high", 150, 15),
  np_tumor = ifelse(condition$tumor_density == "high", 150, 15),
  np_b = ifelse(condition$t_density == "high", 150, 15)
)

# ============================================================================
# RUN SIMULATION
# ============================================================================

cat("Creating patient structure...\n")
structure <- create_patient_structure(
  num_patients = num_patients,
  num_images = condition$num_images,
  seed = 2024 + sim_idx
)

cat("Generating spatial coefficients...\n")
coefficients <- generate_spatial_coefficients(
  structure = structure,
  num_potentials = num_potentials,
  cohort_mean = c(1.5, 1.0, 0.5),
  seed = 2024 + condition$sim_rep
)

cat("Generating spatial patterns with COMPARTMENT confounder...\n")
cat("  Compartments:", condition$n_compartments, "effect:", condition$compartment_effect, "\n")
patterns <- generate_spatial_patterns_compartments(
  structure = structure,
  coefficients = coefficients,
  density_params = density_params,
  size_im = size_im,
  n_compartments = condition$n_compartments,
  compartment_effect = condition$compartment_effect,
  seed = 2024 + sim_idx + 1000
)

# plot_confounder_diagnostics(patterns,pattern_index = 4,confounder_type = "compartment")


# plot(patterns[[1]]$pattern)

cat("Running SHADE analysis (model ignores compartment structure)...\n")
shade_results <- run_shade_analysis(
  patterns,
  structure,
  num_potentials,
  method = condition$method,
  parameterization = condition$parameterization,
  scale_sigma = condition$scale_sigma
)

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

    # Compute SIC using same distance range as global envelopes
    x_seq <- seq(DISTANCE_RANGE_MIN, DISTANCE_RANGE_MAX, 1)
    x_des <- lapply(potentials, \(pot) pot(x_seq)) %>% do.call(cbind,.)

    lp <- as.vector(x_des %*% beta)
    lp_true <- as.vector(x_des %*% beta_true)

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
        lp_true = lp_true
      )

    # Check if simultaneous credible band excludes zero anywhere in the range
    excludes_zero <- any(df_summ$lo_simul > 0 | df_summ$hi_simul < 0, na.rm = TRUE)

    # Check if simultaneous credible band contains true SIC everywhere (coverage)
    coverage <- all(df_summ$lo_simul <= df_summ$lp_true &
                    df_summ$hi_simul >= df_summ$lp_true, na.rm = TRUE)

    return(list(detection = as.numeric(excludes_zero),
                coverage = as.numeric(coverage)))

  }, error = function(e) {
    cat("Error extracting simple SHADE results for image", i, ":", conditionMessage(e), "\n")
    return(NULL)
  })
}, simplify = FALSE)

# Calculate SHADE hierarchical image-level power
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

  # Compute SIC using same distance range as global envelopes
  x_seq <- seq(DISTANCE_RANGE_MIN, DISTANCE_RANGE_MAX, 1)
  x_des <- lapply(potentials, \(pot) pot(x_seq)) %>% do.call(cbind,.)

  lp <- as.vector(x_des %*% beta)
  lp_true <- as.vector(x_des %*% beta_true)

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
      lp_true = lp_true
    )

  # Check if simultaneous credible band excludes zero anywhere in the range
  excludes_zero <- any(df_summ$lo_simul > 0 | df_summ$hi_simul < 0, na.rm = TRUE)

  # Check if simultaneous credible band contains true SIC everywhere (coverage)
  coverage <- all(df_summ$lo_simul <= df_summ$lp_true &
                  df_summ$hi_simul >= df_summ$lp_true, na.rm = TRUE)

  return(list(detection = as.numeric(excludes_zero),
              coverage = as.numeric(coverage)))
}, simplify = FALSE)

# Calculate performance metrics
shade_detections <- sapply(shade_image_results, function(x) if(is.list(x)) x$detection else NA)
shade_coverages <- sapply(shade_image_results, function(x) if(is.list(x)) x$coverage else NA)

simple_shade_detections <- sapply(compact(simple_shade_image_results), function(x) if(is.list(x)) x$detection else NA)
simple_shade_coverages <- sapply(compact(simple_shade_image_results), function(x) if(is.list(x)) x$coverage else NA)

# Power = proportion of images with detection
shade_image_power <- mean(shade_detections, na.rm = TRUE)
simple_shade_image_power <- mean(simple_shade_detections, na.rm = TRUE)

# Coverage = proportion of images where credible band contains true SIC
shade_image_coverage <- mean(shade_coverages, na.rm = TRUE)
simple_shade_image_coverage <- mean(simple_shade_coverages, na.rm = TRUE)

# G-cross and K-cross (detection only)
gx_detect <- unlist(lapply(gcross_results, \(o) o$mean_detect))
gcross_image_power <- mean(gx_detect, na.rm = TRUE)

kx_detect <- unlist(lapply(kcross_results, \(o) o$mean_detect))
kcross_image_power <- mean(kx_detect, na.rm = TRUE)

cat("\nCompartment Confounder Results:\n")
cat(sprintf("  SHADE Hierarchical Power: %.1f%%, Coverage: %.1f%%\n",
            shade_image_power * 100, shade_image_coverage * 100))
cat(sprintf("  SHADE Flat Power: %.1f%%, Coverage: %.1f%%\n",
            simple_shade_image_power * 100, simple_shade_image_coverage * 100))
cat(sprintf("  G-cross Power: %.1f%%\n", gcross_image_power * 100))
cat(sprintf("  K-cross Power: %.1f%%\n", kcross_image_power * 100))

# ============================================================================
# NULL CALIBRATION (Type I Error)
# ============================================================================
# Test: Does SHADE falsely detect interaction when only compartment effect exists
# but NO true source-target interaction?

cat("\n=== NULL CALIBRATION (Type I Error) ===\n")
cat("Generating NULL patterns (compartment effect only, no source-target interaction)...\n")

# Generate NULL coefficients (no spatial interaction)
coefficients_null <- generate_spatial_coefficients(
  structure = structure,
  num_potentials = num_potentials,
  cohort_mean = c(0, 0, 0),  # NULL: no spatial interaction
  sigma_individual = 0,       # NULL: no patient-level variation
  sigma_image = 0,            # NULL: no image-level variation
  seed = 2024 + condition$sim_rep + 5000
)

# Generate NULL patterns with compartment effect but no source-target interaction
patterns_null <- generate_spatial_patterns_compartments(
  structure = structure,
  coefficients = coefficients_null,  # Zero interaction coefficients
  density_params = density_params,
  size_im = size_im,
  n_compartments = condition$n_compartments,
  compartment_effect = condition$compartment_effect,
  seed = 2024 + sim_idx + 10000
)

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
  patient_id <- structure$images_df$patient_id[i]
  num_images_for_patient <- sum(structure$images_df$patient_id == patient_id)

  if(num_images_for_patient == 1) {
    beta <- as.vector(shade_draws_null$beta_indiv[2:4,patient_id])
  } else {
    beta <- as.vector(shade_draws_null$beta_local[i,2:4])
  }

  x_seq <- seq(DISTANCE_RANGE_MIN, DISTANCE_RANGE_MAX, 1)
  x_des <- lapply(potentials, \(pot) pot(x_seq)) %>% do.call(cbind,.)

  lp <- as.vector(x_des %*% beta)
  lp_data <- lp %>% as.matrix() %>% as.data.frame()
  colnames(lp_data) <- "SIC"

  df_summ <- compute_simultaneous_bands(lp_data = lp_data, x_seq = x_seq, alpha = 0.05) %>%
    rename(mn = mean, lo_simul = lower, hi_simul = upper)

  excludes_zero <- any(df_summ$lo_simul > 0 | df_summ$hi_simul < 0, na.rm = TRUE)
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
    rename(mn = mean, lo_simul = lower, hi_simul = upper)

  excludes_zero <- any(df_summ$lo_simul > 0 | df_summ$hi_simul < 0, na.rm = TRUE)
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

cat("\nType I Error Rates (with compartment confounder):\n")
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
  # Type I error metrics (with compartment confounder)
  type_i_error_shade = type_i_error_shade,
  type_i_error_simple = type_i_error_simple,
  type_i_error_gcross = type_i_error_gcross,
  type_i_error_kcross = type_i_error_kcross,
  # Condition parameters
  condition = condition,
  note = "Compartment confounder: discrete spatial regions with different baseline densities"
)

saveRDS(out, file_out)
end <- Sys.time()
cat("\nTotal time:", round(difftime(end, start, units = "secs"), 2), "seconds\n")
