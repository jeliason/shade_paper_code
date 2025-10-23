#!/usr/bin/env Rscript
# Null Calibration Simulation for SHADE vs Classical Methods
#
# Purpose: Test type I error rates under the null hypothesis of no spatial association
# Expected: All methods should have ~5% false positive rate at α=0.05

library(spatstat)
library(tidyverse)
library(Matrix)
library(SHADE)
library(cmdstanr)
library(posterior)

# Load utility functions
source("./utils.R")

# Load shared comparison functions
source("sim_shade_comparison/comparison_functions.R")

start <- Sys.time()

# ============================================================================
# CONFIGURATION
# ============================================================================

set.seed(2025)

# Simulation parameters
# Total images = num_pts_per_group × 2 groups × num_images_per_patient
# = 50 × 2 × 4 = 400 images for stable type I error estimation
num_pts_per_group <- 50
num_images_per_patient <- 4
size_im <- 1500
W <- owin(c(0, size_im), c(0, size_im))
area <- size_im^2

# Cell densities (using high density for better power to detect false positives)
np_t <- 150      # T cells
np_tumor <- 150  # Tumor cells
np_b <- 150      # B cells

# Create RBF potentials (same as main simulation)
num_potentials <- 3
potentials <- make_rbfs(n_basis_functions = num_potentials, max_dist = 75, basis_function_sigma = 15)

# ============================================================================
# GENERATE NULL DATA
# ============================================================================

cat("=== Generating null scenarios (no spatial association) ===\n\n")

# Create patient structure
structure <- create_patient_structure(
  num_pts_per_group = num_pts_per_group,
  num_images = num_images_per_patient,
  seed = 2024
)

# Generate NULL spatial coefficients (all zeros = no spatial interaction)
coefficients <- list(
  group_means = matrix(0, nrow = 2, ncol = num_potentials,
                      dimnames = list(c("responder", "nonresponder"),
                                     paste0("potential_", 1:num_potentials))),
  individual_effects = matrix(0, nrow = structure$total_patients, ncol = num_potentials),
  image_effects = matrix(0, nrow = structure$total_images, ncol = num_potentials)
)

# Generate spatial patterns with NO spatial association
# Target cells are generated uniformly (CSR), independent of source cells
cat("Generating", structure$total_images, "null patterns...\n")

patterns <- vector("list", structure$total_images)

for(i in 1:structure$total_images) {
  if(i %% 20 == 0) cat("  Pattern", i, "of", structure$total_images, "\n")

  # Generate source cells (T cells and B cells) - homogeneous Poisson
  pat_sources <- rmpoispp(lambda = c(np_t/area, np_b/area), win = W)

  # Generate tumor cells UNIFORMLY (NULL: no dependence on source cells)
  tumor_cells <- rpoispp(lambda = np_tumor/area, win = W)
  marks(tumor_cells) <- factor(3)

  # Combine all cell types
  combined_pattern <- superimpose(pat_sources, tumor_cells)

  patterns[[i]] <- list(
    pattern = combined_pattern,
    beta0 = log(np_tumor/area),
    image_coeffs = rep(0, num_potentials),  # NULL coefficients
    patient_id = structure$images_df$patient_id[i],
    group = structure$images_df$group[i]
  )
}

cat("\n✓ Generated", length(patterns), "null patterns\n\n")

# ============================================================================
# RUN ANALYSES
# ============================================================================

cat("=== Running SHADE analysis ===\n")
shade_results <- run_shade_analysis(patterns, structure, num_potentials)

cat("\n=== Running simple SHADE analysis ===\n")
simple_shade_file <- "sim_shade_comparison/data/null_calibration_simple_shade.rds"
if(file.exists(simple_shade_file)) {
  cat("Loading cached simple SHADE results...\n")
  simple_shade_results <- readRDS(simple_shade_file)
} else {
  simple_shade_results <- run_simple_shade_analysis(patterns, potentials)
  cat("Saving simple SHADE results...\n")
  saveRDS(simple_shade_results, simple_shade_file)
}

cat("\n=== Running G-cross analysis ===\n")
gcross_file <- "sim_shade_comparison/data/null_calibration_gcross.rds"
if(file.exists(gcross_file)) {
  cat("Loading cached G-cross results...\n")
  gcross_results <- readRDS(gcross_file)
} else {
  gcross_results <- run_gcross_analysis(patterns, structure)
  cat("Saving G-cross results...\n")
  saveRDS(gcross_results, gcross_file)
}

cat("\n=== Running K-cross analysis ===\n")
kcross_file <- "sim_shade_comparison/data/null_calibration_kcross.rds"
if(file.exists(kcross_file)) {
  cat("Loading cached K-cross results...\n")
  kcross_results <- readRDS(kcross_file)
} else {
  kcross_results <- run_kcross_analysis(patterns, structure)
  cat("Saving K-cross results...\n")
  saveRDS(kcross_results, kcross_file)
}

# ============================================================================
# CALCULATE TYPE I ERROR RATES
# ============================================================================

cat("\n=== Calculating Type I Error Rates ===\n\n")

# Extract SHADE results
shade_draws <- as_draws_rvars(shade_results$fit$draws())

# Check SHADE hierarchical: how often do simultaneous bands exclude zero?
cat("Checking SHADE hierarchical...\n")
shade_false_positives <- sapply(1:structure$total_images, function(i) {
  if(i %% 50 == 0) cat("  Image", i, "of", structure$total_images, "\n")
  beta <- as.vector(shade_draws$beta_local[i,2:4])

  # Use same distance range as global envelopes
  x_seq <- seq(DISTANCE_RANGE_MIN, DISTANCE_RANGE_MAX, 1)
  x_des <- lapply(potentials, \(pot) pot(x_seq)) %>% do.call(cbind,.)
  lp <- as.vector(x_des %*% beta)

  lp_data <- lp %>% as.matrix() %>% as.data.frame()
  colnames(lp_data) <- "SIC"

  df_summ <- compute_simultaneous_bands(
    lp_data = lp_data,
    x_seq = x_seq,
    alpha = 0.05
  )

  # Band excludes zero if lower > 0 OR upper < 0 anywhere in range
  excludes_zero <- any(df_summ$lower > 0) | any(df_summ$upper < 0)

  return(as.numeric(excludes_zero))
})

# Check simple SHADE
cat("\nChecking simple SHADE...\n")
simple_shade_false_positives <- sapply(1:structure$total_images, function(i) {
  if(i %% 10 == 0) cat("  Image", i, "of", structure$total_images, "\n")
  results <- simple_shade_results[[i]]
  if(is.null(results)) return(NA)

  shade_draws <- as_draws_rvars(results$fit$draws())
  beta <- as.vector(shade_draws$beta[1:3])

  # Use same distance range as global envelopes
  x_seq <- seq(DISTANCE_RANGE_MIN, DISTANCE_RANGE_MAX, 1)
  x_des <- lapply(potentials, \(pot) pot(x_seq)) %>% do.call(cbind,.)
  lp <- as.vector(x_des %*% beta)

  lp_data <- lp %>% as.matrix() %>% as.data.frame()
  colnames(lp_data) <- "SIC"

  df_summ <- compute_simultaneous_bands(
    lp_data = lp_data,
    x_seq = x_seq,
    alpha = 0.05
  )

  # Band excludes zero if lower > 0 OR upper < 0 anywhere in range
  excludes_zero <- any(df_summ$lower > 0) | any(df_summ$upper < 0)

  return(as.numeric(excludes_zero))
})

# Check G-cross: how often does global envelope show significant deviation?
# mean_detect is 1 if envelope was violated anywhere in 10-75 μm range, 0 otherwise
gcross_false_positives <- sapply(gcross_results, function(res) {
  if(is.na(res$mean_detect)) return(NA)
  return(res$mean_detect)
})

# Check K-cross: same logic
kcross_false_positives <- sapply(kcross_results, function(res) {
  if(is.na(res$mean_detect)) return(NA)
  return(res$mean_detect)
})

# Calculate type I error rates
type_i_error_shade <- mean(shade_false_positives, na.rm = TRUE)
type_i_error_simple <- mean(simple_shade_false_positives, na.rm = TRUE)
type_i_error_gcross <- mean(gcross_false_positives, na.rm = TRUE)
type_i_error_kcross <- mean(kcross_false_positives, na.rm = TRUE)

# ============================================================================
# REPORT RESULTS
# ============================================================================

cat("\n=== TYPE I ERROR RATES (Expected: ~5%) ===\n\n")
cat(sprintf("SHADE (hierarchical):  %.1f%% (%d/%d)\n",
            type_i_error_shade * 100,
            sum(shade_false_positives, na.rm=TRUE),
            sum(!is.na(shade_false_positives))))
cat(sprintf("SHADE (flat):          %.1f%% (%d/%d)\n",
            type_i_error_simple * 100,
            sum(simple_shade_false_positives, na.rm=TRUE),
            sum(!is.na(simple_shade_false_positives))))
cat(sprintf("G-cross:               %.1f%% (%d/%d)\n",
            type_i_error_gcross * 100,
            sum(gcross_false_positives, na.rm=TRUE),
            sum(!is.na(gcross_false_positives))))
cat(sprintf("K-cross:               %.1f%% (%d/%d)\n",
            type_i_error_kcross * 100,
            sum(kcross_false_positives, na.rm=TRUE),
            sum(!is.na(kcross_false_positives))))

# Save results
results <- list(
  type_i_error_shade = type_i_error_shade,
  type_i_error_simple = type_i_error_simple,
  type_i_error_gcross = type_i_error_gcross,
  type_i_error_kcross = type_i_error_kcross,
  shade_false_positives = shade_false_positives,
  simple_shade_false_positives = simple_shade_false_positives,
  gcross_false_positives = gcross_false_positives,
  kcross_false_positives = kcross_false_positives,
  n_images = structure$total_images,
  config = list(
    num_pts_per_group = num_pts_per_group,
    num_images_per_patient = num_images_per_patient,
    np_t = np_t,
    np_tumor = np_tumor,
    np_b = np_b
  )
)

saveRDS(results, "sim_shade_comparison/data/null_calibration_results.rds")

cat("\n✓ Results saved to sim_shade_comparison/data/null_calibration_results.rds\n")

end <- Sys.time()
cat("\nTotal runtime:", round(difftime(end, start, units = "mins"), 1), "minutes\n")

cat("\n=== Calibration check complete ===\n")
