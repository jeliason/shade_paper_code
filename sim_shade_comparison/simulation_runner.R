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

# Set environment and get parameters
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  path <- "./sim_shade_comparison/data/"
  sim_idx <- 4
} else {
  path <- "./sim_shade_comparison/data/"
  args <- commandArgs(trailingOnly=TRUE)
  sim_idx <- as.numeric(args[1])
  cat("Running simulation", sim_idx, "\n")
}

file_out <- paste0(path,"sim_",sim_idx,".rds")


if (file.exists(file_out) & SYSTEM_ENV == "HPC") {
  cat("Simulation", sim_idx, "already completed. Skipping.\n")
  quit(save="no")
}

grid <- expand.grid(
  # imbalance_level = c("balanced", "moderate", "severe"),
  t_density = c("high", "low"),
  tumor_density = c("high", "low"),
  num_images = c(1,2,3),
  sim_rep = 1:30
)

# Extract current condition parameters
condition <- list(
  # imbalance_level = as.character(grid$imbalance_level[sim_idx]),
  num_images = as.character(grid$num_images[sim_idx]),
  t_density = as.character(grid$t_density[sim_idx]),
  tumor_density = as.character(grid$tumor_density[sim_idx]),
  sim_rep = grid$sim_rep[sim_idx]
)

cat("Condition:", condition$t_density, condition$tumor_density, condition$num_images,
    "Rep:", condition$sim_rep, "\n")

# Fixed parameters
num_pts_per_group <- 20
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
  num_pts_per_group = num_pts_per_group,
  num_images = grid$num_images[sim_idx],
  # imbalance_level = condition$imbalance_level,
  seed = 2024 + sim_idx
)

cat("Generating spatial coefficients...\n")
coefficients <- generate_spatial_coefficients(
  structure = structure,
  num_potentials = num_potentials,
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
shade_results <- run_shade_analysis(patterns, structure, num_potentials)

# Define potentials for simple SHADE
potentials <- make_rbfs(max_dist = 75, n_basis_functions = 3, basis_function_sigma = 15)

simple_shade_results <- run_simple_shade_analysis(patterns, potentials)

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
  true_group <- structure$images_df$group[i]

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
        lo_pw = as.vector(quantile(lp, probs = 0.025)),
        hi_pw = as.vector(quantile(lp, probs = 0.975)),
        lp_true = lp_true
      )

    # Check detection at target distances
    target_x <- c(20, 40, 60)
    target_indices <- sapply(target_x, function(x_val) which.min(abs(df_summ$x - x_val)))
    target_rows <- df_summ[target_indices, ]

    # Determine direction
    is_negative <- with(target_rows, hi_pw < 0)
    is_positive <- with(target_rows, lo_pw > 0)

    # Check if detection matches expected pattern
    if (true_group == "responder") {
      sum_pass <- sum(is_positive, na.rm = TRUE)
    } else {
      sum_pass <- sum(is_negative, na.rm = TRUE)
    }

    return(sum_pass)

  }, error = function(e) {
    cat("Error extracting simple SHADE results for image", i, ":", conditionMessage(e), "\n")
    return(NULL)
  })
})

# Calculate SHADE image-level power
shade_image_results <- sapply(1:structure$total_images, function(i) {
  true_group <- structure$images_df$group[i]

  # Get image-level coefficients for all 3 potentials
  beta <- as.vector(shade_draws$beta_local[i,2:4])
  beta_true <- coefficients$image_effects[i,]

  # Compute SIC using same distance range as global envelopes
  x_seq <- seq(DISTANCE_RANGE_MIN, DISTANCE_RANGE_MAX, 1)
  x_des <- lapply(potentials,\(pot) pot(x_seq)) %>% do.call(cbind,.)

  lp <- as.vector(x_des %*% beta)
  lp_true <- as.vector(x_des %*% beta_true)

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
      lp_true = lp_true
    )

  # Specify target x values
  target_x <- c(20, 40, 60)

  # Find indices of x values in df_summ closest to each target
  target_indices <- sapply(target_x, function(x_val) which.min(abs(df_summ$x - x_val)))

  # Subset the relevant rows
  target_rows <- df_summ[target_indices, ]

  # Determine direction
  is_negative <- with(target_rows, hi_pw < 0)
  is_positive <- with(target_rows, lo_pw > 0)

  # Check if detection matches expected pattern
  if (true_group == "responder") {
    sum_pass <- sum(is_positive,na.rm=TRUE)
  } else {
    sum_pass <- sum(is_negative, na.rm = TRUE)
  }

  return(sum_pass)
})

# Calculate performance metrics
shade_image_power <- mean(shade_image_results, na.rm = TRUE)
gx_detect <- unlist(lapply(gcross_results,\(o) o$mean_detect))
gcross_image_power <- mean(gx_detect, na.rm = TRUE)
kx_detect <- unlist(lapply(kcross_results,\(o) o$mean_detect))
kcross_image_power <- mean(kx_detect, na.rm = TRUE)
simple_shade_image_power <- mean(unlist(compact(simple_shade_image_results)))

out <- list(shade_image_power=shade_image_power,
            simple_shade_image_power=simple_shade_image_power,
            gcross_image_power=gcross_image_power,
            kcross_image_power=kcross_image_power)

saveRDS(out,file_out)

end <- Sys.time()

end - start
