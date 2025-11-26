library(tidyverse)
library(spatstat)
library(Matrix)
library(SHADE)
library(tictoc)

# Load utilities
source("utils.R")

# Set environment and get data path
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  sim_idx <- 7
  n_cores <- 2
} else {
  args <- commandArgs(trailingOnly=TRUE)
  sim_idx <- as.numeric(args[1])
  n_cores <- 1
  print(paste("Simulation index:", sim_idx))
}

# Get data path (handles HPC vs local automatically)
path <- get_data_path("sim_timing")
print(paste("Data path:", path))

draws <- 1e3

# parameters (must match 01_generate_data.R)
grid <- expand.grid(
  total_cells = c(5000, 10000, 25000, 50000, 100000, 250000),
  sim = 1:20
)
total_cells <- grid$total_cells[sim_idx]
sim <- grid$sim[sim_idx]

num_types <- 3
num_pot <- 3
ratio <- 2
cells_per_type <- floor(total_cells / num_types)
n_dummy <- floor(cells_per_type * ratio)

file_pats <- paste0(path, "pats_sim", sim, "_cells_", total_cells, ".rds")
file_data_stan <- paste0(path, "data_stan_sim", sim, "_cells_", total_cells, ".json")
file_fit <- paste0(path, "fit_sim", sim, "_cells_", total_cells, ".rds")
file_timing <- paste0(path, "timing_sim", sim, "_cells_", total_cells, ".rds")

print(paste("Processing sim", sim, "with", total_cells, "total cells"))

# Load pattern data
pats <- readRDS(file_pats)
pat_full <- pats[[1]]$pat

# Prepare coordinates for feature construction timing
coords <- tibble(
  x = pat_full$x,
  y = pat_full$y,
  type = factor(paste0("t", as.character(pat_full$marks)), levels = paste0("t", 1:num_types)),
  image_id = "img_1"
)

image_ids <- "img_1"
patient_metadata <- tibble(
  Spot = image_ids,
  Patient = "pt_1",
  Group = "group_1"
)

# Time feature construction
print("Timing feature construction...")
tic()
prep <- prepare_spatial_model_data(
  x = coords$x,
  y = coords$y,
  cell_type = coords$type,
  image_id = factor(coords$image_id),
  patient_metadata = patient_metadata,
  type_idx = num_types,
  n_dummy = n_dummy,
  n_basis_functions = num_pot
)
time_feature <- toc(quiet = TRUE)
time_feature_sec <- time_feature$toc - time_feature$tic

print(paste("Feature construction time:", round(time_feature_sec, 3), "seconds"))

# Time full model fitting
print("Timing model fitting...")
tic()
fit <- run_SHADE_model(
  file_data_stan,
  method = "variational",
  draws = draws,
  refresh = 100
)
time_fitting <- toc(quiet = TRUE)
time_fitting_sec <- time_fitting$toc - time_fitting$tic

print(paste("Model fitting time:", round(time_fitting_sec, 3), "seconds"))

# Save model fit
# fit$save_object(file_fit)

# Save timing results
timing_results <- list(
  total_cells = total_cells,
  sim = sim,
  time_feature_sec = time_feature_sec,
  time_fitting_sec = time_fitting_sec,
  time_total_sec = time_feature_sec + time_fitting_sec
)
saveRDS(timing_results, file_timing)

print("Done!")
