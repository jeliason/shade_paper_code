library(spatstat)
library(tidyverse)
library(Matrix)
library(SHADE)
library(tictoc)

# Load utilities
source("utils.R")

# Set environment and get data path
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  sim_idx <- 1
} else {
  args <- commandArgs(trailingOnly=TRUE)
  sim_idx <- as.numeric(args[1])
  print(paste("Simulation index:", sim_idx))
}

# Get data path (handles HPC vs local automatically)
path <- get_data_path("sim_timing")
print(paste("Data path:", path))

# parameters to adjust
# cell counts: 1K, 2.5K, 5K, 10K, 25K, 50K, 100K
# 20 replicates each
# SLURM_ARRAY_SIZE: nrow(grid)
grid <- expand.grid(
  total_cells = c(5000, 10000, 25000, 50000, 100000, 250000),
  sim = 1:20
)
total_cells <- grid$total_cells[sim_idx]
sim <- grid$sim[sim_idx]

# fixed parameters
num_types <- 3
num_combos <- num_types - 1
ratio <- 2
size_im <- 1500
num_pot <- 3
potentials <- make_rbfs(n_basis_functions = num_pot, max_dist = 75, basis_function_sigma = 15)

# More realistic hierarchical structure
# 10 patients, 4 images per patient = 40 total images
num_pts <- 40
images_per_pt <- 4
num_images <- num_pts * images_per_pt
num_pt_groups <- 1
sample_to_indiv <- rep(1:num_pts, each = images_per_pt)
indiv_to_group <- rep(1, num_pts)

# cell counts per type per image
cells_per_image <- floor(total_cells / num_images)
cells_per_type <- floor(cells_per_image / num_types)
num_points_per_type <- rep(cells_per_type, num_types)
n_dummy <- floor(cells_per_type * ratio)
num_points_gen <- mean(num_points_per_type)

# intensity parameters
mean_alpha <- log(num_points_gen / size_im^2)
sigma_beta_global <- 0.5
sigma_alpha_global <- 0.5
sigma_beta_indiv <- 0.1
sigma_alpha_indiv <- 0.5
sigma_beta_local <- 0.1
sigma_alpha_local <- 0.5
scale_sigmas <- 5

# derived parameters
seed <- as.integer(3000 + sim_idx)
coef_seed <- as.integer(3000 + sim)
print(paste("Seed:", seed, "Coef seed:", coef_seed))
print(paste("Total cells:", total_cells, "Cells per image:", cells_per_image))

file_data_stan <- paste0(path, "data_stan_sim", sim, "_cells_", total_cells, ".json")
file_ground_truth <- paste0(path, "ground_truth_sim", sim, "_cells_", total_cells, ".rds")
file_pats <- paste0(path, "pats_sim", sim, "_cells_", total_cells, ".rds")

# make logistic parameters
params <- make_simulation_parameters(
  mean_alpha,
  sigma_beta_global,
  sigma_beta_indiv,
  sigma_beta_local,
  scale_sigmas,
  num_pt_groups,
  num_types,
  num_combos,
  num_pot,
  indiv_to_group,
  num_pts,
  num_images,
  sample_to_indiv,
  coef_seed
)
betas_local <- params$betas_local[-1,]

set.seed(seed)
W <- owin(c(0, size_im), c(0, size_im))
area <- size_im^2

print("Generating spatial point patterns...")
tic()
pats <- lapply(1:num_images, \(i) {
  # predict the last type from all of the others
  if(num_types == 2) {
    pat <- rpoispp(lambda = num_points_gen / area, win = W)
    marks(pat) <- factor("t1")
  } else {
    pat <- rmpoispp(lambda = rep(num_points_gen / area, num_types - 1), win = W)
  }

  dens <- lapply(1:(num_types - 1), \(j) {
    start <- (j - 1) * num_pot + 1
    stop <- j * num_pot
    coeffs <- betas_local[start:stop, i]
    custom_kernel <- Vectorize(function(x, y) {
      d <- sqrt(x^2 + y^2)
      kernel_value <- sum(sapply(seq_along(coeffs), function(i) {
        coeffs[i] * potentials[[i]](d)
      }))
      return(kernel_value)
    })

    subs <- unmark(subset(pat, marks == j))
    dens <- smooth_density_fft(subs, custom_kernel, resolution = 128)
  })

  dens <- Reduce("+", dens)
  lambda_integral <- sum(exp(dens$v)) * (dens$xstep * dens$ystep)

  # Compute beta0 to achieve expected cell count
  beta0 <- log(num_points_gen / lambda_integral)
  pat2 <- rpoispp(lambda = exp(dens + beta0))
  marks(pat2) <- factor(num_types)
  print(i)

  pat <- superimpose(pat, pat2)
  list(pat = pat, beta0 = beta0)
})
time_pattern_gen <- toc(quiet = TRUE)
time_pattern_gen_sec <- time_pattern_gen$toc - time_pattern_gen$tic

print(paste("Pattern generation time:", round(time_pattern_gen_sec, 3), "seconds"))
print(paste("Total cells across all images:", sum(sapply(pats, \(p) npoints(p$pat)))))

saveRDS(pats, file_pats)

beta0_local <- sapply(pats, \(o) o$beta0)
params$betas_local[1,] <- beta0_local

beta0_indiv <- sapply(1:num_pts, \(i) mean(beta0_local[which(sample_to_indiv == i)]))
beta0_global <- sapply(1:num_pt_groups, \(i) mean(beta0_indiv[which(indiv_to_group == i)]))

params$betas_indiv[1,] <- beta0_indiv
params$betas_global[1,] <- beta0_global

coords <- do.call(rbind, lapply(seq_along(pats), function(i) {
  pat <- pats[[i]]$pat
  raw_marks <- as.character(pat$marks)
  relabeled_marks <- paste0("t", raw_marks)

  tibble(
    x = pat$x,
    y = pat$y,
    type = factor(relabeled_marks, levels = paste0("t", 1:num_types)),
    image_id = paste0("img_", i)
  )
}))

image_ids <- unique(coords$image_id)
patient_ids <- rep(paste0("pt_", 1:num_pts), each = images_per_pt)
patient_metadata <- tibble(
  Spot = image_ids,
  Patient = patient_ids,
  Group = rep("group_1", length(image_ids))
)

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

print("Writing JSON...")
tic()
write_json_chunked(prep$stan_data, file_data_stan, chunk_size = 1e6)
time_json_write <- toc(quiet = TRUE)
time_json_write_sec <- time_json_write$toc - time_json_write$tic
print(paste("JSON written in", round(time_json_write_sec, 3), "seconds"))

saveRDS(params, file_ground_truth)

# Save timing results for data generationd
file_timing_gen <- paste0(path, "timing_gen_sim", sim, "_cells_", total_cells, ".rds")
timing_gen_results <- list(
  total_cells = total_cells,
  sim = sim,
  time_pattern_gen_sec = time_pattern_gen_sec,
  time_json_write_sec = time_json_write_sec,
  time_total_sec = time_pattern_gen_sec + time_json_write_sec
)
saveRDS(timing_gen_results, file_timing_gen)
