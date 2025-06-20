library(spatstat)
library(tidyverse)
library(Matrix)
library(SHADE)

# set environment
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  path <- "./sim_control_confounding/data/"
  sim_idx <- 7
  # parameters to adjust
  num_pts_per_group <- 5
  images_per_pt <- 2
} else {
  path <- "./sim_control_confounding/data/"
  args <- commandArgs(trailingOnly=TRUE)
  sim_idx <- as.numeric(args[1])
  print(sim_idx)
  # parameters to adjust
  num_pts_per_group <- 20
  images_per_pt <- 2
}

sim <- sim_idx

# other parameters (likely don't need to adjust for now)
num_types <- 3
ratio <- 2
np <- 150
num_points_per_type <- rep(np,num_types)
n_dummy <- floor(np * ratio)
num_combos <- num_types - 1
num_points_per_type <- rep(np,num_types)

num_points_gen <- mean(num_points_per_type)
num_pt_groups <- 1
num_pts <- num_pts_per_group * num_pt_groups
size_im <- 1500
potentials <- make_rbfs(n_basis_functions = 3, max_dist = 75, basis_function_sigma = 15)
mean_alpha <- log(np/size_im^2)
sigma_beta_global <- 0.5
sigma_alpha_global <- 0.5
sigma_beta_indiv <- 0.1
sigma_alpha_indiv <- 0.5
sigma_beta_local <- 0.1
sigma_alpha_local <- 0.5
scale_sigmas <- 5
grainsize <- 1
intensity_resolution <- 128


# derived parameters
seed <- as.integer(2024 + sim_idx)
coef_seed <- as.integer(2024 + sim)
print(seed)
print(coef_seed)
num_images <- images_per_pt * num_pts
sample_to_indiv <- rep(1:num_pts,each = images_per_pt)
num_pts_per_group <- ceiling(num_pts/num_pt_groups)
indiv_to_group <- rep(1:num_pt_groups,each = num_pts_per_group,length.out=num_pts)
num_pot <- length(potentials)
scale_sigma_betas <- seq(5,1,length.out=num_pot)
file_data_stan <- paste0(path,"data_stan_sim",sim,".json")
file_ground_truth <- paste0(path,"ground_truth_sim",sim,".rds")
file_pats <- paste0(path,"pats_sim",sim,".rds")


# make logistic parameters
params <- make_simulation_parameters(mean_alpha,
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
                                  coef_seed)
betas_local <- params$betas_local[-1,]


# Modified simulation for confounding scenario
# T cells -> Tumor cells (TRUE interaction)
# T cells -> B cells (creates confounding)
# B cells appear to interact with Tumor cells (SPURIOUS)

set.seed(seed)
W <- owin(c(0,size_im),c(0,size_im))
area <- size_im^2
num_images <- num_pts * images_per_pt

# Set up 3 cell types: T cells (type 1), B cells (type 2), Tumor cells (type 3 - target)
num_types <- 3
num_points_per_type <- rep(np, num_types)

pats <- lapply(1:num_images, function(i) {
  
  # STEP 1: Generate T cells (independent baseline)
  T_cells <- rpoispp(lambda = num_points_gen/area, win = W)
  marks(T_cells) <- factor("T")  # T cells = source type 1
  
  # STEP 2: Generate B cells that are spatially correlated with T cells
  # Create B cells with some clustering around T cells (confounding correlation)
    
  # Moderate spatial correlation between T and B cells
  coeffs_B_corr <- betas_local[1:3, i] * 0.4  # Weaker than true interactions
  
  custom_kernel <- Vectorize(function(x, y) {
    d <- sqrt(x^2 + y^2)
    kernel_value <- sum(sapply(seq_along(coeffs_B_corr), function(k) {
      coeffs_B_corr[k] * potentials[[k]](d)
    }))
    return(kernel_value)
  })
  
  T_subset <- unmark(subset(T_cells, marks == "T"))
  T_dens_for_B <- smooth_density_fft(T_subset, custom_kernel, resolution = 128)
  # plot(dens)
  
  # Add significant background component so B cells aren't perfectly correlated
  background_rate <- 0.6  # 60% background, 40% T-correlated
  background_lambda <- (num_points_gen/area) * background_rate
  
  # Scale the T-correlated component
  lambda_integral_B <- sum(exp(T_dens_for_B$v)) * (T_dens_for_B$xstep * T_dens_for_B$ystep)
  T_corr_rate <- (num_points_gen) * (1 - background_rate)
  beta0_B <- log(T_corr_rate / lambda_integral_B)
  
  # Generate B cells: background + T-correlated component
  B_background <- rpoispp(lambda = background_lambda, win = W)
  B_Tcorrelated <- rpoispp(lambda = exp(T_dens_for_B + beta0_B))
  B_cells <- superimpose(B_background, B_Tcorrelated)
  marks(B_cells) <- factor("B")  # B cells = source type 2
  
  # STEP 3: Generate Tumor cells (target) - influenced ONLY by T cells, NOT B cells
  # This is the key: tumor cells respond to T cells but not B cells directly
    
  # TRUE interaction: T cells -> tumor cells (use full strength coefficients)
  coeffs_T_to_tumor <- betas_local[4:6, i]
  
  custom_kernel <- Vectorize(function(x, y) {
    d <- sqrt(x^2 + y^2)
    kernel_value <- sum(sapply(seq_along(coeffs_T_to_tumor), function(k) {
      coeffs_T_to_tumor[k] * potentials[[k]](d)
    }))
    return(kernel_value)
  })
  
  T_subset <- unmark(subset(T_cells, marks == "t1"))
  tumor_dens_combined <- smooth_density_fft(T_subset, custom_kernel, resolution = 128)
  
  # Generate tumor cells based ONLY on T cell density (B cells have NO direct effect)
  lambda_integral_tumor <- sum(exp(tumor_dens_combined$v)) * (tumor_dens_combined$xstep * tumor_dens_combined$ystep)
  beta0_tumor <- log(np / lambda_integral_tumor)
  
  Tumor_cells <- rpoispp(lambda = exp(tumor_dens_combined + beta0_tumor))
  marks(Tumor_cells) <- factor("tumor")  # Tumor cells = target type
  
  # Combine all cell types
  pat <- superimpose(T_cells, B_cells, Tumor_cells)
  
  # Store the true interaction structure for validation
  true_interactions <- list(
    T_to_tumor_coeffs = coeffs_T_to_tumor,    # TRUE interaction (should be detected)
    B_to_tumor_coeffs = rep(0, num_pot),      # NO direct interaction (should be ~0)
    T_B_correlation = coeffs_B_corr,          # Creates confounding
    confounding_strength = background_rate
  )
  
  print(paste("Image", i, "- T cells:", sum(marks(pat) == "T"), 
              "B cells:", sum(marks(pat) == "B"), 
              "Tumor cells:", sum(marks(pat) == "tumor")))
  
  list(
    pat = pat,
    beta0_tumor = beta0_tumor,
    true_interactions = true_interactions
  )
})


plot(pats[[8]]$pat)
table(pats[[7]]$pat$marks)

saveRDS(pats,file_pats)

beta0_local <- sapply(pats,\(o) o$beta0)
params$betas_local[1,] <- beta0_local

beta0_indiv <- sapply(1:num_pts,\(i) mean(beta0_local[which(sample_to_indiv == i)]))
beta0_global <- sapply(1:num_pt_groups,\(i) mean(beta0_indiv[which(indiv_to_group == i)]))

params$betas_indiv[1,] <- beta0_indiv
params$betas_global[1,] <- beta0_global

coords <- do.call(rbind, lapply(seq_along(pats), function(i) {
  pat <- pats[[i]]$pat
  
  tibble(
    x = pat$x,
    y = pat$y,
    type = marks(pat),
    image_id = paste0("img_", i)
  )
}))

image_ids <- unique(coords$image_id)
patient_ids <- rep(paste0("pt_", 1:num_pts_per_group), each = images_per_pt)
patient_metadata <- tibble(
  Spot = image_ids,
  Patient = patient_ids,
  Group = rep("group_1", length(image_ids))
)

prep <- prepare_spatial_model_data(
  x = coords$x,
  y = coords$y,
  cell_type = coords$type,
  image_id = coords$image_id,
  patient_metadata = patient_metadata,
  type_idx = num_types,
  n_dummy = n_dummy,
  n_basis_functions = num_pot
)
print("writing json...")
write_json_chunked(prep$stan_data,file_data_stan,chunk_size = 1e6)
print("json written.")
saveRDS(params,file_ground_truth)
