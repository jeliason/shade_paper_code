library(tidyverse)
library(spatstat)
library(cmdstanr)
library(Matrix)
library(posterior)
library(ggdist)
library(survival)
library(SHADE)

# Load utility functions and constants
source("../utils.R")

# Set environment and get data path
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  sim_indices <- c(1,38,50)
} else {
  sim_indices <- 1:135
}

# Get data path (handles HPC vs local automatically)
path <- get_data_path("sim_sample_size")
print(paste("Data path:", path))

num_types <- 3
num_pot <- 3
types_grid <- expand.grid(t1=num_types,t2=1:(num_types-1))
num_combos <- nrow(types_grid)
combo_pot <- num_combos * num_pot
grid <- expand.grid(num_pts_per_group=c(10,20,40),
                    images_per_pt=c(1,2,4),
                    sim=1:15)

grid <- grid[sim_indices,]

# RMSE
rmse_tb <- pmap(grid,\(num_pts_per_group,images_per_pt,sim) {
  file_fit <- paste0(path,"fit_sim",sim,"_pts_",num_pts_per_group,"_ims_",images_per_pt,".rds")
  file_ground_truth <- paste0(path,"ground_truth_sim",sim,"_pts_",num_pts_per_group,"_ims_",images_per_pt,".rds")
  # print(file_fit)
  fit <- readRDS(file_fit)
  ground_truth <- readRDS(file_ground_truth)
  
  # print("read in files")
  # global
  estimate <- as_draws_rvars(fit$draws(variables = "beta_global"))$beta_global
  estimate <-  estimate[-1,] # remove alphas
  true <- ground_truth$betas_global[-1,]
  
  # print(dim(estimate))
  # print(dim(true))
  rmse <- rvar_mean((estimate - true)^2)
  rmse_close <- rvar_mean(((estimate - true)^2)[seq(1,combo_pot,by=num_pot),])
  rmse_mid <- rvar_mean(((estimate - true)^2)[seq(2,combo_pot,by=num_pot),])
  rmse_far <- rvar_mean(((estimate - true)^2)[seq(3,combo_pot,by=num_pot),])
  
  df <- tibble(rmse=c(rmse,rmse_close,rmse_mid,rmse_far),
               level="beta_global",
               scale=c("all","close","mid","far"))
  
  # print("global rmse")
  
  # indiv
  estimate <- as_draws_rvars(fit$draws(variables = "beta_indiv"))$beta_indiv
  estimate <-  estimate[-1,] # remove alphas
  
  true <- ground_truth$betas_indiv[-1,]
  rmse <- rvar_mean((estimate - true)^2)
  rmse_close <- rvar_mean(((estimate - true)^2)[seq(1,combo_pot,by=num_pot),])
  rmse_mid <- rvar_mean(((estimate - true)^2)[seq(2,combo_pot,by=num_pot),])
  rmse_far <- rvar_mean(((estimate - true)^2)[seq(3,combo_pot,by=num_pot),])
  
  df <- df %>%
    bind_rows(tibble(rmse=c(rmse,rmse_close,rmse_mid,rmse_far),
                     level="beta_indiv",
                     scale=c("all","close","mid","far")))
  
  # print("indiv rmse")
  
  # local
  estimate <- t(as_draws_rvars(fit$draws(variables = "beta_local"))$beta_local)
  estimate <-  estimate[-1,] # remove alphas
  
  true <- ground_truth$betas_local[-1,]
  rmse <- rvar_mean((estimate - true)^2)
  rmse_close <- rvar_mean(((estimate - true)^2)[seq(1,combo_pot,by=num_pot),])
  rmse_mid <- rvar_mean(((estimate - true)^2)[seq(2,combo_pot,by=num_pot),])
  rmse_far <- rvar_mean(((estimate - true)^2)[seq(3,combo_pot,by=num_pot),])
  
  df <- df %>%
    bind_rows(tibble(rmse=c(rmse,rmse_close,rmse_mid,rmse_far),
                     level="beta_local",
                     scale=c("all","close","mid","far")))
  
  # print("local rmse")
  
  
  # print("gamma")
  df <- df %>%
    mutate(num_pts_per_group=num_pts_per_group,images_per_pt=images_per_pt,sim=sim)
  
}) %>%
  bind_rows() %>%
  rename(coefficient=level)

### plot example SICs estimate vs ground truth
sic_tb <- grid %>%
  pmap(\(num_pts_per_group,images_per_pt,sim) {
    file_fit <- paste0(path,"fit_sim",sim,"_pts_",num_pts_per_group,"_ims_",images_per_pt,".rds")
    file_ground_truth <- paste0(path,"ground_truth_sim",sim,"_pts_",num_pts_per_group,"_ims_",images_per_pt,".rds")
    potentials <- make_rbfs(n_basis_functions = 3, max_dist = 75, basis_function_sigma = 15)
    # print(file_fit)
    fit <- readRDS(file_fit)
    ground_truth <- readRDS(file_ground_truth)
    
    # global
    estimate <- as_draws_rvars(fit$draws(variables = "beta_global"))$beta_global
    beta_estimate <-  estimate[-1,] # remove alphas
    
    beta_true <- ground_truth$betas_global[-1,]
    # cbind(beta_true[,1],beta_estimate[,1])
    x_seq <- seq(0,100,1)
    x_des <- lapply(potentials,\(pot) pot(x_seq)) %>% do.call(cbind,.)
    ix <- 4:6
    b_t <- as.matrix(beta_true)[ix,]
    b_e <- as.matrix(beta_estimate)[ix,]
    
    lp_t <- x_des %*% b_t
    lp_e <- x_des %*% b_e
    
    
    lp <- bind_cols(true=lp_t[,1],estimate=as.vector(lp_e[,1]))

    # Compute simultaneous 95% credible bands using utility function
    bands <- compute_simultaneous_bands(
      lp_data = lp %>% select(estimate),
      x_seq = x_seq,
      alpha = 0.05
    )

    # Combine with true values and reshape for plotting
    lp <- bands %>%
      rename(`Post. Mean` = mean) %>%
      mutate(True = as.vector(lp[[1]])) %>%
      pivot_longer(c(`Post. Mean`, True), names_to = "Curve", values_to = "value") %>%
      mutate(
        lo = ifelse(Curve == "Post. Mean", lower, NA_real_),
        hi = ifelse(Curve == "Post. Mean", upper, NA_real_)
      ) %>%
      select(-variable, -lower, -upper)
    
    lp <- lp %>%
      mutate(num_pts_per_group=num_pts_per_group,images_per_pt=images_per_pt,sim=sim)
    
  }) %>% bind_rows()


# Save standardized analysis summary for fetching/plotting
analysis_summary <- list(
  rmse_tb = rmse_tb,
  sic_tb = sic_tb,
  metadata = list(
    sim_dir = "sim_sample_size",
    date_analyzed = Sys.time(),
    n_simulations = length(sim_indices)
  )
)

saveRDS(analysis_summary, paste0(path, "analysis_summary.rds"))
print(paste("Saved analysis summary to:", paste0(path, "analysis_summary.rds")))
