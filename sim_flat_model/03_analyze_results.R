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
  sim_indices <- c(1,2,3)
} else {
  sim_indices <- 1:100
}

# Get data path (handles HPC vs local automatically)
path <- get_data_path("sim_flat_model")
print(paste("Data path:", path))

num_types <- 3
num_pot <- 3
types_grid <- expand.grid(t1=num_types,t2=1:(num_types-1))
num_combos <- nrow(types_grid)
combo_pot <- num_combos * num_pot
grid <- expand_grid(sim=1:100)

grid <- grid[sim_indices,]

# RMSE
rmse_tb <- pmap(grid,\(sim) {
  file_fit_hier <- paste0(path,"fit_sim",sim,"_hier.rds")
  file_fit_no_hier <- paste0(path,"fit_sim",sim,"_flat.rds")
  file_ground_truth <- paste0(path,"ground_truth_sim",sim,".rds")
  # print(file_fit)
  fit_hier <- readRDS(file_fit_hier)
  fit_no_hier <- readRDS(file_fit_no_hier)
  ground_truth <- readRDS(file_ground_truth)
  
  # local
  estimate_hier <- t(as_draws_rvars(fit_hier$draws(variables = "beta_local"))$beta_local)[-1,]

  estimate_no_hier <- as_draws_rvars(fit_no_hier$draws(variables = "beta_local"))$beta_local[-1,]

  true <- ground_truth$betas_local[-1,]
  rmse <- rvar_mean((estimate_hier - true)^2)
  rmse_close <- rvar_mean(((estimate_hier - true)^2)[seq(1,combo_pot,by=num_pot),])
  rmse_mid <- rvar_mean(((estimate_hier - true)^2)[seq(2,combo_pot,by=num_pot),])
  rmse_far <- rvar_mean(((estimate_hier - true)^2)[seq(3,combo_pot,by=num_pot),])
  
  df <- tibble(rmse=c(rmse,rmse_close,rmse_mid,rmse_far),
                     scale=c("all","close","mid","far"),
                     hier=TRUE)
  
  rmse <- rvar_mean((estimate_no_hier - true)^2)
  rmse_close <- rvar_mean(((estimate_no_hier - true)^2)[seq(1,combo_pot,by=num_pot),])
  rmse_mid <- rvar_mean(((estimate_no_hier - true)^2)[seq(2,combo_pot,by=num_pot),])
  rmse_far <- rvar_mean(((estimate_no_hier - true)^2)[seq(3,combo_pot,by=num_pot),])
  
  df <- df %>%
    bind_rows(tibble(rmse=c(rmse,rmse_close,rmse_mid,rmse_far),
               scale=c("all","close","mid","far"),
               hier=FALSE))
  
  df <- df %>%
    mutate(sim=sim)
  
}) %>%
  bind_rows()


### plot example SICs estimate vs ground truth
sic_tb <- grid %>%
  filter(sim %in% 1:6) %>%
  pmap(\(sim) {
    file_fit_hier <- paste0(path,"fit_sim",sim,"_hier.rds")
    file_fit_no_hier <- paste0(path,"fit_sim",sim,"_flat.rds")
    file_ground_truth <- paste0(path,"ground_truth_sim",sim,".rds")
    # print(file_fit)
    fit_hier <- readRDS(file_fit_hier)
    fit_no_hier <- readRDS(file_fit_no_hier)
    ground_truth <- readRDS(file_ground_truth)
    
    potentials <- make_rbfs(n_basis_functions = 3, max_dist = 75, basis_function_sigma = 15)
    
    # global
    beta_estimate_hier <- t(as_draws_rvars(fit_hier$draws(variables = "beta_local"))$beta_local)[-1,]
    
    beta_estimate_no_hier <- as_draws_rvars(fit_no_hier$draws(variables = "beta_local"))$beta_local[-1,]
    
    beta_true <- ground_truth$betas_local[-1,]
    # cbind(beta_true[,1],beta_estimate[,1])
    x_seq <- seq(0,100,1)
    x_des <- lapply(potentials,\(pot) pot(x_seq)) %>% do.call(cbind,.)
    ix <- 4:6
    b_t <- as.matrix(beta_true)[ix,]
    b_e_h <- as.matrix(beta_estimate_hier)[ix,]
    b_e_nh <- as.matrix(beta_estimate_no_hier)[ix,]
    
    
    lp_t <- x_des %*% b_t
    lp_e_h <- x_des %*% b_e_h
    lp_e_nh <- x_des %*% b_e_nh
    
    
    lp <- bind_cols(true=lp_t[,1],
                    estimate_hier=as.vector(lp_e_h[,1]),
                    estimate_no_hier=as.vector(lp_e_nh[,1]))

    # Compute simultaneous 95% credible bands using utility function
    bands <- compute_simultaneous_bands(
      lp_data = lp %>% select(estimate_hier, estimate_no_hier),
      x_seq = x_seq,
      alpha = 0.05
    )

    # Combine with true values and reshape for plotting
    lp <- bands %>%
      pivot_wider(
        id_cols = x,
        names_from = variable,
        values_from = c(mean, lower, upper)
      ) %>%
      mutate(True = as.vector(lp[[1]])) %>%
      rename(
        `Post. Mean (Hier)` = mean_estimate_hier,
        `Post. Mean (No Hier)` = mean_estimate_no_hier,
        `Lower Bound (Hier)` = lower_estimate_hier,
        `Upper Bound (Hier)` = upper_estimate_hier,
        `Lower Bound (No Hier)` = lower_estimate_no_hier,
        `Upper Bound (No Hier)` = upper_estimate_no_hier
      ) %>%
      pivot_longer(
        cols = c(`Post. Mean (Hier)`, `Post. Mean (No Hier)`, `True`),
        names_to = "Curve",
        values_to = "Value"
      ) %>%
      mutate(
        Lower = case_when(
          Curve == "Post. Mean (Hier)" ~ `Lower Bound (Hier)`,
          Curve == "Post. Mean (No Hier)" ~ `Lower Bound (No Hier)`,
          TRUE ~ NA_real_
        ),
        Upper = case_when(
          Curve == "Post. Mean (Hier)" ~ `Upper Bound (Hier)`,
          Curve == "Post. Mean (No Hier)" ~ `Upper Bound (No Hier)`,
          TRUE ~ NA_real_
        )
      ) %>%
      select(x, Curve, Value, Lower, Upper) 
    
    lp <- lp %>%
      mutate(sim=sim)
    
  }) %>% bind_rows()


# Save standardized analysis summary for fetching/plotting
analysis_summary <- list(
  rmse_tb = rmse_tb,
  sic_tb = sic_tb,
  metadata = list(
    sim_dir = "sim_flat_model",
    date_analyzed = Sys.time(),
    n_simulations = length(sim_indices)
  )
)

saveRDS(analysis_summary, paste0(path, "analysis_summary.rds"))
print(paste("Saved analysis summary to:", paste0(path, "analysis_summary.rds")))
