library(tidyverse)
library(spatstat)
library(cmdstanr)
library(Matrix)
library(posterior)
library(ggdist)
library(survival)
library(SHADE)

# set environment
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  path <- "./sim_flat_model/data/"
  sim_indices <- c(1,2,3)
} else {
  path <- "./sim_flat_model/data/"
  sim_indices <- 1:100
}
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
    
    lp <- lp %>%
      mutate(x=x_seq) %>%
      mutate(across(c(estimate_hier,estimate_no_hier),list(
        mn = ~as.vector(E(.)),
        lo = ~as.vector(quantile(.,0.1)),
        hi = ~as.vector(quantile(.,0.9))
      ),.names="{.col}_{.fn}"),.keep="unused") %>%
      rename(
        `Post. Mean (Hier)` = estimate_hier_mn,
        `Post. Mean (No Hier)` = estimate_no_hier_mn,
        `Lower Bound (Hier)` = estimate_hier_lo,
        `Upper Bound (Hier)` = estimate_hier_hi,
        `Lower Bound (No Hier)` = estimate_no_hier_lo,
        `Upper Bound (No Hier)` = estimate_no_hier_hi,
        `True` = true
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


out <- list(
  rmse_tb=rmse_tb,
  sic_tb=sic_tb
)

saveRDS(out,paste0(path,"analysis_results.rds"))
