library(tidyverse)
library(spatstat)
library(Matrix)
library(SHADE)

# Load utilities
source("utils.R")

# Set environment and get data path
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  sim_idx <- 2
  sample_or_var <- "var"
} else {
  args <- commandArgs(trailingOnly=TRUE)
  sim_idx <- as.numeric(args[1])
  sample_or_var <- "var"
  print(paste("Simulation index:", sim_idx))
}

# Get data path (handles HPC vs local automatically)
path <- get_data_path("sim_flat_model")
print(paste("Data path:", path))

n_cores <- 2
draws <- 1e3
grid <- expand.grid(sim=1:100)
sim <- grid$sim[sim_idx]

file_data_stan <- paste0(path,"data_stan_sim",sim,".json")
file_fit_hier <- paste0(path,"fit_sim",sim,"_hier.rds")
file_fit_flat <- paste0(path,"fit_sim",sim,"_flat.rds")


# Use run_SHADE_model from SHADE package for hierarchical model
fit_hier <- run_SHADE_model(file_data_stan, 
                           method = "variational", 
                           draws = draws, 
                           refresh = 5,
                           threads = n_cores)
fit_hier$save_object(file_fit_hier)

mod_flat <- cmdstanr::cmdstan_model("sim_flat_model/SHADE_flat.stan",cpp_options = list(stan_threads = TRUE))
fit_flat <- mod_flat$variational(data = file_data_stan, draws = draws, refresh = 5,threads=n_cores)
fit_flat$save_object(file_fit_flat)

