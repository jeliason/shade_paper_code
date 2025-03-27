library(tidyverse)
library(spatstat)
library(Matrix)

# set environment
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV == "laptop") {
  path <- "./no_hier/"
  sim_idx <- 3
  sample_or_var <- "var"
} else {
  path <- "/nfs/turbo/umms-ukarvind/joelne/BISTRO/no_hier/"
  args <- commandArgs(trailingOnly=TRUE)
  sim_idx <- as.numeric(args[1])
  sample_or_var <- "var"
  print(sim_idx)
}

n_cores <- 2
draws <- 1e3
grid <- expand.grid(sim=1:100)
sim <- grid$sim[sim_idx]

file_data_stan <- paste0(path,"data_stan_sim",sim,".json")
file_fit_hier <- paste0(path,"fit_sim",sim,"_hier.rds")
file_fit_no_hier <- paste0(path,"fit_sim",sim,"_no_hier.rds")


mod <- cmdstanr::cmdstan_model("scripts/hier_asymm.stan",cpp_options = list(stan_threads = TRUE))
fit_hier <- mod$variational(data = file_data_stan, draws = draws, refresh = 5,threads=n_cores)
fit_hier$save_object(file_fit_hier)

mod <- cmdstanr::cmdstan_model("scripts/no_hier_asymm.stan",cpp_options = list(stan_threads = TRUE))
fit_no_hier <- mod$variational(data = file_data_stan, draws = draws, refresh = 5,threads=n_cores)
fit_no_hier$save_object(file_fit_no_hier)
