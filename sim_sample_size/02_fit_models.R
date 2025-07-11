library(tidyverse)
library(spatstat)
library(Matrix)
library(SHADE)

# set environment
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  path <- "./sim_control_confounding/data/"
  sim_idx <- 7
  sample_or_var <- "var"
} else {
  path <- "./sim_control_confounding/data/"
  args <- commandArgs(trailingOnly=TRUE)
  sim_idx <- as.numeric(args[1])
  sample_or_var <- "var"
  print(sim_idx)
}

n_cores <- 2
draws <- 1e3

sim <- sim_idx

file_data_stan <- paste0(path,"data_stan_sim",sim,".json")
file_fit <- paste0(path,"fit_sim",sim,".rds")

# Use run_SHADE_model from SHADE package
if(sample_or_var == "sample") {
  warmup <- 1e3
  fit <- run_SHADE_model(file_data_stan, 
                        method = "sample", 
                        iter_warmup = 500, 
                        iter_sampling = 500, 
                        refresh = 20, 
                        chains = 1, 
                        threads_per_chain = 4)
} else if(sample_or_var == "var") {
  fit <- run_SHADE_model(file_data_stan, 
                        method = "variational", 
                        draws = draws, 
                        refresh = 5,
                        threads = n_cores)
}

fit$save_object(file_fit)

library(posterior)
draws <- as_draws_rvars(fit$draws())

draws$beta_global
params$betas_global
# cpp_options <- list(
#   "CXX" = "/usr/bin/clang++",
#   "CXXFLAGS+= -march=native")
# 
# install_cmdstan(overwrite = T,cpp_options = cpp_options)
