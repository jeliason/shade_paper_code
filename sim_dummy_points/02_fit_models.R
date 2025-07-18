library(tidyverse)
library(spatstat)
library(Matrix)
library(SHADE)

# set environment
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  path <- "./sim_dummy_points/data/"
  sim_idx <- 1
  sample_or_var <- "var"
} else {
  path <- "./sim_dummy_points/data/"
  args <- commandArgs(trailingOnly=TRUE)
  sim_idx <- as.numeric(args[1])
  sample_or_var <- "var"
  print(sim_idx)
}

n_cores <- 2
draws <- 1e3
grid <- expand.grid(ratio=c(0.5,1,2,5,10),
                    num_points_per_type=c(20,80,150,300,500),
                    sim=1:5)
ratio <- grid$ratio[sim_idx]
np <- grid$num_points_per_type[sim_idx]
n_dummy <- floor(np * ratio)
sim <- grid$sim[sim_idx]

file_data_stan <- paste0(path,"data_stan_sim",sim,"_ratio_",ratio,"_np_",np,".json")
file_fit <- paste0(path,"fit_sim",sim,"_ratio_",ratio,"_np_",np,".rds")

# Use run_SHADE_model from SHADE package
if(sample_or_var == "sample") {
  warmup <- 1e3
  fit <- run_SHADE_model(file_data_stan, 
                        method = "sample", 
                        iter_warmup = warmup, 
                        iter_sampling = draws, 
                        refresh = 5, 
                        chains = 1, 
                        threads_per_chain = n_cores)
} else if(sample_or_var == "var") {
  fit <- run_SHADE_model(file_data_stan, 
                        method = "variational", 
                        draws = draws, 
                        refresh = 5,
                        threads = n_cores)
}

fit$save_object(file_fit)

# cpp_options <- list(
#   "CXX" = "/usr/bin/clang++",
#   "CXXFLAGS+= -march=native")
# 
# install_cmdstan(overwrite = T,cpp_options = cpp_options)
