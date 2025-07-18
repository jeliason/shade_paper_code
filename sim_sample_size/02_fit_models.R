library(tidyverse)
library(spatstat)
library(Matrix)
library(SHADE)

# set environment
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  path <- "./sim_sample_size/data/"
  sim_idx <- 50
  sample_or_var <- "var"
} else {
  path <- "./sim_sample_size/data/"
  args <- commandArgs(trailingOnly=TRUE)
  sim_idx <- as.numeric(args[1])
  sample_or_var <- "var"
  print(sim_idx)
}

n_cores <- 2
draws <- 1e3
# parameters to adjust
grid <- expand.grid(num_pts_per_group=c(10,20,40),
                    images_per_pt=c(1,2,4),
                    sim=1:15)
num_pts_per_group <- grid$num_pts_per_group[sim_idx]
images_per_pt <- grid$images_per_pt[sim_idx]
sim <- grid$sim[sim_idx]

file_data_stan <- paste0(path,"data_stan_sim",sim,"_pts_",num_pts_per_group,"_ims_",images_per_pt,".json")
file_fit <- paste0(path,"fit_sim",sim,"_pts_",num_pts_per_group,"_ims_",images_per_pt,".rds")

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

# cpp_options <- list(
#   "CXX" = "/usr/bin/clang++",
#   "CXXFLAGS+= -march=native")
# 
# install_cmdstan(overwrite = T,cpp_options = cpp_options)
