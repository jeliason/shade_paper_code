library(tidyverse)
library(spatstat)
library(Matrix)
library(SHADE)

# CRC dataset needs to be loaded here
# Example:
# df_raw <- read_csv("path/to/CRC_data.csv")
cat("NOTE: CRC dataset needs to be loaded here. Replace this with the actual data loading code.\n")

df_raw %>%
  count(Spot,type) %>%
  group_by(type) %>%
  summarise(mn = median(n)) %>%
  filter(mn > 30) %>%
  pull(type) %>%
  as.character() -> types_keep

# set environment
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV == "laptop") {
  path <- "./crc_analysis/data/"
  n_cores <- 2
  type_idx <- 1
} else {
  path <- "./crc_analysis/data/"
  n_cores <- 16
  args <- commandArgs(trailingOnly=TRUE)
  type_idx <- as.numeric(args[1])
}

sample_or_var <- "var"
draws <- 1e3
type <- types_keep[type_idx]

# Load the model with run_SHADE_model from SHADE package
print("model loaded")

file_json <- paste0(path,"data_stan_CRC_type_",make.names(type),".json")
file_fit <- paste0(path,"fit_CRC_type_",make.names(type),".rds")


if(sample_or_var == "sample") {
  warmup <- 1e3
  fit <- run_SHADE_model(file_json, 
                        method = "sample", 
                        iter_warmup = warmup, 
                        iter_sampling = draws, 
                        refresh = 5, 
                        chains = 1, 
                        threads_per_chain = n_cores)
} else if(sample_or_var == "var") {
  print("variational fitting")
  fit <- run_SHADE_model(file_json, 
                        method = "variational", 
                        draws = draws, 
                        refresh = 5,
                        threads = n_cores)
} else if(sample_or_var == "pf") {
  print("pathfinder")
  fit <- run_SHADE_model(file_json,
                        method = "pathfinder",
                        num_threads = n_cores, 
                        single_path_draws = draws,
                        num_paths = 4,
                        show_messages = TRUE)
}

fit$save_object(file_fit)


