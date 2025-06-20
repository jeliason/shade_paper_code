library(tidyverse)
library(spatstat)
library(Matrix)
library(SHADE)


path <- "./crc_analysis/data/"
n_cores <- 2
type_idx <- 1

targets <- c(
  "CTLs",
  "memory CD4+ T",
  "granulocytes"
)

sources <- c(
  "TAMs",
  "CAFs",
  "vasculature",
  "hybrid E/M",
  "tumor cells"
)

type <- targets[type_idx]

sample_or_var <- "var"
draws <- 1e3

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



