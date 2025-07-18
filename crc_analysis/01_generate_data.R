library(tidyverse)
library(spatstat)
library(Matrix)
library(cmdstanr)
library(SHADE)

source("crc_analysis/utils.R")

# CRC dataset needs to be loaded here
# Example:
df_raw <- read_csv("crc_analysis/data/CRC_cleaned.csv")
pt_data <- read_csv("crc_analysis/data/CRC_pt_metadata.csv")

df_raw <- df_raw %>%
  dplyr::mutate(type = as.factor(type)) %>%
  mutate(type = fct_recode(type,"CAFs"="smooth muscle","hybrid E/M"="stroma","TAMs"="CD163+ macros","CTLs"="CD8+ T cells"))

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

spots <- pt_data %>%
  pull(Spot)
  
path <- "./crc_analysis/data/"
n_dummy <- 1e3

type <- targets[type_idx]
types_keep <- c(type,sources)

num_types <- length(types_keep)

num_pot <- 3
potentials <- make_rbfs(n_basis_functions = num_pot, max_dist = 75, basis_function_sigma = 15)
plot_rbfs(seq(0,80,1),potentials)

dats <- df_raw %>%
  dplyr::filter(Spot %in% spots) %>%
  group_by(Spot) %>%
  group_map(~{
    .x <- .x %>%
      filter(type %in% types_keep) %>%
      droplevels()
  })

pats <- lapply(1:length(dats),\(i) {
  df <- dats[[i]]
  tryCatch({
    pat <- make_pat(df$X,df$Y,factor(df$type,levels=types_keep))
    sq_W <- owin(xrange = c(min(df$X),max(df$X)),yrange=c(min(df$Y),max(df$Y)))
    Window(pat) <- sq_W
    pat
  },error=\(e) {
    print(i)
    stop(e)
  })
})

file_json <- paste0(path,"data_stan_CRC_type_",make.names(type),".json")
file_dlist <- paste0(path,"data_lists_CRC_type_",make.names(type),".rds")
file_metadata <- paste0(path,"metadata_CRC_type_",make.names(type),".rds")

# Gather coordinates from all pattern objects
coords <- do.call(rbind, lapply(seq_along(pats), function(i) {
  pat <- pats[[i]]
  tibble(
    x = pat$x,
    y = pat$y,
    type = pat$marks,
    image_id = spots[i]
  )
}))

# Use prepare_spatial_model_data from SHADE package
prep <- prepare_spatial_model_data(
  x = coords$x,
  y = coords$y,
  cell_type = coords$type,
  image_id = factor(coords$image_id),
  patient_metadata = pt_data,
  type_idx = which(types_keep == type),
  n_dummy = n_dummy,
  n_basis_functions = num_pot
)

stan_data <- prep$stan_data

print("writing json...")
write_json_chunked(stan_data, file_json,chunk_size = 1e6)
print("json written.")


saveRDS(prep$metadata,file_metadata)
