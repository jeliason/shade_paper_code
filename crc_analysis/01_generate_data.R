library(tidyverse)
library(spatstat)
library(Matrix)
library(cmdstanr)
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

spots_pt <- pt_data

spots <- spots_pt %>%
  pull(Spot)
  
# set environment
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV == "laptop") {
  path <- "./crc_analysis/data/"
  # types_keep <- types_keep[1:2]
  # scale_factor <- 1
  n_dummy <- 1e3
  spots <- spots[c(1,2,5,6)]
  
} else {
  path <- "./crc_analysis/data/"
  # scale_factor <- 2
  n_dummy <- 1e3
  args <- commandArgs(trailingOnly=TRUE)
  type_idx <- as.numeric(args[1])
}

spots_pt <- spots_pt %>%
  filter(Spot %in% spots)

num_types <- length(types_keep)

num_pot <- 3
potentials <- make_rbfs(n_basis_functions = num_pot, max_dist = 75, basis_function_sigma = 15)
grainsize <- 1
mean_alpha <- -10
scale_sigmas <- 5

sample_to_indiv <- spots_pt %>%
  pull(Patient) %>%
  as.factor() %>%
  as.numeric()
indiv_to_group <- spots_pt %>%
  distinct(Patient,Group) %>%
  pull(Group) %>%
  as.factor() %>%
  as.numeric()

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

type <- types_keep[type_idx]
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

# Create patient metadata
patient_metadata <- spots_pt %>%
  select(Spot, Patient, Group)

# Use prepare_spatial_model_data from SHADE package
prep <- prepare_spatial_model_data(
  x = coords$x,
  y = coords$y,
  cell_type = coords$type,
  image_id = coords$image_id,
  patient_metadata = patient_metadata,
  type_idx = which(types_keep == type),
  n_dummy = n_dummy,
  n_basis_functions = num_pot,
  potentials = potentials
)

data_stan <- prep$data_stan

print("writing json...")
write_json_chunked(data_stan, file_json,chunk_size = 1e6)
print("json written.")

metadata <- list(
  coef_names = colnames(x_cells),
  potentials = potentials,
  spots = unique(pt_data$Spot),
  types = types_keep,
  sample_to_indiv = sample_to_indiv,
  indiv_to_group = indiv_to_group
)

saveRDS(metadata,file_metadata)
