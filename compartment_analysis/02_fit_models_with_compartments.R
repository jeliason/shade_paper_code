# Fit SHADE Models with Binary Compartment Covariates
#
# Purpose: Fit SHADE models including tumor/stromal compartment as spatial covariate
# Approach: Use compartment spatial fields from 01_identify_compartments_binary.R

library(tidyverse)
library(spatstat)
library(Matrix)
library(SHADE)
source("utils.R")  # For MIN_INTERACTION_RADIUS

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("Loading data...\n")

# Load CRC data
df_raw <- read_csv("crc_analysis/data/CRC_cleaned.csv")
pt_data <- read_csv("crc_analysis/data/CRC_pt_metadata.csv")

# Recode cell types
df_raw <- df_raw %>%
  mutate(type = as.factor(type)) %>%
  mutate(type = fct_recode(
    type,
    "CAFs" = "smooth muscle",
    "hybrid E/M" = "stroma",
    "TAMs" = "CD163+ macros",
    "CTLs" = "CD8+ T cells"
  ))

spots <- pt_data %>% pull(Spot)

# Load binary compartment fields from step 01
compartment_fields <- readRDS("./compartment_analysis/data/binary_compartment_fields.rds")

cat(sprintf("Loaded compartment fields for %d images\n", length(compartment_fields)))

# ============================================================================
# 2. CONFIGURATION
# ============================================================================

path <- "./compartment_analysis/data/"
n_dummy <- 1e3

# Define target and source types (following crc_analysis convention)
targets <- c("CTLs", "memory CD4+ T", "granulocytes")
sources <- c("TAMs", "CAFs", "vasculature", "hybrid E/M", "tumor cells")

# Basis function configuration
n_basis <- 5  # Increased from 3 for more flexibility
max_dist <- 75
basis_sigma <- 12  # Slightly narrower for more local features

cat(sprintf("Configuration: %d basis functions, sigma=%d μm, truncated at %d μm\n",
            n_basis, basis_sigma, MIN_INTERACTION_RADIUS))

# Create truncated RBFs for plotting
truncated_rbfs <- make_truncated_rbfs(n_basis, max_dist, basis_sigma, MIN_INTERACTION_RADIUS)

# Evaluate over distance range
x_seq <- seq(0, max_dist, length.out = 200)
rbf_data <- map_dfr(seq_along(truncated_rbfs), function(i) {
  tibble(
    distance = x_seq,
    value = truncated_rbfs[[i]](x_seq),
    rbf = paste0("RBF ", i)
  )
})

# Plot
p <- ggplot(rbf_data, aes(x = distance, y = value, color = rbf)) +
  geom_vline(xintercept = MIN_INTERACTION_RADIUS, linetype = "dashed",
             color = "gray30", alpha = 0.5) +
  geom_line(linewidth = 1) +
  labs(
    x = "Distance (μm)",
    y = "Basis function value",
    color = "Basis"
  ) +
  theme_minimal()

print(p)

# ============================================================================
# 3. FIT MODELS WITH COMPARTMENT COVARIATE
# ============================================================================

for (type_idx in 1:length(targets)) {

  type <- targets[type_idx]
  types_keep <- c(type, sources)

  cat("\nProcessing:", type, "...")

  # Filter and create point patterns
  dats <- df_raw %>%
    filter(Spot %in% spots) %>%
    group_by(Spot) %>%
    group_map(~ {
      .x %>%
        filter(type %in% types_keep) %>%
        droplevels()
    })

  pats <- lapply(1:length(dats), function(i) {
    df <- dats[[i]]
    tryCatch({
      pat <- make_pat(df$X, df$Y, factor(df$type, levels = types_keep))
      sq_W <- owin(xrange = c(min(df$X), max(df$X)), yrange = c(min(df$Y), max(df$Y)))
      Window(pat) <- sq_W
      pat
    }, error = function(e) {
      print(i)
      stop(e)
    })
  })

  names(pats) <- spots

  # Create quadrature schemes
  Qs <- lapply(pats, make_quadrature, n_dummy = n_dummy)

  # Create compartment covariate for each image
  focal_type <- types_keep[1]

  covariate_list <- lapply(spots, function(spot) {
    Q <- Qs[[spot]]

    # Subset to focal cell type
    Q_focal <- quadscheme.logi(
      Q$data[marks(Q$data) == focal_type],
      Q$dummy[marks(Q$dummy) == focal_type]
    )

    # Drop unused levels
    marks(Q_focal$data)  <- droplevels(marks(Q_focal$data))
    marks(Q_focal$dummy) <- droplevels(marks(Q_focal$dummy))

    # Get counts
    n_data <- npoints(Q_focal$data)
    n_dummy <- npoints(Q_focal$dummy)
    n_total <- n_data + n_dummy

    # Get coordinates for data and dummy points
    coords_data <- coords(Q_focal$data)
    coords_dummy <- coords(Q_focal$dummy)
    all_x <- c(coords_data$x, coords_dummy$x)
    all_y <- c(coords_data$y, coords_dummy$y)

    # Debug: check dimensions
    if (length(all_x) != n_total) {
      cat(sprintf("  WARNING spot %s: n_total=%d but length(all_x)=%d\n",
                  spot, n_total, length(all_x)))
    }

    # Look up compartment at each quadrature point location
    # Check if this spot has a compartment field
    if (!(spot %in% names(compartment_fields))) {
      cat(sprintf("  Warning: no compartment field for spot %s, using 0\n", spot))
      compartment_values <- rep(0, length(all_x))
    } else {
      comp_im <- compartment_fields[[spot]]$compartment_im

      # Look up compartment values using safelookup (nearest-neighbor, no interpolation)
      # Create ppp with window encompassing ALL quadrature points (won't drop any)
      quad_window <- owin(xrange = range(all_x), yrange = range(all_y))
      quad_ppp <- ppp(all_x, all_y, window = quad_window)

      # safelookup returns NA for points outside comp_im, which we'll impute below
      compartment_values <- safelookup(comp_im, quad_ppp)

      # Check for NAs and impute
      n_na <- sum(is.na(compartment_values))
      if (n_na > 0) {
        cat(sprintf("  Note: %d/%d points outside compartment field, assigning to stromal (0)\n",
                    n_na, length(compartment_values)))
        # Impute NA values as stromal compartment (0)
        # Rationale: points outside tumor density field are by definition low tumor density
        compartment_values[is.na(compartment_values)] <- 0
      }
    }

    # Create covariate matrix
    # Column name: "compartment" (0 = stromal, 1 = tumor)
    cov_matrix <- cbind(
      compartment = compartment_values
    )

    # Debug: final check
    # cat(sprintf("  Spot %s: n_data=%d, n_dummy=%d, n_total=%d, cov_rows=%d\n",
    #             spot, n_data, n_dummy, n_total, nrow(cov_matrix)))

    return(cov_matrix)
  })

  names(covariate_list) <- spots

  # Gather coordinates
  coords <- do.call(rbind, lapply(seq_along(pats), function(i) {
    pat <- pats[[i]]
    tibble(
      x = pat$x,
      y = pat$y,
      type = pat$marks,
      image_id = spots[i]
    )
  }))

  # Create truncated RBFs that are zero below MIN_INTERACTION_RADIUS
  truncated_potentials <- make_truncated_rbfs(
    n_basis_functions = n_basis,
    max_dist = max_dist,
    basis_function_sigma = basis_sigma,
    min_dist = MIN_INTERACTION_RADIUS
  )

  prep <- prepare_spatial_model_data(
    x = coords$x,
    y = coords$y,
    cell_type = coords$type,
    image_id = factor(coords$image_id),
    patient_metadata = pt_data,
    type_idx = 1,
    n_dummy = n_dummy,
    potentials = truncated_potentials,
    covariate_list = covariate_list,
    quadrature_list = Qs
  )

  # Save prepared data
  file_json_comp <- paste0(path, "data_stan_", make.names(type), "_with_compartment.json")
  file_metadata_comp <- paste0(path, "metadata_", make.names(type), "_with_compartment.rds")

  write_json_chunked(prep$stan_data, file_json_comp, chunk_size = 1e6)
  saveRDS(prep$metadata, file_metadata_comp)

  # Fit model
  fit_comp <- run_SHADE_model(
    file_json_comp,
    method = "variational",
    draws = 1000,
    refresh = 100,
    threads = 2
  )

  # Save fit
  file_fit_comp <- paste0(path, "fit_", make.names(type), "_with_compartment.rds")
  fit_comp$save_object(file_fit_comp)

  cat(" done\n")
}

cat("\nAll models fitted with truncated basis functions!\n")
cat(sprintf("  %d basis functions, sigma=%d μm, truncated at %d μm\n",
            n_basis, basis_sigma, MIN_INTERACTION_RADIUS))
cat("Next: Run 03_compare_with_without.R\n")
