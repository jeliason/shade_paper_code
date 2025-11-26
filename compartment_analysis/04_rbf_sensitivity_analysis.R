# RBF Sensitivity Analysis for Compartment Models
#
# Purpose: Test sensitivity of SIC estimates to basis function choice
# Explores: different RBF configurations (n_basis, sigma, truncation) and B-splines

library(tidyverse)
library(spatstat)
library(Matrix)
library(SHADE)
library(posterior)
library(patchwork)
source("utils.R")

# ============================================================================
# 1. CONFIGURATION
# ============================================================================

cat("RBF Sensitivity Analysis\n")
cat("========================\n\n")

# Load data
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

# Load compartment fields
compartment_fields <- readRDS("./compartment_analysis/data/binary_compartment_fields.rds")

# Analysis configuration
path <- "./compartment_analysis/rbf_sensitivity/"
dir.create(path, recursive = TRUE, showWarnings = FALSE)

# For quick testing, just use one target type
# Change to targets <- c("CTLs", "memory CD4+ T", "granulocytes") for full analysis
targets <- c("CTLs")
sources <- c("TAMs", "CAFs", "vasculature", "hybrid E/M", "tumor cells")

n_dummy <- 1e3
max_dist <- 75

# ============================================================================
# 2. DEFINE BASIS FUNCTION CONFIGURATIONS
# ============================================================================

#' Create B-spline basis functions
#' @param n_basis_functions Number of B-spline basis functions
#' @param max_dist Maximum distance for support
#' @param degree Degree of B-spline (default 3 for cubic)
#' @param min_dist Minimum distance (default 0, or can truncate)
make_bspline_basis <- function(n_basis_functions,
                                max_dist,
                                degree = 3,
                                min_dist = 0) {

  # Create knot sequence
  # For B-splines, we need (n_basis + degree + 1) knots
  # We'll use equally spaced interior knots
  n_interior_knots <- n_basis_functions - degree - 1

  if (n_interior_knots < 0) {
    stop("Need at least degree + 1 basis functions for B-splines")
  }

  # Create knot vector with repeated endpoints for clamped B-splines
  interior_knots <- seq(min_dist, max_dist, length.out = n_interior_knots + 2)

  # Add repeated knots at boundaries
  knots <- c(
    rep(min_dist, degree),
    interior_knots,
    rep(max_dist, degree)
  )

  # Define B-spline basis functions
  bsplines <- lapply(1:n_basis_functions, function(i) {
    function(x) {
      # Use splines::bs() to evaluate the i-th B-spline
      # This is a simplified version - for production use splines package
      result <- numeric(length(x))

      # Only evaluate for x in valid range
      valid_idx <- x >= min_dist & x <= max_dist

      if (any(valid_idx)) {
        x_valid <- x[valid_idx]

        # Evaluate all B-splines at once
        basis_matrix <- splines::bs(
          x_valid,
          knots = interior_knots[2:(length(interior_knots)-1)],
          degree = degree,
          Boundary.knots = c(min_dist, max_dist)
        )

        # Extract the i-th basis function
        if (i <= ncol(basis_matrix)) {
          result[valid_idx] <- basis_matrix[, i]
        }
      }

      return(result)
    }
  })

  names(bsplines) <- paste0("bspline", 1:n_basis_functions)

  return(bsplines)
}

# Define configurations to test
configs <- list(
  # Standard RBFs (no truncation)
  list(
    name = "RBF_3_std",
    label = "RBF: 3 basis, σ=15 μm",
    type = "rbf",
    n_basis = 3,
    sigma = 15,
    truncate = FALSE
  ),
  list(
    name = "RBF_5_std",
    label = "RBF: 5 basis, σ=12 μm",
    type = "rbf",
    n_basis = 5,
    sigma = 12,
    truncate = FALSE
  ),
  list(
    name = "RBF_7_std",
    label = "RBF: 7 basis, σ=10 μm",
    type = "rbf",
    n_basis = 7,
    sigma = 10,
    truncate = FALSE
  ),

  # Truncated RBFs
  list(
    name = "RBF_3_trunc",
    label = "RBF: 3 basis, σ=15 μm, truncated",
    type = "rbf",
    n_basis = 3,
    sigma = 15,
    truncate = TRUE
  ),
  list(
    name = "RBF_5_trunc",
    label = "RBF: 5 basis, σ=12 μm, truncated",
    type = "rbf",
    n_basis = 5,
    sigma = 12,
    truncate = TRUE
  ),
  list(
    name = "RBF_7_trunc",
    label = "RBF: 7 basis, σ=10 μm, truncated",
    type = "rbf",
    n_basis = 7,
    sigma = 10,
    truncate = TRUE
  ),

  # B-splines
  list(
    name = "Bspline_5",
    label = "B-spline: 5 basis, cubic",
    type = "bspline",
    n_basis = 5,
    degree = 3,
    truncate = FALSE
  ),
  list(
    name = "Bspline_7",
    label = "B-spline: 7 basis, cubic",  # Will be filtered in comparison plots
    type = "bspline",
    n_basis = 7,
    degree = 3,
    truncate = FALSE
  ),
  list(
    name = "Bspline_5_trunc",
    label = "B-spline: 5 basis, cubic, truncated",
    type = "bspline",
    n_basis = 5,
    degree = 3,
    truncate = TRUE
  )
)

cat("Testing", length(configs), "basis function configurations:\n")
for (cfg in configs) {
  cat("  -", cfg$label, "\n")
}
cat("\n")

# ============================================================================
# 3. FIT MODELS WITH EACH CONFIGURATION
# ============================================================================

cat("Fitting models with each configuration...\n")
cat("==========================================\n\n")

all_fits <- list()
all_metadata <- list()

for (config in configs) {

  cfg_name <- config$name
  cat("\nConfiguration:", config$label, "\n")
  cat(strrep("-", 60), "\n")

  # Create basis functions
  if (config$type == "rbf") {
    if (config$truncate) {
      potentials <- make_truncated_rbfs(
        config$n_basis,
        max_dist,
        config$sigma,
        MIN_INTERACTION_RADIUS
      )
    } else {
      potentials <- make_rbfs(
        max_dist,
        config$n_basis,
        config$sigma
      )
    }
  } else if (config$type == "bspline") {
    min_dist <- if (config$truncate) MIN_INTERACTION_RADIUS else 0
    potentials <- make_bspline_basis(
      config$n_basis,
      max_dist,
      config$degree,
      min_dist
    )
  }

  # Visualize basis functions
  x_seq_viz <- seq(0, max_dist, length.out = 200)
  basis_data <- map_dfr(seq_along(potentials), function(i) {
    tibble(
      distance = x_seq_viz,
      value = potentials[[i]](x_seq_viz),
      basis = paste0("Basis ", i)
    )
  })

  p_basis <- ggplot(basis_data, aes(x = distance, y = value, color = basis)) +
    geom_vline(xintercept = MIN_INTERACTION_RADIUS, linetype = "dashed",
               color = "gray30", alpha = 0.5) +
    geom_line(linewidth = 1) +
    labs(
      title = config$label,
      x = "Distance (μm)",
      y = "Basis function value",
      color = "Basis"
    ) +
    theme_minimal()

  ggsave(
    paste0(path, "basis_", cfg_name, ".pdf"),
    p_basis,
    width = 8,
    height = 5
  )

  # Fit models for each target type
  for (type in targets) {

    types_keep <- c(type, sources)

    cat("  Processing:", type, "...\n")

    # Check if model already exists
    file_fit <- paste0(path, "fit_", make.names(type), "_", cfg_name, ".rds")
    file_metadata <- paste0(path, "metadata_", make.names(type), "_", cfg_name, ".rds")

    if (file.exists(file_fit) && file.exists(file_metadata)) {
      cat("    Model already exists, loading from disk...\n")
      fit <- readRDS(file_fit)
      metadata <- readRDS(file_metadata)

      # Store for later extraction
      all_fits[[paste0(cfg_name, "_", type)]] <- fit
      all_metadata[[paste0(cfg_name, "_", type)]] <- metadata

      next
    }

    cat("    Fitting new model...\n")

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

      # Get coordinates for data and dummy points
      coords_data <- coords(Q_focal$data)
      coords_dummy <- coords(Q_focal$dummy)
      all_x <- c(coords_data$x, coords_dummy$x)
      all_y <- c(coords_data$y, coords_dummy$y)

      # Look up compartment at each quadrature point location
      if (!(spot %in% names(compartment_fields))) {
        compartment_values <- rep(0, length(all_x))
      } else {
        comp_im <- compartment_fields[[spot]]$compartment_im

        quad_window <- owin(xrange = range(all_x), yrange = range(all_y))
        quad_ppp <- ppp(all_x, all_y, window = quad_window)

        compartment_values <- safelookup(comp_im, quad_ppp)

        # Impute NA values as stromal (0)
        if (any(is.na(compartment_values))) {
          compartment_values[is.na(compartment_values)] <- 0
        }
      }

      # Create covariate matrix
      cov_matrix <- cbind(compartment = compartment_values)

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

    # Prepare data for SHADE
    prep <- prepare_spatial_model_data(
      x = coords$x,
      y = coords$y,
      cell_type = coords$type,
      image_id = factor(coords$image_id),
      patient_metadata = pt_data,
      type_idx = 1,
      n_dummy = n_dummy,
      potentials = potentials,
      covariate_list = covariate_list,
      quadrature_list = Qs
    )

    # Save prepared data
    file_json <- paste0(path, "data_stan_", make.names(type), "_", cfg_name, ".json")
    file_metadata <- paste0(path, "metadata_", make.names(type), "_", cfg_name, ".rds")

    write_json_chunked(prep$stan_data, file_json, chunk_size = 1e6)
    saveRDS(prep$metadata, file_metadata)

    # Fit model
    fit <- run_SHADE_model(
      file_json,
      method = "variational",
      draws = 1000,
      refresh = 500,
      threads = 2
    )

    # Save fit (file paths already defined above)
    fit$save_object(file_fit)

    # Store for later extraction
    all_fits[[paste0(cfg_name, "_", type)]] <- fit
    all_metadata[[paste0(cfg_name, "_", type)]] <- prep$metadata

    cat("    Done.\n")
  }
}

cat("\n\nAll models fitted!\n")

# ============================================================================
# 4. EXTRACT COHORT SICs FROM ALL CONFIGURATIONS
# ============================================================================

cat("\n\nExtracting cohort SICs from all configurations...\n")
cat("==================================================\n\n")

#' Extract cohort-level SIC curves
extract_cohort_sics_sensitivity <- function(fit, metadata, source_types, target_type,
                                            pt_data, potentials, config_name,
                                            alpha = 0.1) {

  # Extract draws
  draws <- posterior::as_draws_rvars(fit$draws())
  beta_global <- draws$beta_global
  coef_names <- metadata$coef_names
  rownames(beta_global) <- coef_names

  # Get groups
  pt_df <- pt_data %>% filter(Spot %in% metadata$spots)
  groups <- levels(factor(pt_df$Group))

  # Create design matrix
  # Only for distances where potentials are non-zero
  x_seq <- seq(MIN_INTERACTION_RADIUS, 100, 1)
  x_des <- lapply(potentials, \(pot) pot(x_seq)) %>% do.call(cbind, .)

  # Extract SICs for each source type
  sics_list <- list()

  for (source in source_types) {
    # Find coefficients
    ix <- grep(paste0("_", target_type, "_", source), rownames(beta_global), fixed = TRUE)

    if (length(ix) == 0) {
      warning(paste0("No coefficients found for ", target_type, " -> ", source))
      next
    }

    # Extract coefficients
    b_g <- as.matrix(beta_global)[ix, ]

    # Compute linear predictor
    lp_g <- x_des %*% b_g

    # Extract cohort-level columns
    num_groups <- length(groups)
    cohort_cols <- (ncol(lp_g) - num_groups + 1):ncol(lp_g)
    lp_cohort <- lp_g[, cohort_cols, drop = FALSE]
    colnames(lp_cohort) <- groups

    # Compute simultaneous credible bands
    bands <- compute_simultaneous_bands(
      lp_data = as.data.frame(lp_cohort),
      x_seq = x_seq,
      alpha = alpha
    )

    # Format output
    bands_formatted <- bands %>%
      rename(
        distance = x,
        group = variable,
        mean = mean,
        lower = lower,
        upper = upper
      ) %>%
      mutate(
        source = source,
        config = config_name
      )

    sics_list[[source]] <- bands_formatted
  }

  bind_rows(sics_list)
}

# Extract SICs for all configurations
all_sics <- list()

for (config in configs) {
  cfg_name <- config$name

  # Recreate potentials
  if (config$type == "rbf") {
    if (config$truncate) {
      potentials <- make_truncated_rbfs(config$n_basis, max_dist, config$sigma, MIN_INTERACTION_RADIUS)
    } else {
      potentials <- make_rbfs(max_dist, config$n_basis, config$sigma)
    }
  } else if (config$type == "bspline") {
    min_dist <- if (config$truncate) MIN_INTERACTION_RADIUS else 0
    potentials <- make_bspline_basis(config$n_basis, max_dist, config$degree, min_dist)
  }

  for (type in targets) {
    cat("  Extracting:", config$label, "-", type, "\n")

    # Load fit and metadata
    file_fit <- paste0(path, "fit_", make.names(type), "_", cfg_name, ".rds")
    file_metadata <- paste0(path, "metadata_", make.names(type), "_", cfg_name, ".rds")

    fit <- readRDS(file_fit)
    metadata <- readRDS(file_metadata)

    sics <- extract_cohort_sics_sensitivity(
      fit, metadata, sources, type, pt_data, potentials,
      config$label,
      alpha = 0.1
    ) %>%
      mutate(target = type)

    all_sics[[paste0(cfg_name, "_", type)]] <- sics
  }
}

sics_data <- bind_rows(all_sics)

# Save
saveRDS(sics_data, paste0(path, "all_sics_sensitivity.rds"))
cat("\nSaved SIC data to:", paste0(path, "all_sics_sensitivity.rds\n"))

# ============================================================================
# 5. CREATE COMPARISON VISUALIZATIONS
# ============================================================================

cat("\n\nCreating comparison visualizations...\n")
cat("=====================================\n\n")

# For each configuration, create CLR vs DII facet_grid plot
for (config in configs) {
  cfg_name <- config$name
  cfg_label <- config$label

  cat("  Creating plot for:", cfg_label, "\n")

  plot_data <- sics_data %>%
    filter(config == cfg_label) %>%
    filter(distance >= MIN_INTERACTION_RADIUS)

  p <- ggplot(plot_data, aes(x = distance, y = mean, color = group, fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(
      values = c(CLR = "#D95F02", DII = "#1B9E77"),
      labels = c(CLR = "CLR", DII = "DII")
    ) +
    scale_fill_manual(
      values = c(CLR = "#D95F02", DII = "#1B9E77"),
      labels = c(CLR = "CLR", DII = "DII")
    ) +
    facet_grid(target ~ source, scales = "free_y") +
    labs(
      title = paste("CLR vs DII:", cfg_label),
      subtitle = "Compartment-adjusted models (90% simultaneous CI)",
      x = "Distance (μm)",
      y = "SIC",
      color = "Group",
      fill = "Group"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 9),
      plot.title = element_text(size = 12, face = "bold")
    )

  ggsave(
    paste0(path, "sics_", cfg_name, ".pdf"),
    p,
    width = 14,
    height = 8
  )
}

# Create mega comparison plot: one interaction across all configurations
cat("\n  Creating configuration comparison plots...\n")

sics_data %>%
  filter(str_detect(config,"B-spline: 5 basis"))

  fit <- readRDS("compartment_analysis/rbf_sensitivity/fit_CTLs_Bspline_7.rds")
  draws <- posterior::as_draws_rvars(fit$draws())

  # Check if beta_global has variation
  beta_global <- draws$beta_global
  sd(beta_global[1,1])  # Should be non-zero

# Pick a few representative interactions
example_interactions <- list(
  list(target = "CTLs", source = "TAMs", group = "CLR"),
  list(target = "CTLs", source = "tumor cells", group = "DII"),
  list(target = "CTLs", source = "CAFs", group = "CLR")
)

for (interaction in example_interactions) {

  plot_data <- sics_data %>%
    filter(
      target == interaction$target,
      source == interaction$source,
      group == interaction$group
    )%>%
    # Filter out 7 basis configurations (unstable)
    filter(!grepl("5 basis", config))

  # Debug: print which configs are being plotted
  cat(sprintf("    %s → %s (%s): %d configs - %s\n",
              interaction$source,
              interaction$target,
              interaction$group,
              length(unique(plot_data$config)),
              paste(unique(plot_data$config), collapse = ", ")))

  p <- ggplot(plot_data, aes(x = distance, y = mean, color = config, fill = config)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    labs(
      title = paste0(interaction$source, " → ", interaction$target, " (", interaction$group, ")"),
      subtitle = "Sensitivity to basis function choice (90% simultaneous CI, 7-basis excluded)",
      x = "Distance (μm)",
      y = "SIC",
      color = "Configuration",
      fill = "Configuration"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 12, face = "bold")
    )

  filename <- paste0(
    "comparison_",
    make.names(interaction$source), "_",
    make.names(interaction$target), "_",
    interaction$group,
    ".pdf"
  )

  ggsave(
    paste0(path, filename),
    p,
    width = 10,
    height = 6
  )
}

# ============================================================================
# 6. SUMMARY STATISTICS
# ============================================================================

cat("\n\nComputing summary statistics...\n")

# Compare SIC estimates across configurations
summary_stats <- sics_data %>%
  group_by(target, source, group, distance) %>%
  summarize(
    n_configs = n(),
    mean_sic = mean(mean),
    sd_sic = sd(mean),
    range_sic = max(mean) - min(mean),
    .groups = "drop"
  ) %>%
  group_by(target, source, group) %>%
  summarize(
    avg_sd = mean(sd_sic),
    max_sd = max(sd_sic),
    avg_range = mean(range_sic),
    max_range = max(range_sic),
    .groups = "drop"
  ) %>%
  arrange(desc(max_range))

cat("\nVariability across configurations:\n")
print(summary_stats, n = 20)

write_csv(summary_stats, paste0(path, "sensitivity_summary.csv"))

# ============================================================================
# 7. SUMMARY
# ============================================================================

cat("\n\n========================================\n")
cat("Sensitivity analysis complete!\n")
cat("========================================\n\n")

cat("Tested", length(configs), "basis function configurations\n")
cat("Fitted", length(targets), "target cell type(s)\n")
cat("\nOutput files in:", path, "\n")
cat("  - Basis function plots: basis_*.pdf\n")
cat("  - SIC plots by config:  sics_*.pdf\n")
cat("  - Comparison plots:     comparison_*.pdf\n")
cat("  - Summary statistics:   sensitivity_summary.csv\n")
cat("  - All SIC data:         all_sics_sensitivity.rds\n")
cat("\n")
