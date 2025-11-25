# Analyze Compartment-Adjusted Models
#
# Purpose: Extract and visualize SICs from compartment-adjusted SHADE models
# Focus: CLR vs DII group comparisons with compartment covariates

library(tidyverse)
library(posterior)
library(patchwork)
source("utils.R")  # For compute_simultaneous_bands(), compute_pointwise_bands()

# ============================================================================
# 1. CONFIGURATION
# ============================================================================

# Load patient metadata
pt_data <- read_csv("crc_analysis/data/CRC_pt_metadata.csv")

targets <- c("CTLs", "memory CD4+ T", "granulocytes")
sources <- c("TAMs", "CAFs", "vasculature", "hybrid E/M", "tumor cells")

path_compartment <- "./compartment_analysis/data/"
path_figures <- "./compartment_analysis/compartment_analysis_figures/"

# Create figures directory if needed
dir.create(path_figures, recursive = TRUE, showWarnings = FALSE)

# Define truncated RBF basis functions
n_basis <- 5
max_dist <- 75
basis_sigma <- 12

# Create truncated RBFs for analysis (using function from utils.R)
truncated_potentials <- make_truncated_rbfs(n_basis, max_dist, basis_sigma, MIN_INTERACTION_RADIUS)

cat("Using truncated RBFs:", n_basis, "basis functions, sigma =", basis_sigma, "μm\n")
cat("  Truncated at", MIN_INTERACTION_RADIUS, "μm (zero below this distance)\n\n")

# ============================================================================
# 2. EXTRACT COHORT-LEVEL SICs
# ============================================================================

#' Extract cohort-level SIC curves from a fitted SHADE model
extract_cohort_sics <- function(fit, metadata, source_types, target_type, pt_data,
                                potentials,
                                band_type = "simultaneous", alpha = 0.1) {

  # Extract draws using rvar format
  draws <- posterior::as_draws_rvars(fit$draws())

  beta_global <- draws$beta_global
  coef_names <- metadata$coef_names

  # Set rownames for coefficient extraction
  rownames(beta_global) <- coef_names

  # Get patient metadata for group information
  pt_df <- pt_data %>% filter(Spot %in% metadata$spots)
  groups <- levels(factor(pt_df$Group))

  # Create design matrix from potentials
  # Only compute for distances >= MIN_INTERACTION_RADIUS (truncated RBFs are zero below this)
  x_seq <- seq(MIN_INTERACTION_RADIUS, 100, 1)
  x_des <- lapply(potentials, \(pot) pot(x_seq)) %>% do.call(cbind, .)

  # Extract SICs for each source type
  sics_list <- list()

  for (source in source_types) {
    # Find coefficients for this target-source pair
    ix <- grep(paste0("_", target_type, "_", source), rownames(beta_global), fixed = TRUE)

    if (length(ix) == 0) {
      warning(paste0("No coefficients found for ", target_type, " -> ", source))
      next
    }

    # Extract coefficients
    b_g <- as.matrix(beta_global)[ix, ]

    # Compute linear predictor via matrix multiplication
    lp_g <- x_des %*% b_g

    # Extract cohort-level columns (last num_groups columns)
    num_groups <- length(groups)
    cohort_cols <- (ncol(lp_g) - num_groups + 1):ncol(lp_g)
    lp_cohort <- lp_g[, cohort_cols, drop = FALSE]

    # Set column names to group names
    colnames(lp_cohort) <- groups

    # Compute credible bands
    if (band_type == "simultaneous") {
      bands <- compute_simultaneous_bands(
        lp_data = as.data.frame(lp_cohort),
        x_seq = x_seq,
        alpha = alpha
      )
    } else if (band_type == "pointwise") {
      bands <- compute_pointwise_bands(
        lp_data = as.data.frame(lp_cohort),
        x_seq = x_seq,
        alpha = alpha
      )
    } else {
      stop("band_type must be 'simultaneous' or 'pointwise'")
    }

    # Format output
    bands_formatted <- bands %>%
      rename(
        distance = x,
        group = variable,
        mean = mean,
        lower = lower,
        upper = upper
      ) %>%
      mutate(source = source, band_type = band_type)

    sics_list[[source]] <- bands_formatted
  }

  bind_rows(sics_list)
}

# ============================================================================
# 3. LOAD ALL COMPARTMENT MODELS AND EXTRACT SICS
# ============================================================================

cat("Extracting SICs from compartment-adjusted models...\n")
cat("===================================================\n\n")

# Extract with simultaneous bands
all_sics_simul <- list()

for (type in targets) {
  cat("Processing:", type, "(simultaneous bands)\n")

  file_fit <- paste0(path_compartment, "fit_", make.names(type), "_with_compartment.rds")
  file_metadata <- paste0(path_compartment, "metadata_", make.names(type), "_with_compartment.rds")

  if (!file.exists(file_fit)) {
    cat("  Warning: Fit not found:", file_fit, "\n")
    next
  }

  fit <- readRDS(file_fit)
  metadata <- readRDS(file_metadata)

  sics <- extract_cohort_sics(fit, metadata, sources, type, pt_data,
                               truncated_potentials,
                               band_type = "simultaneous", alpha = 0.1) %>%
    mutate(target = type)

  all_sics_simul[[type]] <- sics
}

sics_data_simul <- bind_rows(all_sics_simul)

# Extract with pointwise bands
all_sics_point <- list()

for (type in targets) {
  cat("Processing:", type, "(pointwise bands)\n")

  file_fit <- paste0(path_compartment, "fit_", make.names(type), "_with_compartment.rds")
  file_metadata <- paste0(path_compartment, "metadata_", make.names(type), "_with_compartment.rds")

  if (!file.exists(file_fit)) next

  fit <- readRDS(file_fit)
  metadata <- readRDS(file_metadata)

  sics <- extract_cohort_sics(fit, metadata, sources, type, pt_data,
                               truncated_potentials,
                               band_type = "pointwise", alpha = 0.1) %>%
    mutate(target = type)

  all_sics_point[[type]] <- sics
}

sics_data_point <- bind_rows(all_sics_point)

# Save extracted data
saveRDS(sics_data_simul, paste0(path_compartment, "sics_compartment_simultaneous.rds"))
saveRDS(sics_data_point, paste0(path_compartment, "sics_compartment_pointwise.rds"))

cat("\nSaved SIC data to:\n")
cat("  -", paste0(path_compartment, "sics_compartment_simultaneous.rds\n"))
cat("  -", paste0(path_compartment, "sics_compartment_pointwise.rds\n"))

# ============================================================================
# 4. CLR vs DII COMPARISON (SIMULTANEOUS BANDS)
# ============================================================================

cat("\n\nCreating CLR vs DII comparison (simultaneous bands)...\n")

plot_data_simul <- sics_data_simul %>%
  filter(distance >= MIN_INTERACTION_RADIUS)

clr_dii_simul_plot <- ggplot(plot_data_simul,
                              aes(x = distance, y = mean, color = group, fill = group)) +
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
    title = "CLR vs DII: Spatial Interaction Curves (Compartment-Adjusted)",
    subtitle = "90% simultaneous credible bands, truncated RBFs (5 basis, σ=12 μm, min=30 μm)",
    x = "Distance (μm)",
    y = "SIC",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9),
    plot.title = element_text(size = 14, face = "bold")
  )

print(clr_dii_simul_plot)

ggsave(
  paste0(path_figures, "clr_dii_simultaneous.pdf"),
  clr_dii_simul_plot,
  width = 14,
  height = 8
)

cat("Saved figure to:", paste0(path_figures, "clr_dii_simultaneous.pdf\n"))

# ============================================================================
# 5. CLR vs DII COMPARISON (POINTWISE BANDS)
# ============================================================================

cat("\nCreating CLR vs DII comparison (pointwise bands)...\n")

plot_data_point <- sics_data_point %>%
  filter(distance >= MIN_INTERACTION_RADIUS)

clr_dii_point_plot <- ggplot(plot_data_point,
                              aes(x = distance, y = mean, color = group, fill = group)) +
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
    title = "CLR vs DII: Spatial Interaction Curves (Compartment-Adjusted)",
    subtitle = "90% pointwise credible bands, truncated RBFs (5 basis, σ=12 μm, min=30 μm)",
    x = "Distance (μm)",
    y = "SIC",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9),
    plot.title = element_text(size = 14, face = "bold")
  )

print(clr_dii_point_plot)

ggsave(
  paste0(path_figures, "clr_dii_pointwise.pdf"),
  clr_dii_point_plot,
  width = 14,
  height = 8
)

cat("Saved figure to:", paste0(path_figures, "clr_dii_pointwise.pdf\n"))

# ============================================================================
# 6. SIMULTANEOUS vs POINTWISE COMPARISON
# ============================================================================

cat("\n\nComparing simultaneous vs pointwise bands...\n")

# Combine both band types
bands_comparison <- bind_rows(
  sics_data_simul %>% mutate(band_label = "Simultaneous (90%)"),
  sics_data_point %>% mutate(band_label = "Pointwise (90%)")
)

# Compute width statistics
band_widths <- bands_comparison %>%
  mutate(width = upper - lower) %>%
  group_by(target, source, group, band_label) %>%
  summarize(
    mean_width = mean(width),
    max_width = max(width),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = band_label,
    values_from = c(mean_width, max_width)
  ) %>%
  mutate(
    mean_width_ratio = `mean_width_Simultaneous (90%)` / `mean_width_Pointwise (90%)`,
    max_width_ratio = `max_width_Simultaneous (90%)` / `max_width_Pointwise (90%)`
  ) %>%
  arrange(desc(mean_width_ratio))

cat("\nBand width comparison (simultaneous / pointwise):\n")
print(band_widths, n = 30)

write_csv(band_widths, paste0(path_compartment, "band_width_comparison.csv"))
cat("\nSaved band width comparison to:", paste0(path_compartment, "band_width_comparison.csv\n"))

# Create comparison plot showing both band types for selected interactions
cat("\nCreating band comparison visualization...\n")

# Select a few interactions to display
example_interactions <- band_widths %>%
  head(6) %>%
  select(target, source, group)

plot_comparison_data <- bands_comparison %>%
  semi_join(example_interactions, by = c("target", "source", "group")) %>%
  filter(distance >= MIN_INTERACTION_RADIUS) %>%
  mutate(facet_label = paste0(source, " → ", target, " (", group, ")"))

band_comparison_plot <- ggplot(plot_comparison_data,
                               aes(x = distance, y = mean, color = band_label, fill = band_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = c("Simultaneous (90%)" = "#2166AC", "Pointwise (90%)" = "#D6604D")
  ) +
  scale_fill_manual(
    values = c("Simultaneous (90%)" = "#2166AC", "Pointwise (90%)" = "#D6604D")
  ) +
  facet_wrap(~ facet_label, ncol = 2, scales = "free_y") +
  labs(
    title = "Simultaneous vs Pointwise Credible Bands",
    subtitle = "Compartment-adjusted models with truncated RBFs",
    x = "Distance (μm)",
    y = "SIC",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9),
    plot.title = element_text(size = 14, face = "bold")
  )

print(band_comparison_plot)

ggsave(
  paste0(path_figures, "band_comparison.pdf"),
  band_comparison_plot,
  width = 10,
  height = 10
)

cat("Saved figure to:", paste0(path_figures, "band_comparison.pdf\n"))

# ============================================================================
# 7. EXTRACT COMPARTMENT EFFECT SIZES
# ============================================================================

cat("\n\nExtracting compartment effect sizes...\n")

compartment_effects <- list()

for (type in targets) {

  file_fit <- paste0(path_compartment, "fit_", make.names(type), "_with_compartment.rds")

  if (!file.exists(file_fit)) next

  fit <- readRDS(file_fit)
  draws <- fit$draws(format = "draws_df")

  # Find compartment coefficient (beta_cov)
  comp_vars <- draws %>%
    select(starts_with("beta_cov")) %>%
    names()

  if (length(comp_vars) > 0) {
    comp_effects <- draws %>%
      select(all_of(comp_vars)) %>%
      summarize(across(everything(), list(
        mean = mean,
        sd = sd,
        q025 = ~ quantile(.x, 0.025),
        q975 = ~ quantile(.x, 0.975)
      ))) %>%
      pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
      separate(variable, into = c("param", "stat"), sep = "_(?=[^_]+$)") %>%
      pivot_wider(names_from = stat, values_from = value) %>%
      mutate(target = type)

    compartment_effects[[type]] <- comp_effects
  }
}

if (length(compartment_effects) > 0) {
  comp_effects_df <- bind_rows(compartment_effects)
  write_csv(comp_effects_df, paste0(path_compartment, "compartment_effects.csv"))

  cat("\n\nCompartment effect sizes:\n")
  cat("=========================\n")
  print(comp_effects_df, n = 20)

  cat("\nSaved to:", paste0(path_compartment, "compartment_effects.csv\n"))
}

# ============================================================================
# 8. SUMMARY
# ============================================================================

cat("\n\n========================================\n")
cat("Compartment model analysis complete!\n")
cat("========================================\n\n")

cat("Output files:\n")
cat("\nData files:\n")
cat("  - SICs (simultaneous):     ", paste0(path_compartment, "sics_compartment_simultaneous.rds\n"))
cat("  - SICs (pointwise):        ", paste0(path_compartment, "sics_compartment_pointwise.rds\n"))
cat("  - Compartment effects:     ", paste0(path_compartment, "compartment_effects.csv\n"))
cat("  - Band width stats:        ", paste0(path_compartment, "band_width_comparison.csv\n"))
cat("\nFigures:\n")
cat("  - CLR vs DII (simul):      ", paste0(path_figures, "clr_dii_simultaneous.pdf\n"))
cat("  - CLR vs DII (pointwise):  ", paste0(path_figures, "clr_dii_pointwise.pdf\n"))
cat("  - Band comparison:         ", paste0(path_figures, "band_comparison.pdf\n"))
cat("\n")
