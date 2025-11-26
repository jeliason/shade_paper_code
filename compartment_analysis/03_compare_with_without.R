# Compare Original vs Compartment-Adjusted Models
#
# Purpose: Compare SHADE estimates with and without compartment covariate
# Shows how SIC estimates change when accounting for spatial compartments

library(tidyverse)
library(posterior)
library(patchwork)
source("utils.R")  # For compute_simultaneous_bands()

# ============================================================================
# 1. CONFIGURATION
# ============================================================================

# Load patient metadata
pt_data <- read_csv("crc_analysis/data/CRC_pt_metadata.csv")

targets <- c("CTLs", "memory CD4+ T", "granulocytes")
sources <- c("TAMs", "CAFs", "vasculature", "hybrid E/M", "tumor cells")

path_original <- "./crc_analysis/data/"
path_compartment <- "./compartment_analysis/data/"

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

  # Extract draws using rvar format (matches crc_analysis approach)
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
    # Find coefficients for this target-source pair (matches crc_analysis grep pattern)
    ix <- grep(paste0("_", target_type, "_", source), rownames(beta_global), fixed = TRUE)

    if (length(ix) == 0) {
      warning(paste0("No coefficients found for ", target_type, " -> ", source))
      next
    }

    # Extract coefficients
    b_g <- as.matrix(beta_global)[ix, ]

    # Compute linear predictor via matrix multiplication
    lp_g <- x_des %*% b_g

    # The global level estimates are in the last columns (one per group)
    # Extract only the cohort-level columns
    num_groups <- length(groups)
    cohort_cols <- (ncol(lp_g) - num_groups + 1):ncol(lp_g)
    lp_cohort <- lp_g[, cohort_cols, drop = FALSE]

    # Set column names to group names
    colnames(lp_cohort) <- groups

    # Compute credible bands (simultaneous or pointwise)
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
# 3. LOAD AND COMPARE ALL MODELS
# ============================================================================

cat("Comparing original vs compartment-adjusted models...\n")
cat("====================================================\n\n")

all_comparisons <- list()

for (type in targets) {

  cat("Processing:", type, "\n")

  # Load original fit (from crc_analysis)
  file_fit_orig <- paste0(path_original, "fit_CRC_type_", make.names(type), ".rds")
  file_metadata_orig <- paste0(path_original, "metadata_CRC_type_", make.names(type), ".rds")

  # Load compartment fit (from compartment_analysis)
  file_fit_comp <- paste0(path_compartment, "fit_", make.names(type), "_with_compartment.rds")
  file_metadata_comp <- paste0(path_compartment, "metadata_", make.names(type), "_with_compartment.rds")

  if (!file.exists(file_fit_orig)) {
    cat("  Warning: Original fit not found:", file_fit_orig, "\n")
    next
  }

  if (!file.exists(file_fit_comp)) {
    cat("  Warning: Compartment fit not found:", file_fit_comp, "\n")
    next
  }

  # Load fits
  fit_orig <- readRDS(file_fit_orig)
  metadata_orig <- readRDS(file_metadata_orig)

  fit_comp <- readRDS(file_fit_comp)
  metadata_comp <- readRDS(file_metadata_comp)

  # Extract SICs
  cat("  Extracting SICs from original model...\n")
  sics_orig <- extract_cohort_sics(fit_orig, metadata_orig, sources, type, pt_data,
                                   truncated_potentials) %>%
    mutate(model = "original")

  cat("  Extracting SICs from compartment-adjusted model...\n")
  sics_comp <- extract_cohort_sics(fit_comp, metadata_comp, sources, type, pt_data,
                                   truncated_potentials) %>%
    mutate(model = "compartment_adjusted")

  # Combine
  comparison <- bind_rows(sics_orig, sics_comp) %>%
    mutate(target = type)

  all_comparisons[[type]] <- comparison
}

comparison_data <- bind_rows(all_comparisons)

# Save for later use
saveRDS(comparison_data, paste0(path_compartment, "comparison_sics.rds"))

cat("\nSaved SIC comparisons to:", paste0(path_compartment, "comparison_sics.rds\n"))

# ============================================================================
# 4. SUMMARY STATISTICS
# ============================================================================

cat("\n\nComputing summary statistics...\n")

summary_stats <- comparison_data %>%
  # Use 'mean' column from simultaneous bands (already computed)
  select(target, source, group, distance, model, mean) %>%
  pivot_wider(names_from = model, values_from = mean) %>%
  mutate(
    change = compartment_adjusted - original,
    abs_change = abs(change),
    pct_change = 100 * change / (abs(original) + 1e-6)
  ) %>%
  group_by(target, source, group) %>%
  summarize(
    max_original = max(abs(original)),
    max_adjusted = max(abs(compartment_adjusted)),
    mean_abs_change = mean(abs_change),
    max_abs_change = max(abs_change),
    mean_pct_change = mean(abs(pct_change)),
    .groups = "drop"
  ) %>%
  arrange(desc(max_abs_change))

cat("\n\nSummary: Changes in SIC estimates\n")
cat("==================================\n\n")
print(summary_stats, n = 30)

# Save table
write_csv(summary_stats, paste0(path_compartment, "comparison_summary.csv"))

cat("\nSaved summary table to:", paste0(path_compartment, "comparison_summary.csv\n"))

# ============================================================================
# 5. VISUALIZATIONS
# ============================================================================

cat("\n\nCreating visualizations...\n")

# Filter to top changing pairs for visualization
top_pairs <- summary_stats %>%
  head(9) %>%
  select(target, source, group)

plot_data <- comparison_data %>%
  semi_join(top_pairs, by = c("target", "source", "group")) %>%
  filter(distance >= MIN_INTERACTION_RADIUS) %>%
  mutate(facet_label = paste0(source, " → ", target, "\n(", group, ")"))

# Create single plot with facet_wrap
combined <- ggplot(plot_data, aes(x = distance, y = mean, color = model, fill = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = c(original = "black", compartment_adjusted = "#2166AC"),
    labels = c(original = "Original", compartment_adjusted = "Compartment-adjusted")
  ) +
  scale_fill_manual(
    values = c(original = "black", compartment_adjusted = "#2166AC"),
    labels = c(original = "Original", compartment_adjusted = "Compartment-adjusted")
  ) +
  facet_wrap(~ facet_label, ncol = 3, scales = "free_y") +
  labs(
    title = "Effect of compartment adjustment on SIC estimates",
    subtitle = "Binary compartments (tumor/stromal) defined by tumor cell density threshold",
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
combined

# ============================================================================
# CREATE CLR vs DII COMPARISON PLOT (COMPARTMENT-ADJUSTED MODELS ONLY)
# ============================================================================

cat("\nCreating CLR vs DII comparison for compartment-adjusted models...\n")

# Filter to compartment-adjusted models only
clr_dii_plot_data <- comparison_data %>%
  filter(model == "compartment_adjusted") %>%
  filter(distance >= MIN_INTERACTION_RADIUS)

# Create facet_grid plot: target x source
clr_dii_plot <- ggplot(clr_dii_plot_data, aes(x = distance, y = mean, color = group, fill = group)) +
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
    title = "CLR vs DII group comparison (compartment-adjusted models)",
    subtitle = "Spatial interaction curves with binary compartment covariates",
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

clr_dii_plot

ggsave(
  "./compartment_analysis/compartment_analysis_figures/clr_dii_comparison.pdf",
  clr_dii_plot,
  width = 14,
  height = 8
)

cat("Saved figure to: ./compartment_analysis/compartment_analysis_figures/clr_dii_comparison.pdf\n")

# ============================================================================
# CREATE CLR vs DII COMPARISON PLOT (ORIGINAL MODELS)
# ============================================================================

cat("\nCreating CLR vs DII comparison for original models...\n")

# Filter to original models only
clr_dii_orig_plot_data <- comparison_data %>%
  filter(model == "original") %>%
  filter(distance >= MIN_INTERACTION_RADIUS)

# Create facet_grid plot: target x source
clr_dii_orig_plot <- ggplot(clr_dii_orig_plot_data, aes(x = distance, y = mean, color = group, fill = group)) +
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
    title = "CLR vs DII group comparison (original models)",
    subtitle = "Spatial interaction curves without compartment covariates",
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

clr_dii_orig_plot

ggsave(
  "./compartment_analysis/compartment_analysis_figures/clr_dii_original.pdf",
  clr_dii_orig_plot,
  width = 14,
  height = 8
)

cat("Saved figure to: ./compartment_analysis/compartment_analysis_figures/clr_dii_original.pdf\n")

# ============================================================================
# 6. EXTRACT COMPARTMENT EFFECT SIZES
# ============================================================================

cat("\n\nExtracting compartment effect sizes...\n")

compartment_effects <- list()

for (type in targets) {

  file_fit_comp <- paste0(path_compartment, "fit_", make.names(type), "_with_compartment.rds")

  if (!file.exists(file_fit_comp)) next

  fit_comp <- readRDS(file_fit_comp)
  draws <- fit_comp$draws(format = "draws_df")

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
# 7. COMPARE SIMULTANEOUS VS POINTWISE BANDS
# ============================================================================

cat("\n\n========================================\n")
cat("Comparing simultaneous vs pointwise bands\n")
cat("========================================\n\n")

# Extract SICs with pointwise bands for one target type as example
cat("Extracting pointwise bands for CTLs...\n")

fit_comp_ctls <- readRDS(paste0(path_compartment, "fit_CTLs_with_compartment.rds"))
metadata_comp_ctls <- readRDS(paste0(path_compartment, "metadata_CTLs_with_compartment.rds"))

# Extract with both band types
sics_simul <- extract_cohort_sics(fit_comp_ctls, metadata_comp_ctls, sources, "CTLs",
                                  pt_data, truncated_potentials,
                                  band_type = "simultaneous", alpha = 0.1) %>%
  mutate(target = "CTLs")

sics_point <- extract_cohort_sics(fit_comp_ctls, metadata_comp_ctls, sources, "CTLs",
                                  pt_data, truncated_potentials,
                                  band_type = "pointwise", alpha = 0.1) %>%
  mutate(target = "CTLs")

# Combine
bands_comparison <- bind_rows(sics_simul, sics_point)

# Save comparison data
saveRDS(bands_comparison, paste0(path_compartment, "bands_comparison_ctls.rds"))

# Compute width statistics
band_widths <- bands_comparison %>%
  mutate(width = upper - lower) %>%
  group_by(target, source, group, band_type) %>%
  summarize(
    mean_width = mean(width),
    max_width = max(width),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = band_type,
    values_from = c(mean_width, max_width)
  ) %>%
  mutate(
    mean_width_ratio = mean_width_simultaneous / mean_width_pointwise,
    max_width_ratio = max_width_simultaneous / max_width_pointwise
  ) %>%
  arrange(desc(mean_width_ratio))

cat("\nBand width comparison (simultaneous / pointwise):\n")
print(band_widths, n = 30)

write_csv(band_widths, paste0(path_compartment, "band_width_comparison.csv"))
cat("\nSaved band width comparison to:", paste0(path_compartment, "band_width_comparison.csv\n"))

# Create comparison plot showing both band types
cat("\nCreating band comparison visualization...\n")

# Select a few interactions to display
example_interactions <- band_widths %>%
  head(6) %>%
  select(source, group)

plot_comparison_data <- bands_comparison %>%
  semi_join(example_interactions, by = c("source", "group")) %>%
  filter(distance >= MIN_INTERACTION_RADIUS) %>%
  mutate(
    facet_label = paste0(source, " → CTLs (", group, ")"),
    band_label = ifelse(band_type == "simultaneous", "Simultaneous (90%)", "Pointwise (90%)")
  )

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
    subtitle = "Simultaneous bands are wider to ensure entire curve coverage",
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

band_comparison_plot

ggsave(
  "./compartment_analysis/compartment_analysis_figures/band_comparison.pdf",
  band_comparison_plot,
  width = 10,
  height = 10
)

cat("Saved figure to: ./compartment_analysis/compartment_analysis_figures/band_comparison.pdf\n")

# ============================================================================
# 8. CLR vs DII COMPARISON WITH POINTWISE BANDS
# ============================================================================

cat("\n\n========================================\n")
cat("Creating CLR vs DII comparison with pointwise bands\n")
cat("========================================\n\n")

# Extract all targets with pointwise bands (compartment models)
all_comparisons_pointwise <- list()

for (type in targets) {
  cat("Processing:", type, "\n")

  file_fit_comp <- paste0(path_compartment, "fit_", make.names(type), "_with_compartment.rds")
  file_metadata_comp <- paste0(path_compartment, "metadata_", make.names(type), "_with_compartment.rds")

  if (!file.exists(file_fit_comp)) {
    cat("  Warning: Compartment fit not found:", file_fit_comp, "\n")
    next
  }

  fit_comp <- readRDS(file_fit_comp)
  metadata_comp <- readRDS(file_metadata_comp)

  # Extract SICs with pointwise bands
  cat("  Extracting SICs with pointwise bands...\n")
  sics_point <- extract_cohort_sics(fit_comp, metadata_comp, sources, type, pt_data,
                                    truncated_potentials,
                                    band_type = "pointwise", alpha = 0.1) %>%
    mutate(target = type)

  all_comparisons_pointwise[[type]] <- sics_point
}

comparison_data_pointwise <- bind_rows(all_comparisons_pointwise)

# Save pointwise comparison data
saveRDS(comparison_data_pointwise, paste0(path_compartment, "comparison_sics_pointwise.rds"))
cat("\nSaved pointwise SIC data to:", paste0(path_compartment, "comparison_sics_pointwise.rds\n"))

# Create CLR vs DII facet_grid plot with pointwise bands
cat("\nCreating CLR vs DII facet_grid with pointwise bands...\n")

clr_dii_pointwise_data <- comparison_data_pointwise %>%
  filter(distance >= MIN_INTERACTION_RADIUS)

clr_dii_pointwise_plot <- ggplot(clr_dii_pointwise_data,
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
    title = "CLR vs DII group comparison (compartment-adjusted, pointwise bands)",
    subtitle = "Spatial interaction curves with binary compartment covariates (90% pointwise CI)",
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

clr_dii_pointwise_plot

ggsave(
  "./compartment_analysis/compartment_analysis_figures/clr_dii_pointwise.pdf",
  clr_dii_pointwise_plot,
  width = 14,
  height = 8
)

cat("Saved figure to: ./compartment_analysis/compartment_analysis_figures/clr_dii_pointwise.pdf\n")

# ============================================================================
# 9. SUMMARY
# ============================================================================

cat("\n\n========================================\n")
cat("Compartment comparison complete!\n")
cat("========================================\n\n")

cat("Output files:\n")
cat("\nData files:\n")
cat("  - SIC comparisons (simultaneous): ./compartment_analysis/data/comparison_sics.rds\n")
cat("  - SIC comparisons (pointwise):    ./compartment_analysis/data/comparison_sics_pointwise.rds\n")
cat("  - Summary table:                  ./compartment_analysis/data/comparison_summary.csv\n")
cat("  - Compartment effects:            ./compartment_analysis/data/compartment_effects.csv\n")
cat("  - Band comparison:                ./compartment_analysis/data/bands_comparison_ctls.rds\n")
cat("  - Band width stats:               ./compartment_analysis/data/band_width_comparison.csv\n")
cat("\nFigures:\n")
cat("  - Original vs compartment:            ./compartment_analysis/compartment_analysis_figures/compartment_comparison.pdf\n")
cat("  - CLR vs DII (compartment, simul):    ./compartment_analysis/compartment_analysis_figures/clr_dii_comparison.pdf\n")
cat("  - CLR vs DII (compartment, pointwise): ./compartment_analysis/compartment_analysis_figures/clr_dii_pointwise.pdf\n")
cat("  - CLR vs DII (original, simul):       ./compartment_analysis/compartment_analysis_figures/clr_dii_original.pdf\n")
cat("  - Simultaneous vs pointwise:          ./compartment_analysis/compartment_analysis_figures/band_comparison.pdf\n")
cat("\n")
