# ============================================================================
# VISUALIZATIONS: Create All Figures for Method Comparison
# ============================================================================
# Module 4 of method comparison analysis
# Creates SHADE SIC plots, group difference plots, mFPCA examples, and comparison heatmaps

cat("\n=== Module 4: Visualizations ===\n")

# ============================================================================
# SHADE SIC VISUALIZATIONS (BY GROUP)
# ============================================================================

cat("\n=== Creating SHADE SIC visualizations ===\n")

# For each target, create a plot showing all source SICs
for (targ in targets) {

  cat("Plotting SHADE SICs for target:", targ, "\n")

  file_fit <- paste0(path, "fit_CRC_type_", make.names(targ), ".rds")
  file_metadata <- paste0(path, "metadata_CRC_type_", make.names(targ), ".rds")

  metadata <- readRDS(file_metadata)
  fit <- readRDS(file_fit)

  types <- metadata$types
  potentials <- metadata$potentials
  coef_names <- metadata$coef_names

  draws <- as_draws_rvars(fit$draws())
  beta_global <- draws$beta_global
  rownames(beta_global) <- coef_names

  # Get group names
  pt_data %>%
    filter(Spot %in% metadata$spots) %>%
    pull(Group) %>%
    unique() %>%
    sort() -> group_names

  # Create design matrix
  x_seq <- seq(MIN_INTERACTION_RADIUS, 75, 1)
  x_des <- lapply(potentials, \(pot) pot(x_seq)) %>% do.call(cbind, .)

  # Compute SICs for all sources and both groups
  sic_data <- lapply(sources, \(src) {

    if(!(src %in% types)) {
      return(NULL)
    }

    # Extract coefficients for this source-target pair
    ix <- grep(paste0("_", targ, "_", src), rownames(beta_global), fixed = TRUE)

    # Compute for each group
    lapply(1:ncol(beta_global), \(col_idx) {

      b_grp <- as.matrix(beta_global)[ix, col_idx]
      lp_grp <- x_des %*% b_grp

      bands_grp <- compute_simultaneous_bands(
        lp_data = as.data.frame(lp_grp),
        x_seq = x_seq,
        alpha = 0.05
      )

      # Add detection info
      detected <- any(bands_grp$lower > 0 | bands_grp$upper < 0)

      bands_grp %>%
        mutate(
          source = src,
          Group = group_names[col_idx],
          detected = detected
        )

    }) %>% bind_rows()

  }) %>%
    bind_rows()

  # Create plot
  p <- sic_data %>%
    ggplot(aes(x = x, color = Group, fill = Group)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.5) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    geom_line(aes(y = mean), linewidth = 1) +
    facet_wrap(~source, ncol = 2, scales = "free_y") +
    scale_color_manual(values = c("CLR" = "#0072B2", "DII" = "#D55E00")) +
    scale_fill_manual(values = c("CLR" = "#0072B2", "DII" = "#D55E00")) +
    labs(
      title = paste("SHADE Group-Level SICs:", targ, "← Sources"),
      subtitle = "Ribbons show 95% simultaneous credible bands",
      x = "Distance (μm)",
      y = "Spatial Interaction Curve"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA),
      strip.background = element_rect(fill = "gray90", color = "black")
    )

  print(p)
  fsave(paste0("shade_sics_", make.names(targ)), height = 8, width = 10)

}

# ============================================================================
# METHOD COMPARISON HEATMAP
# ============================================================================

cat("\n=== Creating method comparison heatmap ===\n")

comparison_plot <- all_detections %>%
  mutate(detected = ifelse(detected, "Yes", "No")) %>%
  ggplot(aes(x = source, y = target, fill = detected)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_wrap(~method, ncol = 3) +
  scale_fill_manual(values = c("No" = "gray90", "Yes" = "#D55E00"), name = "Detected") +
  labs(
    x = "Source cell type",
    y = "Target cell type",
    title = "Method Comparison: Group-Level Detection"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_rect(fill = "gray90", color = "black")
  )

print(comparison_plot)
fsave("method_comparison_sofr", height = 6, width = 12)

# ============================================================================
# EXAMPLE mFPCA/SOFR PLOT (CHOOSE ONE SIGNIFICANT PAIR)
# ============================================================================

cat("\n=== Creating example mFPCA/SOFR visualization ===\n")

# Find a significant pair from G-cross results
sig_pair <- gcross_sofr_results %>%
  filter(detected == TRUE, !is.na(p_value)) %>%
  arrange(p_value) %>%
  slice(1)

if (nrow(sig_pair) > 0) {
  cat("Plotting mFPCA example for:", sig_pair$source, "->", sig_pair$target, "\n")

  beta_curve_data <- sig_pair$beta_curve[[1]]

  p_mfpca <- ggplot(beta_curve_data, aes(x = r, y = beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "steelblue") +
    geom_line(color = "steelblue", linewidth = 1) +
    labs(
      x = "Distance (μm)",
      y = "Functional coefficient β(r)",
      title = paste("G-cross mFPCA/SOFR:", sig_pair$source, "→", sig_pair$target),
      subtitle = paste0("p = ", round(sig_pair$p_value, 4), ", edf = ", round(sig_pair$edf, 2))
    ) +
    theme_bw()

  print(p_mfpca)
  fsave("mfpca_example_gcross", height = 5, width = 7)
}

cat("✓ Module 4 complete\n")
