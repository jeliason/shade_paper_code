# ============================================================================
# VISUALIZATIONS: Create All Figures for Method Comparison
# ============================================================================
# Module 4 of method comparison analysis
# Creates SHADE SIC plots and group difference plots

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

cat("✓ Module 4 complete\n")
