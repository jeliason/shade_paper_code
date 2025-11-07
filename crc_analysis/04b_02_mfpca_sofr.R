# ============================================================================
# MFPCA/SOFR ANALYSIS: G-cross and L-cross Group Comparisons
# ============================================================================
# Module 2 of method comparison analysis
# Runs multilevel functional PCA and scalar-on-function regression for
# G-cross and L-cross to test CLR vs DII group differences

cat("\n=== Module 2: mFPCA/SOFR Analysis ===\n")

# Create mxfda object (same for all analyses)
fda_base <- make_mxfda(
  metadata = pt_data %>% mutate(Group = as.numeric(as.factor(Group)) - 1),
  spatial = bind_rows(dats) %>% rename(x=X, y=Y),
  subject_key = "Patient",
  sample_key = "Spot"
)

# ============================================================================
# G-CROSS mFPCA/SOFR ANALYSIS (COMMENTED OUT - NOT USING)
# ============================================================================

if(FALSE) {
cat("\n=== G-cross mFPCA/SOFR Analysis ===\n")

gcross_sofr_results <- expand_grid(
  target = targets,
  source = sources
) %>%
  pmap_dfr(\(target, source) {

    cat("G-cross mFPCA/SOFR:", source, "->", target, "\n")

    tryCatch({
      # Extract G-cross for this pair
      # Note: mxfda requires r_vec to start at 0
      fda <- extract_summary_functions(
        fda_base,
        summary_func = Gcross,
        extract_func = bivariate,
        r_vec = seq(0, 75, by = 1),
        edge_correction = "km",
        markvar = "type",
        mark1 = target,
        mark2 = source
      )

      # Run multilevel FPCA
      fda <- run_mfpca(
        fda,
        metric = "bi g",
        r = "r",
        value = "fundiff",
        pve = 0.99
      )

      # Run scalar-on-function regression
      fda <- run_sofr(
        fda,
        model_name = "group_test",
        formula = Group ~ s(Patient, bs="re"),
        family = "binomial",
        metric = "bi g",
        r = "r",
        value = "fundiff"
      )

      # Extract model and p-value
      model <- extract_model(fda, 'bi g', type = 'sofr', model_name = 'group_test')
      p_value <- summary(model)$s.table["s(xmat.tmat):L.xmat", "p-value"]
      edf <- summary(model)$s.table["s(xmat.tmat):L.xmat", "edf"]

      # Extract functional coefficient beta(r) for visualization
      plot_obj <- mgcv::plot.gam(model, select = 2, shade = TRUE)
      plot_data <- plot_obj[[2]]

      # Convert to actual distance scale (r_vec goes from 0 to 75)
      beta_curve <- data.frame(
        r_scaled = plot_data$x,
        beta = plot_data$fit,
        se = plot_data$se
      ) %>%
        mutate(
          r = r_scaled * 75,
          lower = beta - 1.96 * se,
          upper = beta + 1.96 * se
        )

      data.frame(
        target = target,
        source = source,
        method = "G-cross (mFPCA)",
        p_value = p_value,
        edf = edf,
        detected = p_value < 0.05,
        beta_curve = I(list(beta_curve)),
        model = I(list(model))
      )

    }, error = \(e) {
      cat("  Error:", conditionMessage(e), "\n")
      data.frame(
        target = target,
        source = source,
        method = "G-cross (mFPCA)",
        p_value = NA,
        edf = NA,
        detected = NA,
        beta_curve = I(list(NULL)),
        model = I(list(NULL))
      )
    })
  })

cat("✓ G-cross mFPCA/SOFR complete\n")
print(gcross_sofr_results %>%
        select(target, source, detected, p_value, edf))
}

# ============================================================================
# G-CROSS mFPCA BY GROUP (FOR VISUAL COMPARISON)
# ============================================================================

cat("\n=== G-cross mFPCA by Group ===\n")

# Create separate mxfda objects for each group
fda_base_CLR <- make_mxfda(
  metadata = pt_data %>% filter(Group == "CLR") %>% mutate(Group = as.numeric(as.factor(Group)) - 1),
  spatial = bind_rows(dats) %>%
    filter(Spot %in% (pt_data %>% filter(Group == "CLR") %>% pull(Spot))) %>%
    rename(x=X, y=Y),
  subject_key = "Patient",
  sample_key = "Spot"
)

fda_base_DII <- make_mxfda(
  metadata = pt_data %>% filter(Group == "DII") %>% mutate(Group = as.numeric(as.factor(Group)) - 1),
  spatial = bind_rows(dats) %>%
    filter(Spot %in% (pt_data %>% filter(Group == "DII") %>% pull(Spot))) %>%
    rename(x=X, y=Y),
  subject_key = "Patient",
  sample_key = "Spot"
)

# Run mFPCA for each target-source pair, extract data and combine
all_mfpca_data <- expand_grid(
  target = targets,
  source = sources
) %>%
  pmap_dfr(function(target, source) {

    cat("G-cross mFPCA:", source, "->", target, "\n")

    # Extract data for both groups
    plot_data_list <- lapply(c("CLR", "DII"), function(grp) {

      fda_base_grp <- if(grp == "CLR") fda_base_CLR else fda_base_DII

      tryCatch({
        fda <- extract_summary_functions(
          fda_base_grp,
          summary_func = Gcross,
          extract_func = bivariate,
          r_vec = seq(0, 75, by = 1),
          edge_correction = "km",
          markvar = "type",
          mark1 = target,
          mark2 = source
        )

        fda <- run_mfpca(fda, metric = "bi g", r = "r", value = "fundiff", pve = 0.99)

        # Extract plot data from the mFPCA plot
        p <- plot(fda, what = "bi g mfpca", level1 = 1, level2 = 1)
        built_plot <- ggplot_build(p[[1]])

        # The mxfda plot has multiple layers - extract data from all layers
        plot_layers <- built_plot$data

        # Combine data from layers
        plot_data <- plot_layers[[1]] %>% select(x, y)
        plot_data$ymin <- plot_layers[[2]]$y
        plot_data$ymax <- plot_layers[[3]]$y

        # Add group, target, source labels
        plot_data$Group <- grp
        plot_data$target <- target
        plot_data$source <- source

        return(plot_data)

      }, error = function(e) {
        cat("  ", grp, "Error:", conditionMessage(e), "\n")
        return(NULL)
      })
    })

    # Combine both groups and return
    bind_rows(plot_data_list)
  })

cat("✓ Extracted mFPCA data for", nrow(all_mfpca_data), "observations\n")

# Create faceted grid plot
cat("\n=== Creating mFPCA faceted grid ===\n")

mfpca_grid <- ggplot(all_mfpca_data, aes(x = x, y = y, color = Group, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(target ~ source) +
  scale_color_manual(values = c("CLR" = "#0072B2", "DII" = "#D55E00")) +
  scale_fill_manual(values = c("CLR" = "#0072B2", "DII" = "#D55E00")) +
  labs(
    x = "Distance (μm)",
    y = "Mean G-cross function",
    # subtitle = "Bands: mean ± 1 SD from first functional PC"
  ) +
  theme(legend.position = "bottom",
legend.title = element_blank())

print(mfpca_grid)
fsave("mfpca_grid_supplement", height = 10, width = 16)

cat("✓ G-cross mFPCA by group complete\n")

# ============================================================================
# L-CROSS mFPCA BY GROUP (FOR VISUAL COMPARISON)
# ============================================================================

cat("\n=== L-cross mFPCA by Group ===\n")

# Run mFPCA for each target-source pair, extract data and combine
all_lcross_mfpca_data <- expand_grid(
  target = targets,
  source = sources
) %>%
  pmap_dfr(function(target, source) {

    cat("L-cross mFPCA:", source, "->", target, "\n")

    # Extract data for both groups
    plot_data_list <- lapply(c("CLR", "DII"), function(grp) {

      fda_base_grp <- if(grp == "CLR") fda_base_CLR else fda_base_DII

      tryCatch({
        fda <- extract_summary_functions(
          fda_base_grp,
          summary_func = Lcross,
          extract_func = bivariate,
          r_vec = seq(0, 75, by = 1),
          edge_correction = "iso",
          markvar = "type",
          mark1 = target,
          mark2 = source
        )

        fda <- run_mfpca(fda, metric = "bi l", r = "r", value = "fundiff", pve = 0.99)

        # Extract plot data from the mFPCA plot
        p <- plot(fda, what = "bi l mfpca", level1 = 1, level2 = 1)
        built_plot <- ggplot_build(p[[1]])

        # The mxfda plot has multiple layers - extract data from all layers
        plot_layers <- built_plot$data

        # Combine data from layers
        plot_data <- plot_layers[[1]] %>% select(x, y)
        plot_data$ymin <- plot_layers[[2]]$y
        plot_data$ymax <- plot_layers[[3]]$y

        # Add group, target, source labels
        plot_data$Group <- grp
        plot_data$target <- target
        plot_data$source <- source

        return(plot_data)

      }, error = function(e) {
        cat("  ", grp, "Error:", conditionMessage(e), "\n")
        return(NULL)
      })
    })

    # Combine both groups and return
    bind_rows(plot_data_list)
  })

cat("✓ Extracted L-cross mFPCA data for", nrow(all_lcross_mfpca_data), "observations\n")

# Create faceted grid plot
cat("\n=== Creating L-cross mFPCA faceted grid ===\n")

lcross_mfpca_grid <- ggplot(all_lcross_mfpca_data, aes(x = x, y = y, color = Group, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(target ~ source) +
  scale_color_manual(values = c("CLR" = "#0072B2", "DII" = "#D55E00")) +
  scale_fill_manual(values = c("CLR" = "#0072B2", "DII" = "#D55E00")) +
  labs(
    x = "Distance (μm)",
    y = "Mean L-cross function"
  ) +
  theme(legend.position = "bottom",
        legend.title = element_blank())

print(lcross_mfpca_grid)
fsave("lcross_mfpca_grid_supplement", height = 10, width = 16)

cat("✓ L-cross mFPCA by group complete\n")

# ============================================================================
# SUMMARY TABLES AT KEY DISTANCE INTERVALS
# ============================================================================

cat("\n=== Creating Summary Tables ===\n")

# Define key distances for summary
key_distances <- c(10, 25, 50, 75)

# G-cross summary
gcross_summary <- all_mfpca_data %>%
  filter(x %in% key_distances) %>%
  mutate(CI_excludes_zero = (ymin > 0) | (ymax < 0)) %>%
  select(target, source, Group, distance = x, mean = y, lower = ymin, upper = ymax, CI_excludes_zero) %>%
  arrange(target, source, Group, distance)

cat("\n### G-cross mFPCA Summary at Key Distances ###\n")
print(gcross_summary, n = 100)

# L-cross summary
lcross_summary <- all_lcross_mfpca_data %>%
  filter(x %in% key_distances) %>%
  mutate(CI_excludes_zero = (ymin > 0) | (ymax < 0)) %>%
  select(target, source, Group, distance = x, mean = y, lower = ymin, upper = ymax, CI_excludes_zero) %>%
  arrange(target, source, Group, distance)

cat("\n### L-cross mFPCA Summary at Key Distances ###\n")
print(lcross_summary, n = 100)

# Find ranges where CI excludes zero
cat("\n### G-cross: Distance Ranges Where CI Excludes Zero ###\n")
gcross_exclusion_ranges <- all_mfpca_data %>%
  mutate(CI_excludes_zero = (ymin > 0) | (ymax < 0)) %>%
  filter(CI_excludes_zero) %>%
  group_by(target, source, Group) %>%
  summarise(
    min_dist = min(x),
    max_dist = max(x),
    n_points = n(),
    .groups = "drop"
  )
print(gcross_exclusion_ranges)

cat("\n### L-cross: Distance Ranges Where CI Excludes Zero ###\n")
lcross_exclusion_ranges <- all_lcross_mfpca_data %>%
  mutate(CI_excludes_zero = (ymin > 0) | (ymax < 0)) %>%
  filter(CI_excludes_zero) %>%
  group_by(target, source, Group) %>%
  summarise(
    min_dist = min(x),
    max_dist = max(x),
    n_points = n(),
    .groups = "drop"
  )
print(lcross_exclusion_ranges)

# Save summaries
saveRDS(list(
  gcross_summary = gcross_summary,
  lcross_summary = lcross_summary,
  gcross_exclusion_ranges = gcross_exclusion_ranges,
  lcross_exclusion_ranges = lcross_exclusion_ranges
), paste0(path, "mfpca_summaries.rds"))

cat("\n✓ Summary tables saved\n")

# ============================================================================
# L-CROSS (K-FUNCTION) mFPCA/SOFR ANALYSIS (COMMENTED OUT - NOT USING)
# ============================================================================

if(FALSE) {
cat("\n=== L-cross mFPCA/SOFR Analysis ===\n")

lcross_sofr_results <- expand_grid(
  target = targets,
  source = sources
) %>%
  pmap_dfr(\(target, source) {

    cat("L-cross mFPCA/SOFR:", source, "->", target, "\n")

    tryCatch({
      # Extract L-cross for this pair
      # Note: mxfda requires r_vec to start at 0
      fda <- extract_summary_functions(
        fda_base,
        summary_func = Lcross,
        extract_func = bivariate,
        r_vec = seq(0, 75, by = 1),
        edge_correction = "iso",
        markvar = "type",
        mark1 = target,
        mark2 = source
      )

      # Run multilevel FPCA
      fda <- run_mfpca(
        fda,
        metric = "bi l",
        r = "r",
        value = "fundiff",
        pve = 0.99
      )

      # Run scalar-on-function regression
      fda <- run_sofr(
        fda,
        model_name = "group_test",
        formula = Group ~ s(Patient, bs="re"),
        family = "binomial",
        metric = "bi l",
        r = "r",
        value = "fundiff"
      )

      # Extract model and p-value
      model <- extract_model(fda, 'bi l', type = 'sofr', model_name = 'group_test')
      p_value <- summary(model)$s.table["s(xmat.tmat):L.xmat", "p-value"]
      edf <- summary(model)$s.table["s(xmat.tmat):L.xmat", "edf"]

      # Extract functional coefficient beta(r) for visualization
      plot_obj <- mgcv::plot.gam(model, select = 2, shade = TRUE)
      plot_data <- plot_obj[[2]]

      # Convert to actual distance scale (r_vec goes from 0 to 75)
      beta_curve <- data.frame(
        r_scaled = plot_data$x,
        beta = plot_data$fit,
        se = plot_data$se
      ) %>%
        mutate(
          r = r_scaled * 75,
          lower = beta - 1.96 * se,
          upper = beta + 1.96 * se
        )

      data.frame(
        target = target,
        source = source,
        method = "L-cross (mFPCA)",
        p_value = p_value,
        edf = edf,
        detected = p_value < 0.05,
        beta_curve = I(list(beta_curve)),
        model = I(list(NULL))
      )

    }, error = \(e) {
      cat("  Error:", conditionMessage(e), "\n")
      data.frame(
        target = target,
        source = source,
        method = "L-cross (mFPCA)",
        p_value = NA,
        edf = NA,
        detected = NA,
        beta_curve = I(list(NULL)),
        model = I(list(NULL))
      )
    })
  })

cat("✓ L-cross mFPCA/SOFR complete\n")
print(lcross_sofr_results %>%
        select(target, source, detected, p_value, edf))
}  # End if(FALSE) for L-cross SOFR

cat("✓ Module 2 complete\n")
