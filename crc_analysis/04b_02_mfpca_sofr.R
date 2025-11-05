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
# G-CROSS mFPCA/SOFR ANALYSIS
# ============================================================================

cat("\n=== G-cross mFPCA/SOFR Analysis ===\n")

gcross_sofr_results <- expand_grid(
  target = targets,
  source = sources
) %>%
  pmap_dfr(\(target, source) {

    cat("G-cross mFPCA/SOFR:", source, "->", target, "\n")

    tryCatch({
      # Extract G-cross for this pair
      fda <- extract_summary_functions(
        fda_base,
        summary_func = Gcross,
        extract_func = bivariate,
        r_vec = seq(MIN_INTERACTION_RADIUS, 75, by = 1),
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

      # Convert to actual distance scale
      beta_curve <- data.frame(
        r_scaled = plot_data$x,
        beta = plot_data$fit,
        se = plot_data$se
      ) %>%
        mutate(
          r = r_scaled * (75 - MIN_INTERACTION_RADIUS) + MIN_INTERACTION_RADIUS,
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

# ============================================================================
# L-CROSS (K-FUNCTION) mFPCA/SOFR ANALYSIS
# ============================================================================

cat("\n=== L-cross mFPCA/SOFR Analysis ===\n")

lcross_sofr_results <- expand_grid(
  target = targets,
  source = sources
) %>%
  pmap_dfr(\(target, source) {

    cat("L-cross mFPCA/SOFR:", source, "->", target, "\n")

    tryCatch({
      # Extract L-cross for this pair
      fda <- extract_summary_functions(
        fda_base,
        summary_func = Lcross,
        extract_func = bivariate,
        r_vec = seq(MIN_INTERACTION_RADIUS, 75, by = 1),
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

      # Convert to actual distance scale
      beta_curve <- data.frame(
        r_scaled = plot_data$x,
        beta = plot_data$fit,
        se = plot_data$se
      ) %>%
        mutate(
          r = r_scaled * (75 - MIN_INTERACTION_RADIUS) + MIN_INTERACTION_RADIUS,
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

cat("✓ Module 2 complete\n")
