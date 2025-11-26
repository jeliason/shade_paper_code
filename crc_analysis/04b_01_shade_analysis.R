# ============================================================================
# SHADE ANALYSIS: Group-Level Detection and Group Differences
# ============================================================================
# Module 1 of method comparison analysis
# Computes SHADE group-level detections and group differences (CLR vs DII)

cat("\n=== Module 1: SHADE Analysis ===\n")

# ============================================================================
# SHADE DETECTION (SIMULTANEOUS CREDIBLE BANDS)
# ============================================================================

cat("\n=== SHADE Detection using Simultaneous Credible Bands ===\n")

shade_detections <- lapply(targets, \(target) {

  cat("Processing target:", target, "\n")

  file_fit <- paste0(path, "fit_CRC_type_", make.names(target), ".rds")
  file_metadata <- paste0(path, "metadata_CRC_type_", make.names(target), ".rds")

  metadata <- readRDS(file_metadata)
  fit <- readRDS(file_fit)

  spots <- metadata$spots
  types <- metadata$types
  potentials <- metadata$potentials
  coef_names <- metadata$coef_names

  pt_data %>%
    filter(Spot %in% spots) -> pt_df

  draws <- as_draws_rvars(fit$draws())
  beta_global <- draws$beta_global
  rownames(beta_global) <- coef_names

  # Get group names from pt_df (should be 2 groups: CLR and DII)
  group_names <- pt_df %>%
    pull(Group) %>%
    unique() %>%
    sort()

  # For each source, compute simultaneous bands and check for detection
  source_results <- lapply(sources, \(source) {

    if(!(source %in% types)) {
      return(NULL)
    }

    # Extract coefficients for this source-target pair
    ix <- grep(paste0("_", target, "_", source), rownames(beta_global), fixed = TRUE)

    # Create design matrix from potentials
    x_seq <- seq(MIN_INTERACTION_RADIUS, 75, 1)
    x_des <- lapply(potentials, \(pot) pot(x_seq)) %>% do.call(cbind, .)

    # Group-level detections (one for each column of beta_global)
    group_results <- lapply(1:ncol(beta_global), \(col_idx) {

      # Extract group-specific coefficients (column col_idx)
      b_grp <- as.matrix(beta_global)[ix, col_idx]

      # Compute linear predictor
      lp_grp <- x_des %*% b_grp

      bands_grp <- compute_simultaneous_bands(
        lp_data = as.data.frame(lp_grp),
        x_seq = x_seq,
        alpha = 0.05
      )

      grp_detect <- bands_grp %>%
        mutate(excludes_zero = (lower > 0) | (upper < 0)) %>%
        summarise(detected = any(excludes_zero)) %>%
        pull(detected)

      data.frame(
        target = target,
        source = source,
        level = "Group",
        Group = group_names[col_idx],
        detected = grp_detect
      )
    }) %>%
      bind_rows()

    return(group_results)

  }) %>%
    bind_rows()

}) %>%
  bind_rows()

cat("✓ SHADE detections computed\n")
print(shade_detections %>%
        mutate(method = "SHADE") %>%
        select(method, target, source, level, Group, detected))

# ============================================================================
# SHADE GROUP DIFFERENCE ANALYSIS
# ============================================================================

cat("\n=== Computing SHADE Group Differences (CLR vs DII) ===\n")

shade_group_differences <- lapply(targets, \(target) {

  cat("Processing target:", target, "\n")

  file_fit <- paste0(path, "fit_CRC_type_", make.names(target), ".rds")
  file_metadata <- paste0(path, "metadata_CRC_type_", make.names(target), ".rds")

  metadata <- readRDS(file_metadata)
  fit <- readRDS(file_fit)

  spots <- metadata$spots
  types <- metadata$types
  potentials <- metadata$potentials
  coef_names <- metadata$coef_names

  pt_data %>%
    filter(Spot %in% spots) -> pt_df

  draws <- as_draws_rvars(fit$draws())
  beta_global <- draws$beta_global
  rownames(beta_global) <- coef_names

  # Get group names (should be CLR and DII)
  group_names <- pt_df %>%
    pull(Group) %>%
    unique() %>%
    sort()

  # For each source, compute group difference in SICs
  source_results <- lapply(sources, \(source) {

    if(!(source %in% types)) {
      return(NULL)
    }

    # Extract coefficients for this source-target pair
    ix <- grep(paste0("_", target, "_", source), rownames(beta_global), fixed = TRUE)

    # Create design matrix from potentials
    x_seq <- seq(MIN_INTERACTION_RADIUS, 75, 1)
    x_des <- lapply(potentials, \(pot) pot(x_seq)) %>% do.call(cbind, .)

    # Extract group-specific coefficients (columns of beta_global)
    # Assuming column 1 = CLR, column 2 = DII
    b_CLR <- as.matrix(beta_global)[ix, 1]
    b_DII <- as.matrix(beta_global)[ix, 2]

    # Compute group difference: DII - CLR
    b_diff <- b_DII - b_CLR

    # Compute SIC difference
    lp_diff <- x_des %*% b_diff

    # Compute simultaneous 95% credible bands on the difference
    bands_diff <- compute_simultaneous_bands(
      lp_data = as.data.frame(lp_diff),
      x_seq = x_seq,
      alpha = 0.05
    )

    # Detection: does the difference exclude zero anywhere?
    diff_detected <- bands_diff %>%
      mutate(excludes_zero = (lower > 0) | (upper < 0)) %>%
      summarise(detected = any(excludes_zero)) %>%
      pull(detected)

    data.frame(
      target = target,
      source = source,
      detected = diff_detected,
      # Store the bands for later visualization
      bands = I(list(bands_diff))
    )

  }) %>%
    bind_rows()

}) %>%
  bind_rows()

cat("✓ SHADE group differences computed\n")
print(shade_group_differences %>%
        select(target, source, detected))

cat("✓ Module 1 complete\n")
