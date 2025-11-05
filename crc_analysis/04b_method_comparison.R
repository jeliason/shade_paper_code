library(tidyverse)
library(spatstat)
library(Matrix)
library(cmdstanr)
library(SHADE)
library(posterior)
library(mxfda)

source("crc_analysis/utils.R")
source("utils.R")

theme_set(theme_bw(base_size=14, base_family='Helvetica')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

fsave <- \(fname,height=5,width=8) {
  ggsave(paste0(figures_folder,fname,".pdf"),device=cairo_pdf, height=height, width=width, units="in")
}
figures_folder <- "./manuscript/images/CRC_analysis_paper/"

# ============================================================================
# LOAD DATA
# ============================================================================

df_raw <- read_csv("crc_analysis/data/CRC_cleaned.csv")
pt_data <- read_csv("crc_analysis/data/CRC_pt_metadata.csv")

df_raw <- df_raw %>%
  dplyr::mutate(type = as.factor(type)) %>%
  mutate(type = fct_recode(type,"CAFs"="smooth muscle","hybrid E/M"="stroma","TAMs"="CD163+ macros","CTLs"="CD8+ T cells"))

targets <- c(
  "CTLs",
  "memory CD4+ T",
  "granulocytes"
)

sources <- c(
  "TAMs",
  "CAFs",
  "vasculature",
  "hybrid E/M",
  "tumor cells"
)

path <- "./crc_analysis/data/"
types_keep <- c(targets, sources)
num_types <- length(types_keep)

# ============================================================================
# PREPARE SPATIAL POINT PATTERNS
# ============================================================================

spots <- pt_data %>% pull(Spot)

dats <- df_raw %>%
  dplyr::filter(Spot %in% spots) %>%
  group_by(Spot) %>%
  group_map(~{
    .x <- .x %>%
      filter(type %in% types_keep) %>%
      droplevels() %>%
      mutate(Spot = .y$Spot)
  })

pats <- lapply(1:length(dats), \(i) {
  df <- dats[[i]]
  tryCatch({
    pat <- make_pat(df$X, df$Y, factor(df$type, levels=types_keep))
    sq_W <- owin(xrange = c(min(df$X), max(df$X)), yrange=c(min(df$Y), max(df$Y)))
    Window(pat) <- sq_W
    pat
  }, error=\(e) {
    print(paste("Error in pattern", i))
    stop(e)
  })
})

# Add metadata to patterns
for(i in seq_along(pats)) {
  spot <- unique(dats[[i]]$Spot)
  pt_info <- pt_data %>% filter(Spot == spot)
  attr(pats[[i]], "Spot") <- spot
  attr(pats[[i]], "Patient") <- pt_info$Patient
  attr(pats[[i]], "Group") <- pt_info$Group
}

cat("Loaded", length(pats), "spatial patterns\n")

# ============================================================================
# PART 1: SHADE DETECTION (SIMULTANEOUS CREDIBLE BANDS)
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

# ============================================================================
# SHADE VISUALIZATION: SIC CURVES WITH DETECTIONS
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
# PART 2: G-CROSS GLOBAL ENVELOPE TEST
# ============================================================================

fda <- make_mxfda(metadata = pt_data %>% mutate(Group = as.numeric(as.factor(Group)) - 1),
                         spatial = bind_rows(dats) %>% rename(x=X,y=Y),
                         subject_key = "Patient",
                         sample_key = "Spot")

fda <- extract_summary_functions(fda,
                summary_func = Gcross,
                extract_func = bivariate,
                r_vec = seq(0, 100, by = 1),
                edge_correction = "km",
                markvar = "type",
                mark1 = "CTLs",
                mark2 = "TAMs")

plot(fda, y = "fundiff", what = "bi g", sampleID = "Spot") +
  geom_hline(yintercept = 0, color = "red", linetype = 2)

fda <- run_mfpca(fda, 
                         metric = "bi g", 
                         r = "r", 
                         value = "fundiff",
                         pve = .99)

fda <- run_sofr(
    fda,
    model_name = "group_test",
    formula = Group ~ s(Patient, bs="re"),  # Random effect for Patient
    family = "binomial",
    metric = "bi g",
    r = "r",
    value = "fundiff"
  )

model = extract_model(fda, 'bi g', type = 'sofr', model_name = 'group_test')
# plot(model)

  p <- mgcv::plot.gam(model, select  =2, shade = TRUE,
       xlab = "Distance (μm)",
       ylab = "Functional coefficient β(r)",
       main = "Effect of G-cross on log-odds of Group")
summary(model)

  # Simpler approach: extract from the plot object
  plot_data <- p[[2]]  # This contains the plotting data

  # Convert to dataframe for ggplot
  df_plot <- data.frame(
    r_scaled = plot_data$x,
    beta = plot_data$fit,
    se = plot_data$se
  ) %>%
    mutate(
      # Map scaled [0,1] back to actual distances
      r = r_scaled * (75 - MIN_INTERACTION_RADIUS) +
  MIN_INTERACTION_RADIUS,
      lower = beta - 1.96 * se,
      upper = beta + 1.96 * se
    )

  # Create ggplot
  ggplot(df_plot, aes(x = r, y = beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill =
  "steelblue") +
    geom_line(color = "steelblue", linewidth = 1) +
    labs(
      x = "Distance (μm)",
      y = "Functional coefficient β(r)",
      title = "Effect of G-cross on log-odds of Group",
      subtitle = paste("p =", round(summary(model)$s.table[2, "p-value"],
   3))
    ) +
    theme_bw()
cat("\n=== G-cross Global Envelope Tests ===\n")

gcross_detections <- expand_grid(
  target = targets,
  source = sources
) %>%
  pmap_dfr(\(target, source) {

    cat("G-cross:", source, "->", target, "\n")

    # Test at image level
    image_results <- lapply(seq_along(pats), \(i) {

      pat <- pats[[i]]
      spot <- attr(pat, "Spot")
      patient <- attr(pat, "Patient")
      group <- attr(pat, "Group")

      tryCatch({
        # Compute G-cross with global envelope
        # Match parameters from sim_shade_comparison
        env <- envelope(
          pat,
          fun = Gcross,
          nsim = 99,
          nrank = 10,
          i = source,
          j = target,
          correction = "km",
          r = seq(0, 75, by = 1),
          global = TRUE,
          savefuns = FALSE,
          verbose = FALSE,
          silent=TRUE
        )

        # Detection: does observed function exit envelope anywhere?
        # Only check in same distance range as SHADE (MIN_INTERACTION_RADIUS to 75)
        r_check <- env$r >= MIN_INTERACTION_RADIUS
        detected <- any((env$obs > env$hi | env$obs < env$lo)[r_check], na.rm = TRUE)

        data.frame(
          Spot = spot,
          Patient = patient,
          Group = group,
          detected = detected
        )

      }, error = \(e) {
        # Return NA if error, just note the pattern index
        cat("  Error in pattern", i, "\n")
        data.frame(
          Spot = spot,
          Patient = patient,
          Group = group,
          detected = NA
        )
      })

    }) %>%
      bind_rows()

    # Summarize detections
    # Image-level: proportion of images with detection
    image_summary <- image_results %>%
      summarise(
        n_images = n(),
        n_detected = sum(detected, na.rm = TRUE),
        prop_detected = n_detected / n_images,
        detected = prop_detected > 0.5  # Majority rule
      ) %>%
      mutate(
        target = target,
        source = source,
        level = "Image",
        Group = "All"
      )

    # Group-level: proportion within each group
    group_summary <- image_results %>%
      group_by(Group) %>%
      summarise(
        n_images = n(),
        n_detected = sum(detected, na.rm = TRUE),
        prop_detected = n_detected / n_images,
        detected = prop_detected > 0.5
      ) %>%
      mutate(
        target = target,
        source = source,
        level = "Group"
      )

    bind_rows(image_summary, group_summary)

  })

cat("✓ G-cross detections computed\n")
print(gcross_detections %>%
        mutate(method = "G-cross") %>%
        select(method, target, source, level, Group, detected, prop_detected) %>%
        filter(level == "Group"))

# ============================================================================
# PART 3: K-CROSS GLOBAL ENVELOPE TEST
# ============================================================================

cat("\n=== K-cross Global Envelope Tests ===\n")

kcross_detections <- expand_grid(
  target = targets,
  source = sources
) %>%
  pmap_dfr(\(target, source) {

    cat("K-cross:", source, "->", target, "\n")

    # Test at image level
    image_results <- lapply(seq_along(pats), \(i) {

      pat <- pats[[i]]
      spot <- attr(pat, "Spot")
      patient <- attr(pat, "Patient")
      group <- attr(pat, "Group")

      tryCatch({
        # Compute L-cross with global envelope
        # Match parameters from sim_shade_comparison
        env <- envelope(
          pat,
          fun = Lcross,
          nsim = 99,
          nrank = 8,
          i = source,
          j = target,
          correction = "iso",
          r = seq(0, 75, by = 1),
          global = TRUE,
          savefuns = FALSE,
          verbose = FALSE,
          silent = TRUE
        )

        # Detection: does observed function exit envelope anywhere?
        # Only check in same distance range as SHADE (MIN_INTERACTION_RADIUS to 75)
        r_check <- env$r >= MIN_INTERACTION_RADIUS
        detected <- any((env$obs > env$hi | env$obs < env$lo)[r_check], na.rm = TRUE)

        data.frame(
          Spot = spot,
          Patient = patient,
          Group = group,
          detected = detected
        )

      }, error = \(e) {
        # Return NA if error, just note the pattern index
        cat("  Error in pattern", i, "\n")
        data.frame(
          Spot = spot,
          Patient = patient,
          Group = group,
          detected = NA
        )
      })

    }) %>%
      bind_rows()

    # Summarize detections
    # Image-level: proportion of images with detection
    image_summary <- image_results %>%
      summarise(
        n_images = n(),
        n_detected = sum(detected, na.rm = TRUE),
        prop_detected = n_detected / n_images,
        detected = prop_detected > 0.5  # Majority rule
      ) %>%
      mutate(
        target = target,
        source = source,
        level = "Image",
        Group = "All"
      )

    # Group-level: proportion within each group
    group_summary <- image_results %>%
      group_by(Group) %>%
      summarise(
        n_images = n(),
        n_detected = sum(detected, na.rm = TRUE),
        prop_detected = n_detected / n_images,
        detected = prop_detected > 0.5
      ) %>%
      mutate(
        target = target,
        source = source,
        level = "Group"
      )

    bind_rows(image_summary, group_summary)

  })

cat("✓ L-cross detections computed\n")
print(kcross_detections %>%
        mutate(method = "L-cross") %>%
        select(method, target, source, level, Group, detected, prop_detected))

# ============================================================================
# COMBINE ALL DETECTIONS
# ============================================================================

all_detections <- bind_rows(
  shade_detections %>% mutate(method = "SHADE", prop_detected = as.numeric(detected)),
  gcross_detections %>% mutate(method = "G-cross"),
  kcross_detections %>% mutate(method = "L-cross")
) %>%
  filter(level == "Group")  # Only keep group-level detections

# Save results
saveRDS(all_detections, paste0(path, "method_comparison_detections.rds"))
cat("\n✓ Saved all detections to:", paste0(path, "method_comparison_detections.rds"), "\n")

# ============================================================================
# VISUALIZATION: GROUP-LEVEL DETECTION COMPARISON
# ============================================================================

cat("\n=== Creating visualization: Group-level detection comparison ===\n")

group_comparison <- all_detections %>%
  select(method, target, source, Group, detected) %>%
  mutate(detected = ifelse(detected, "Yes", "No"))

group_comparison %>%
  ggplot(aes(x = source, y = target, fill = detected)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_grid(Group ~ method) +
  scale_fill_manual(values = c("No" = "gray90", "Yes" = "#D55E00"), name = "Detected") +
  labs(x = "Source cell type", y = "Target cell type",
       title = "Method Comparison by Patient Group (CLR vs DII)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA))

fsave("method_comparison_groups", height = 7, width = 12)

# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("\n=== Summary Statistics ===\n")

summary_stats <- all_detections %>%
  group_by(method, Group) %>%
  summarise(
    total_pairs = n(),
    n_detected = sum(detected, na.rm = TRUE),
    prop_detected = n_detected / total_pairs,
    .groups = "drop"
  )

print(summary_stats)

cat("\n=== Detection Agreement by Group ===\n")

# Compare SHADE vs G-cross vs L-cross at group level
for (grp in unique(all_detections$Group)) {
  cat("\nGroup:", grp, "\n")

  group_wide <- group_comparison %>%
    filter(Group == grp) %>%
    pivot_wider(names_from = method, values_from = detected, id_cols = c(target, source)) %>%
    mutate(pair = paste(source, "->", target))

  # Count agreement patterns
  agreement <- group_wide %>%
    count(SHADE, `G-cross`, `L-cross`) %>%
    arrange(desc(n))

  print(agreement)
}

cat("\n✓ Method comparison analysis complete!\n")
