library(tidyverse)

# Load utility functions
source("utils.R")

# Get data path (handles HPC vs local automatically)
path <- get_data_path("sim_shade_comparison")
compartment_path <- paste0(path, "compartments/")
print(paste("Data path:", compartment_path))

# ============================================================================
# LOAD SIMULATION RESULTS
# ============================================================================

# COMPARTMENT GRID
grid <- expand.grid(
  method = c("variational"),
  parameterization = c("centered"),
  scale_sigma = c(5),
  t_density = c("high", "low"),
  tumor_density = c("high", "low"),
  num_images = 3,  # Stable regime from main results
  n_compartments = 3,  # Moderate complexity
  compartment_effect = c(0.8, 1.2, 1.5),  # Weak, moderate, strong confounding
  sim_rep = 1:100
) %>%
  mutate(sim_idx = 1:nrow(.))

cat("Loading", nrow(grid), "compartment confounder simulation results...\n")

df <- pmap(grid, function(...) {
  row <- list(...)
  sim_idx <- row$sim_idx
  file_out <- paste0(compartment_path, "sim_", sim_idx, ".rds")

  out <- tryCatch(readRDS(file_out), error = function(e) NULL)
  if(is.null(out)) return(NULL)

  tibble(
    # Power metrics
    shade_image_power = out$shade_image_power,
    simple_shade_image_power = out$simple_shade_image_power,
    gcross_image_power = out$gcross_image_power,
    kcross_image_power = out$kcross_image_power,
    # Coverage metrics
    shade_image_coverage = out$shade_image_coverage %||% NA,
    simple_shade_image_coverage = out$simple_shade_image_coverage %||% NA,
    # Type I error metrics (with compartment confounder)
    type_i_error_shade = out$type_i_error_shade %||% NA,
    type_i_error_simple = out$type_i_error_simple %||% NA,
    type_i_error_gcross = out$type_i_error_gcross %||% NA,
    type_i_error_kcross = out$type_i_error_kcross %||% NA,
    # Extract condition from saved results
    method = out$condition$method,
    parameterization = out$condition$parameterization,
    scale_sigma = out$condition$scale_sigma,
    t_density = out$condition$t_density,
    tumor_density = out$condition$tumor_density,
    num_images = as.numeric(out$condition$num_images),
    n_compartments = out$condition$n_compartments,
    compartment_effect = out$condition$compartment_effect,
    sim_rep = out$condition$sim_rep,
    sim_idx = sim_idx
  )
}) %>%
  bind_rows()

cat("Loaded", nrow(df), "compartment confounder simulation results\n")

# ============================================================================
# SUMMARY TABLES
# ============================================================================

# Power summary (grouped by density and compartment effect)
power_summary <- df %>%
  group_by(t_density, tumor_density, compartment_effect) %>%
  summarise(
    n = n(),
    shade_power_median = median(shade_image_power, na.rm=TRUE),
    shade_power_iqr = IQR(shade_image_power, na.rm=TRUE),
    flat_power_median = median(simple_shade_image_power, na.rm=TRUE),
    flat_power_iqr = IQR(simple_shade_image_power, na.rm=TRUE),
    gcross_power_median = median(gcross_image_power, na.rm=TRUE),
    gcross_power_iqr = IQR(gcross_image_power, na.rm=TRUE),
    kcross_power_median = median(kcross_image_power, na.rm=TRUE),
    kcross_power_iqr = IQR(kcross_image_power, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_median"), ~round(.x * 100, 1))) %>%
  mutate(across(ends_with("_iqr"), ~round(.x * 100, 1)))

cat("\n=== DETECTION POWER (%) BY COMPARTMENT EFFECT ===\n")
print(power_summary, n = Inf)

# Coverage summary (SHADE only)
coverage_summary <- df %>%
  group_by(t_density, tumor_density, compartment_effect) %>%
  summarise(
    n = n(),
    shade_coverage_median = median(shade_image_coverage, na.rm=TRUE),
    shade_coverage_iqr = IQR(shade_image_coverage, na.rm=TRUE),
    flat_coverage_median = median(simple_shade_image_coverage, na.rm=TRUE),
    flat_coverage_iqr = IQR(simple_shade_image_coverage, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_median"), ~round(.x * 100, 1))) %>%
  mutate(across(ends_with("_iqr"), ~round(.x * 100, 1)))

cat("\n=== COVERAGE (%) BY COMPARTMENT EFFECT - Expected: ~95% ===\n")
print(coverage_summary, n = Inf)

# ============================================================================
# TYPE I ERROR SUMMARY
# ============================================================================

# Type I error summary by condition
# This tests: Does SHADE falsely detect interaction when only compartment exists?
type_i_summary <- df %>%
  group_by(t_density, tumor_density, compartment_effect) %>%
  summarise(
    n = n(),
    shade_type_i_median = median(type_i_error_shade, na.rm=TRUE),
    shade_type_i_iqr = IQR(type_i_error_shade, na.rm=TRUE),
    flat_type_i_median = median(type_i_error_simple, na.rm=TRUE),
    flat_type_i_iqr = IQR(type_i_error_simple, na.rm=TRUE),
    gcross_type_i_median = median(type_i_error_gcross, na.rm=TRUE),
    gcross_type_i_iqr = IQR(type_i_error_gcross, na.rm=TRUE),
    kcross_type_i_median = median(type_i_error_kcross, na.rm=TRUE),
    kcross_type_i_iqr = IQR(type_i_error_kcross, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_median"), ~round(.x * 100, 1))) %>%
  mutate(across(ends_with("_iqr"), ~round(.x * 100, 1)))

cat("\n=== TYPE I ERROR RATE (%) WITH COMPARTMENT CONFOUNDER - Expected: ~5% ===\n")
cat("Tests: False detection rate when compartment effect exists but NO source-target interaction\n")
print(type_i_summary, n = Inf)

# Overall Type I error (aggregated across densities)
overall_type_i <- df %>%
  group_by(compartment_effect) %>%
  summarise(
    shade_median = median(type_i_error_shade, na.rm=TRUE),
    shade_iqr = IQR(type_i_error_shade, na.rm=TRUE),
    flat_median = median(type_i_error_simple, na.rm=TRUE),
    flat_iqr = IQR(type_i_error_simple, na.rm=TRUE),
    gcross_median = median(type_i_error_gcross, na.rm=TRUE),
    gcross_iqr = IQR(type_i_error_gcross, na.rm=TRUE),
    kcross_median = median(type_i_error_kcross, na.rm=TRUE),
    kcross_iqr = IQR(type_i_error_kcross, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_median"), ~round(.x * 100, 1))) %>%
  mutate(across(ends_with("_iqr"), ~round(.x * 100, 1)))

cat("\n=== OVERALL TYPE I ERROR BY COMPARTMENT EFFECT STRENGTH ===\n")
print(overall_type_i, n = Inf)

# ============================================================================
# RENAME FOR CONSISTENCY WITH MAIN ANALYSIS
# ============================================================================
# Use source/target terminology instead of t_cell/tumor

df_renamed <- df %>%
  rename(
    source_density = t_density,
    target_density = tumor_density
  )

power_summary_renamed <- power_summary %>%
  rename(
    source_density = t_density,
    target_density = tumor_density
  )

coverage_summary_renamed <- coverage_summary %>%
  rename(
    source_density = t_density,
    target_density = tumor_density
  )

type_i_summary_renamed <- type_i_summary %>%
  rename(
    source_density = t_density,
    target_density = tumor_density
  )

# ============================================================================
# SAVE ANALYSIS SUMMARY
# ============================================================================

analysis_summary <- list(
  df = df_renamed,
  grid = grid,
  power_summary = power_summary_renamed,
  coverage_summary = coverage_summary_renamed,
  type_i_summary = type_i_summary_renamed,
  overall_type_i = overall_type_i
)

saveRDS(analysis_summary, paste0(compartment_path, "analysis_summary.rds"))
cat("\n✓ Compartment analysis summary saved to", paste0(compartment_path, "analysis_summary.rds"), "\n")
