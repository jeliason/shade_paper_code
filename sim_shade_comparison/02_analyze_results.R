library(tidyverse)

# Load utility functions
source("utils.R")

# Get data path (handles HPC vs local automatically)
path <- get_data_path("sim_shade_comparison")
print(paste("Data path:", path))

# ============================================================================
# LOAD SIMULATION RESULTS
# ============================================================================

# ORIGINAL GRID
# grid <- expand.grid(
#   t_density = c("high", "low"),
#   tumor_density = c("high", "low"),
#   num_images = c(1,2,3),
#   sim_rep = 1:30
# ) %>%
#   mutate(sim_idx = 1:nrow(.))

# EXPERIMENTAL GRID
grid <- expand.grid(
  method = c("variational"),# "sampling"),
  parameterization = c("centered"),
  scale_sigma = c(5),
  t_density = c("high", "low"),        # Test both to see if data quantity matters
  tumor_density = c("high","low"),              # Keep tumor density high
  num_images = c(1,2,3),                      # Best case: 3 images per patient
  sim_rep = 1:50
) %>%
  mutate(sim_idx = 1:nrow(.))

cat("Loading", nrow(grid), "simulation results...\n")

df <- pmap(grid, function(...) {
  row <- list(...)
  sim_idx <- row$sim_idx
  file_out <- paste0(path, "sim_", sim_idx, ".rds")

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
    # Type I error metrics
    type_i_error_shade = out$type_i_error_shade %||% NA,
    type_i_error_simple = out$type_i_error_simple %||% NA,
    type_i_error_gcross = out$type_i_error_gcross %||% NA,
    type_i_error_kcross = out$type_i_error_kcross %||% NA,
    # Extract condition from saved results (handles both grid types)
    method = out$condition$method %||% "variational",
    parameterization = out$condition$parameterization %||% "centered",
    scale_sigma = out$condition$scale_sigma %||% 5,
    t_density = out$condition$t_density,
    tumor_density = out$condition$tumor_density,
    num_images = as.numeric(out$condition$num_images),
    sim_rep = out$condition$sim_rep,
    sim_idx = sim_idx
  )
}) %>%
  bind_rows()

cat("Loaded", nrow(df), "simulation results (power + null calibration)\n")

# ============================================================================
# SUMMARY TABLES
# ============================================================================

# Power summary (grouped by experimental conditions)
power_summary <- df %>%
  group_by(method, parameterization, scale_sigma, t_density) %>%
  summarise(
    n = n(),
    shade_power_mean = mean(shade_image_power, na.rm=TRUE),
    flat_power_mean = mean(simple_shade_image_power, na.rm=TRUE),
    gcross_power_mean = mean(gcross_image_power, na.rm=TRUE),
    kcross_power_mean = mean(kcross_image_power, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_mean"), ~round(.x * 100, 1)))

cat("\n=== DETECTION POWER (%) ===\n")
print(power_summary, n = Inf)

# Coverage summary (SHADE only)
coverage_summary <- df %>%
  group_by(method, parameterization, scale_sigma, t_density) %>%
  summarise(
    n = n(),
    shade_coverage_mean = mean(shade_image_coverage, na.rm=TRUE),
    flat_coverage_mean = mean(simple_shade_image_coverage, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_mean"), ~round(.x * 100, 1)))

cat("\n=== COVERAGE (%) - Expected: ~95% ===\n")
print(coverage_summary, n = Inf)

# ============================================================================
# TYPE I ERROR SUMMARY
# ============================================================================

# Type I error summary by condition
type_i_summary <- df %>%
  group_by(method, parameterization, scale_sigma, t_density) %>%
  summarise(
    n = n(),
    shade_type_i_mean = mean(type_i_error_shade, na.rm=TRUE),
    flat_type_i_mean = mean(type_i_error_simple, na.rm=TRUE),
    gcross_type_i_mean = mean(type_i_error_gcross, na.rm=TRUE),
    kcross_type_i_mean = mean(type_i_error_kcross, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_mean"), ~round(.x * 100, 1)))

cat("\n=== TYPE I ERROR RATE BY CONDITION (%) - Expected: ~5% ===\n")
print(type_i_summary, n = Inf)

# Overall Type I error (aggregated across all conditions)
overall_type_i <- df %>%
  group_by(method, parameterization, scale_sigma) %>%
  summarise(
    shade = mean(type_i_error_shade, na.rm=TRUE),
    flat = mean(type_i_error_simple, na.rm=TRUE),
    gcross = mean(type_i_error_gcross, na.rm=TRUE),
    kcross = mean(type_i_error_kcross, na.rm=TRUE),
    .groups = "drop"
  )

cat("\n=== OVERALL TYPE I ERROR RATE (%) BY MODEL CONFIG ===\n")
print(overall_type_i %>% mutate(across(c(shade, flat, gcross, kcross), ~round(.x * 100, 1))), n = Inf)

# ============================================================================
# SAVE ANALYSIS SUMMARY
# ============================================================================

analysis_summary <- list(
  df = df,
  grid = grid,
  power_summary = power_summary,
  coverage_summary = coverage_summary,
  type_i_summary = type_i_summary,
  overall_type_i = overall_type_i
)

saveRDS(analysis_summary, paste0(path, "analysis_summary.rds"))
cat("\n✓ Analysis summary saved to", paste0(path, "analysis_summary.rds"), "\n")
