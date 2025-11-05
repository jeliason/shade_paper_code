library(tidyverse)
library(spatstat)
library(Matrix)
library(SHADE)
library(cmdstanr)
library(posterior)
library(patchwork)

# Load utility functions and constants
source("./utils.R")
source("sim_shade_comparison/comparison_functions.R")

# Set theme
theme_set(theme_bw(base_size=14, base_family='Helvetica')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

# ============================================================================
# CONFIGURATION
# ============================================================================

# Set random seed for reproducibility
set.seed(2024)

# Output directory
fig_dir <- "manuscript/images/comparison_figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Simulation parameters (matching your current grid)
num_patients <- 10
num_images_per_patient <- 2
num_potentials <- 3
size_im <- 1500

# Use high density for a good example
density_params <- list(
  np_t = 50,      # High T cell density
  np_tumor = 50,  # High tumor density
  np_b = 150       # High B cell density
)

# Model configuration
method <- "variational"
parameterization <- "centered"
scale_sigma <- 5

# ============================================================================
# GENERATE EXAMPLE DATA
# ============================================================================

cat("Creating patient structure...\n")
structure <- create_patient_structure(
  num_patients = num_patients,
  num_images = num_images_per_patient,
  seed = 2024
)

cat("Generating spatial coefficients...\n")
coefficients <- generate_spatial_coefficients(
  structure = structure,
  num_potentials = num_potentials,
  cohort_mean = c(1.5, 1.0, 0.5),
  seed = 2024
)

cat("Generating spatial patterns...\n")
patterns <- generate_spatial_patterns(
  structure = structure,
  coefficients = coefficients,
  density_params = density_params,
  size_im = size_im,
  seed = 2024 + 1000
)

# Select first pattern as example
example_idx <- 1
example_pattern <- patterns[[example_idx]]

cat(sprintf("\nExample image %d (Patient %d):\n",
            example_idx, example_pattern$patient_id))
cat("  T cells:", sum(example_pattern$pattern$marks == "1"), "\n")
cat("  B cells:", sum(example_pattern$pattern$marks == "2"), "\n")
cat("  Tumor cells:", sum(example_pattern$pattern$marks == "3"), "\n")
cat("  True coefficients:", paste(round(example_pattern$image_coeffs, 3), collapse = ", "), "\n")

# ============================================================================
# RUN ANALYSES ON EXAMPLE IMAGE
# ============================================================================

# Define potentials
potentials <- make_rbfs(max_dist = 75, n_basis_functions = 3, basis_function_sigma = 15)

# Run SHADE hierarchical (on all images to get hierarchical structure)
cat("\nRunning SHADE hierarchical analysis...\n")
shade_results <- run_shade_analysis(
  patterns,
  structure,
  num_potentials,
  method = method,
  parameterization = parameterization,
  scale_sigma = scale_sigma
)
shade_draws <- as_draws_rvars(shade_results$fit$draws())

# Run simple SHADE on example image only
cat("Running SHADE flat analysis on example...\n")
mod <- cmdstan_model("sim_shade_comparison/super_simple_shade.stan", quiet = TRUE)
simple_shade_result <- simple_shade(example_pattern$pattern, potentials, mod)
simple_shade_draws <- as_draws_rvars(simple_shade_result$fit$draws())

# Run G-cross on example image
cat("Running G-cross on example...\n")
pat <- example_pattern$pattern
r_seq <- seq(DISTANCE_RANGE_MIN, DISTANCE_RANGE_MAX, by = 1)
gcross_env <- envelope(pat, Gcross, i = "1", j = "3",
                      nsim = 99, nrank = 10,
                      correction = "km", global = TRUE,
                      r = r_seq, verbose = FALSE, savefuns = TRUE)

# Run K-cross (L-cross) on example image
cat("Running K-cross on example...\n")
kcross_env <- envelope(pat, Lcross, i = "1", j = "3",
                      nsim = 99, nrank = 8,
                      correction = "iso", global = TRUE,
                      r = r_seq, verbose = FALSE, savefuns = TRUE)

# ============================================================================
# EXTRACT SIC ESTIMATES
# ============================================================================

patient_id <- example_pattern$patient_id
num_images_for_patient <- sum(structure$images_df$patient_id == patient_id)

# For multi-image patient: use image-level estimates
if(num_images_for_patient > 1) {
  beta_hier <- as.vector(shade_draws$beta_local[example_idx, 2:4])
  estimate_level <- "Image"
} else {
  beta_hier <- as.vector(shade_draws$beta_indiv[2:4, patient_id])
  estimate_level <- "Patient"
}

# Simple SHADE: image-level only
beta_flat <- as.vector(simple_shade_draws$beta[1:3])

# True coefficients
beta_true <- example_pattern$image_coeffs

# Compute SICs
x_seq <- seq(DISTANCE_RANGE_MIN, DISTANCE_RANGE_MAX, 1)
x_des <- lapply(potentials, \(pot) pot(x_seq)) %>% do.call(cbind, .)

# Hierarchical SHADE SIC
lp_hier <- as.vector(x_des %*% beta_hier)
lp_data_hier <- lp_hier %>% as.matrix() %>% as.data.frame()
colnames(lp_data_hier) <- "SIC"
df_hier <- compute_simultaneous_bands(lp_data_hier, x_seq, alpha = 0.05) %>%
  rename(mn = mean, lo = lower, hi = upper)

# Flat SHADE SIC
lp_flat <- as.vector(x_des %*% beta_flat)
lp_data_flat <- lp_flat %>% as.matrix() %>% as.data.frame()
colnames(lp_data_flat) <- "SIC"
df_flat <- compute_simultaneous_bands(lp_data_flat, x_seq, alpha = 0.05) %>%
  rename(mn = mean, lo = lower, hi = upper)

# True SIC
lp_true <- as.vector(x_des %*% beta_true)

# ============================================================================
# CREATE FIGURE PANELS
# ============================================================================

# Panel A: Simulated pattern
p_pattern <- ggplot() +
  geom_point(data = data.frame(x = pat$x, y = pat$y, type = marks(pat)),
            aes(x = x, y = y, color = type, shape = type),
            size = 2, alpha = 0.6) +
  scale_color_manual(
    values = c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A"),
    labels = c("1" = "T cells", "2" = "B cells", "3" = "Tumor cells"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("1" = 16, "2" = 17, "3" = 15),
    labels = c("1" = "T cells", "2" = "B cells", "3" = "Tumor cells"),
    name = NULL
  ) +
  coord_fixed() +
  labs(
    title = "(a) Simulated Pattern",
    x = expression(paste("x (", mu, "m)")),
    y = expression(paste("y (", mu, "m)"))
  ) +
  theme(
    legend.position = c(0.85, 0.2),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 2)))
p_pattern
# Panel B: G-cross
gcross_df <- data.frame(
  r = gcross_env$r,
  obs = gcross_env$obs,
  theo = gcross_env$theo,
  lo = gcross_env$lo,
  hi = gcross_env$hi
)

p_gcross <- ggplot(gcross_df, aes(x = r)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "gray80", alpha = 0.5) +
  geom_line(aes(y = theo, color = "CSR", linetype = "CSR"), size = 0.8) +
  geom_line(aes(y = obs, color = "Observed", linetype = "Observed"), size = 1) +
  scale_color_manual(
    name = NULL,
    values = c("Observed" = "black", "CSR" = "red"),
    breaks = c("Observed", "CSR")
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("Observed" = "solid", "CSR" = "dashed"),
    breaks = c("Observed", "CSR")
  ) +
  labs(
    title = "(b) G-cross Analysis",
    x = expression(paste("Distance r (", mu, "m)")),
    y = "G-cross(r)"
  ) +
  theme(
    legend.position = c(0.2, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.title = element_blank()
  )
p_gcross
# Panel C: K-cross (L-cross)
kcross_df <- data.frame(
  r = kcross_env$r,
  obs = kcross_env$obs,
  theo = kcross_env$theo,
  lo = kcross_env$lo,
  hi = kcross_env$hi
)

p_kcross <- ggplot(kcross_df, aes(x = r)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "gray80", alpha = 0.5) +
  geom_line(aes(y = theo, color = "CSR", linetype = "CSR"), size = 0.8) +
  geom_line(aes(y = obs, color = "Observed", linetype = "Observed"), size = 1) +
  scale_color_manual(
    name = NULL,
    values = c("Observed" = "black", "CSR" = "red"),
    breaks = c("Observed", "CSR")
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("Observed" = "solid", "CSR" = "dashed"),
    breaks = c("Observed", "CSR")
  ) +
  labs(
    title = "(c) L-cross Analysis",
    x = expression(paste("Distance r (", mu, "m)")),
    y = "L-cross(r)"
  ) +
  theme(
    legend.position = c(0.2, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.title = element_blank()
  )
p_kcross
# Panel D: SIC comparison
sic_df <- data.frame(
  x = x_seq,
  true = lp_true,
  hier_mn = df_hier$mn,
  hier_lo = df_hier$lo,
  hier_hi = df_hier$hi,
  flat_mn = df_flat$mn,
  flat_lo = df_flat$lo,
  flat_hi = df_flat$hi
)

p_sic <- ggplot(sic_df, aes(x = x)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50", size = 0.3) +
  # Hierarchical SHADE
  geom_ribbon(aes(ymin = hier_lo, ymax = hier_hi, fill = "SHADE Hierarchical"),
              alpha = 0.3) +
  geom_line(aes(y = hier_mn, color = "SHADE Hierarchical", linetype = "SHADE Hierarchical"), size = 1) +
  # Flat SHADE
  geom_ribbon(aes(ymin = flat_lo, ymax = flat_hi, fill = "SHADE Flat"),
              alpha = 0.3) +
  geom_line(aes(y = flat_mn, color = "SHADE Flat", linetype = "SHADE Flat"), size = 1) +
  # Ground truth
  geom_line(aes(y = true, color = "Ground Truth", linetype = "Ground Truth"), size = 1) +
  scale_color_manual(
    name = NULL,
    values = c(
      "SHADE Hierarchical" = "#E41A1C",
      "SHADE Flat" = "#377EB8",
      "Ground Truth" = "black"
    ),
    breaks = c("SHADE Hierarchical", "SHADE Flat", "Ground Truth")
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "SHADE Hierarchical" = "solid",
      "SHADE Flat" = "solid",
      "Ground Truth" = "dashed"
    ),
    breaks = c("SHADE Hierarchical", "SHADE Flat", "Ground Truth")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "SHADE Hierarchical" = "#E41A1C",
      "SHADE Flat" = "#377EB8"
    ),
    breaks = c("SHADE Hierarchical", "SHADE Flat"),
    guide = "none"
  ) +
  labs(
    title = "(d) Spatial Interaction Curves",
    x = expression(paste("Distance r (", mu, "m)")),
    y = "SIC"
  ) +
  theme(
    legend.position = c(0.82, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.title = element_blank()
  )
p_sic
# ============================================================================
# COMBINE PANELS
# ============================================================================

# 2x2 layout
p_combined <- (p_pattern + p_gcross) / (p_kcross + p_sic)
p_combined
# Save
ggsave(
  paste0(fig_dir, "/example_comparison_2x2.pdf"),
  p_combined,
  width = 12,
  height = 10
)

cat("\n✓ Saved 2x2 comparison figure\n")

# Also save 1x4 horizontal version
p_horizontal <- p_pattern + p_gcross + p_kcross + p_sic +
  plot_layout(nrow = 1)

ggsave(
  paste0(fig_dir, "/example_comparison_1x4.pdf"),
  p_horizontal,
  width = 16,
  height = 4.5
)

cat("✓ Saved 1x4 horizontal comparison figure\n")

# ============================================================================
# DETECTION STATISTICS
# ============================================================================

cat("\n=== Detection Statistics for Example Image ===\n")

# SHADE hierarchical
excludes_zero_hier <- any(df_hier$lo > 0 | df_hier$hi < 0, na.rm = TRUE)
cat("SHADE Hierarchical detection:", excludes_zero_hier, "\n")

# SHADE flat
excludes_zero_flat <- any(df_flat$lo > 0 | df_flat$hi < 0, na.rm = TRUE)
cat("SHADE Flat detection:", excludes_zero_flat, "\n")

# G-cross
gcross_detect <- any(gcross_env$obs > gcross_env$hi | gcross_env$obs < gcross_env$lo, na.rm = TRUE)
cat("G-cross detection:", gcross_detect, "\n")

# K-cross
kcross_detect <- any(kcross_env$obs > kcross_env$hi | kcross_env$obs < kcross_env$lo, na.rm = TRUE)
cat("K-cross detection:", kcross_detect, "\n")

cat("\n=== Figure saved to:", fig_dir, "===\n")
