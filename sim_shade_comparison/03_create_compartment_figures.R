library(tidyverse)
library(spatstat)
library(patchwork)

# Load utility functions
source("utils.R")
source("sim_shade_comparison/comparison_functions.R")
source("sim_shade_comparison/plot_confounder_diagnostics.R")

# Set theme
theme_set(theme_bw(base_size=14, base_family='Helvetica')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

# Figure saving helper
fsave <- \(fname,height=5,width=8) {
  ggsave(paste0(figures_folder,fname,".pdf"),device=cairo_pdf, height=height, width=width, units="in")
}
figures_folder <- "./manuscript/images/comparison_figures/"

# ============================================================================
# LOAD ANALYSIS SUMMARY
# ============================================================================

path <- get_data_path("sim_shade_comparison")
compartment_path <- paste0(path, "compartments/")
print(paste("Data path:", compartment_path))

analysis_file <- paste0(compartment_path, "analysis_summary.rds")
if(!file.exists(analysis_file)) {
  stop("Analysis summary not found. Run 02_analyze_compartments.R first.")
}

analysis_summary <- readRDS(analysis_file)
df <- analysis_summary$df
power_summary <- analysis_summary$power_summary
coverage_summary <- analysis_summary$coverage_summary
type_i_summary <- analysis_summary$type_i_summary
overall_type_i <- analysis_summary$overall_type_i

cat("Loaded compartment analysis summary with", nrow(df), "simulation results\n")

# ============================================================================
# COMPREHENSIVE SUMMARY TABLE
# ============================================================================

# Create comprehensive summary stratified by density and compartment effect
comprehensive_summary <- df %>%
  group_by(target_density, source_density, compartment_effect) %>%
  summarise(
    n = n(),
    across(c(shade_image_coverage, simple_shade_image_coverage,
             shade_image_power, simple_shade_image_power,
             gcross_image_power, kcross_image_power,
             type_i_error_shade, type_i_error_simple,
             type_i_error_gcross, type_i_error_kcross),
           list(median = ~median(., na.rm = TRUE),
                iqr = ~IQR(., na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = -c(target_density, source_density, compartment_effect, n),
               names_to = c("metric", "stat"),
               names_pattern = "(.*)_(median|iqr)$",
               values_to = "value") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    method = case_when(
      grepl("simple_shade", metric) ~ "SHADE Flat",
      grepl("shade_image", metric) ~ "SHADE Hierarchical",
      grepl("type_i_error_simple", metric) ~ "SHADE Flat",
      grepl("type_i_error_shade", metric) ~ "SHADE Hierarchical",
      grepl("gcross", metric) ~ "G-cross",
      grepl("kcross", metric) ~ "K-cross"
    ),
    measure = case_when(
      grepl("coverage", metric) ~ "Coverage",
      grepl("power", metric) ~ "Power",
      grepl("type_i", metric) ~ "Type I Error"
    )
  ) %>%
  select(target_density, source_density, compartment_effect, n, method, measure, median, iqr) %>%
  arrange(target_density, source_density, compartment_effect, measure, method) %>%
  mutate(across(c(median, iqr), ~round(.x, 3)))

cat("\n=== COMPREHENSIVE SUMMARY TABLE ===\n")
cat("Note: 'Source' refers to the conditioning cell type (T cells in simulation)\n")
cat("      'Target' refers to the cell type being modeled (tumor cells in simulation)\n")
cat("      'Compartment effect' strength: 0.8 (weak), 1.2 (moderate), 1.5 (strong)\n\n")
print(comprehensive_summary, n = Inf)

write_csv(comprehensive_summary, paste0(figures_folder, "compartment_comprehensive_summary.csv"))
cat("\n✓ Saved comprehensive summary table\n")

# ============================================================================
# MAIN MANUSCRIPT FIGURE: DIAGNOSTIC + POWER
# ============================================================================

cat("\n=== Creating combined figure for main manuscript ===\n")

# Generate example pattern for diagnostic visualization
# Use high-density scenario with moderate compartment effect as representative example
structure <- create_patient_structure(num_patients = 1, num_images = 1, seed = 2024)
coefficients <- generate_spatial_coefficients(
  structure = structure,
  num_potentials = 3,
  cohort_mean = c(1.5, 1.0, 0.5),
  seed = 2024
)
density_params <- list(np_t = 150, np_tumor = 150, np_b = 150)
patterns <- generate_spatial_patterns_compartments(
  structure = structure,
  coefficients = coefficients,
  density_params = density_params,
  n_compartments = 3,
  compartment_effect = 1.2,  # Moderate effect
  seed = 2024
)

# Create simplified diagnostic (just compartment structure + simulated pattern)
pattern_info <- patterns[[1]]
combined_pattern <- pattern_info$pattern
W <- Window(combined_pattern)

# Extract cell types
all_points <- as.data.frame(combined_pattern)
t_cells <- all_points %>% filter(marks == 1)
tumor_cells <- all_points %>% filter(marks == 3)

# Panel A: Compartment structure
compartment_image <- pattern_info$compartment_image
compartment_effects <- pattern_info$compartment_effects

effect_mat <- matrix(log(compartment_effects[compartment_image$v]),
                    nrow = nrow(compartment_image$v),
                    ncol = ncol(compartment_image$v))
confounder_density <- im(effect_mat,
                        xcol = compartment_image$xcol,
                        yrow = compartment_image$yrow)
confounder_df <- as.data.frame(confounder_density)

# Calculate centroids for each compartment for text annotations
compartment_mat <- as.matrix(compartment_image)
compartment_centroids <- data.frame()
for(comp_id in unique(as.vector(compartment_mat))) {
  coords <- which(compartment_mat == comp_id, arr.ind = TRUE)
  x_idx <- mean(coords[, "col"])
  y_idx <- mean(coords[, "row"])
  # Convert indices to actual coordinates
  x_coord <- compartment_image$xcol[round(x_idx)]
  y_coord <- compartment_image$yrow[round(y_idx)]
  effect_val <- log(compartment_effects[comp_id])
  compartment_centroids <- rbind(compartment_centroids,
                                 data.frame(x = x_coord, y = y_coord,
                                           label = sprintf("%.1f", effect_val)))
}

pA_compartment <- ggplot(confounder_df, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "Log\neffect") +
  geom_text(data = compartment_centroids, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 6, fontface = "bold", color = "white") +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(0, 0, 0, 0),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  guides(fill = "none")

pA_compartment

pA_pattern <- ggplot() +
  geom_point(data = t_cells, aes(x = x, y = y),
             color = "blue", size = 1.2, alpha = 0.5) +
  geom_point(data = tumor_cells, aes(x = x, y = y),
             color = "red", size = 1.2, alpha = 0.5) +
  coord_fixed(ratio = 1, xlim = c(0, max(W$xrange)), ylim = c(0, max(W$yrange)), expand = FALSE) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

pA <- pA_compartment / pA_pattern
# Panel C: Power results
pC <- df %>%
  mutate(compartment_effect = factor(compartment_effect,
                                    levels = c("0.8", "1.2", "1.5"),
                                    labels = c("Weak", "Moderate", "Strong"))) %>%
  rename(Gcross=gcross_image_power, Kcross=kcross_image_power,
         SHADE=shade_image_power, Flat=simple_shade_image_power) %>%
  pivot_longer(c(SHADE, Flat, Gcross, Kcross), names_to = "Method") %>%
  mutate(Method = factor(Method, levels = c("Gcross", "Kcross", "Flat", "SHADE"))) %>%
  rename(`Target density` = target_density, `Source density` = source_density) %>%
  ggplot(aes(x = compartment_effect, y = value, fill = Method, color = Method)) +
  geom_boxplot(alpha = 0.7) +
  facet_grid(`Source density` ~ `Target density`, labeller = label_both) +
  labs(x = "Compartment effect strength",
       y = "Detection power") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1))
pC
# Combine panels with A, B, C labeling
# Use plot_layout to control relative widths - give more space to left panels
combined_main <- pA | pC
combined_main <- combined_main +
  # plot_layout(widths = c(2.5, 2.5)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 16),
        plot.margin = margin(0, 0, 0, 0))
combined_main
# Save for main manuscript
ggsave(paste0(figures_folder, "compartment_main_figure.pdf"),
       plot = combined_main,
       device = cairo_pdf,
       height = 4.5,
       width = 7.5,
       units = "in",
       bg = "white")

cat("✓ Created combined figure: compartment_main_figure.pdf\n")

# ============================================================================
# SUPPLEMENT FIGURES
# ============================================================================

# Figure 1: Power by compartment effect strength
df %>%
  mutate(compartment_effect = factor(compartment_effect)) %>%
  rename(Gcross=gcross_image_power,Kcross=kcross_image_power,SHADE=shade_image_power,Flat=simple_shade_image_power) %>%
  pivot_longer(c(SHADE,Flat,Gcross,Kcross),names_to="Method") %>%
  mutate(Method=factor(Method,levels=c("Gcross","Kcross","Flat","SHADE"))) %>%
  rename(`Target density`=target_density,`Source density`=source_density) %>%
  ggplot(aes(x=compartment_effect,y=value,fill=Method,color=Method)) +
  geom_boxplot() +
  facet_grid(`Source density` ~ `Target density`,labeller=label_both) +
  labs(x="Compartment effect strength",y="Detection power (proportion)",
       title="Power: Can methods detect true interaction despite compartment confounder?")
fsave("compartment_power_by_effect")

# Figure 2: Coverage by compartment effect (SHADE only)
df %>%
  select(source_density, target_density, compartment_effect, shade_image_coverage, simple_shade_image_coverage) %>%
  pivot_longer(c(shade_image_coverage, simple_shade_image_coverage),
               names_to = "Method", values_to = "coverage") %>%
  mutate(
    Method = recode(Method,
                   shade_image_coverage = "SHADE Hierarchical",
                   simple_shade_image_coverage = "SHADE Flat"),
    compartment_effect = factor(compartment_effect)
  ) %>%
  rename(`Target density`=target_density,`Source density`=source_density) %>%
  ggplot(aes(x=compartment_effect,y=coverage,fill=Method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.95, linetype="dashed", color="red") +
  facet_grid(`Source density` ~ `Target density`,labeller=label_both) +
  labs(x="Compartment effect strength",y="Coverage (proportion)",
       title="95% Credible Band Coverage (Expected: 0.95)")
fsave("compartment_coverage_by_effect")

# Figure 3: Type I error rate by compartment effect
# This is the KEY robustness test: false positive rate when only compartment exists
df %>%
  mutate(compartment_effect = factor(compartment_effect)) %>%
  rename(Gcross=type_i_error_gcross,Kcross=type_i_error_kcross,SHADE=type_i_error_shade,Flat=type_i_error_simple) %>%
  pivot_longer(c(SHADE,Flat,Gcross,Kcross),names_to="Method") %>%
  mutate(Method=factor(Method,levels=c("Gcross","Kcross","Flat","SHADE"))) %>%
  rename(`Target density`=target_density,`Source density`=source_density) %>%
  ggplot(aes(x=compartment_effect,y=value,fill=Method,color=Method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05, linetype="dashed", color="red") +
  facet_grid(`Source density` ~ `Target density`,labeller=label_both) +
  labs(x="Compartment effect strength",y="Type I error rate (proportion)",
       title="False Positive Rate with Compartment Confounder (Expected: 0.05)",
       subtitle="Null pattern: Compartment effect exists but NO source-target interaction")
fsave("compartment_type_i_error_by_effect")

# Figure 4: Overall Type I error by compartment effect (aggregated across densities)
overall_type_i %>%
  pivot_longer(cols = -compartment_effect,
               names_to = c("Method", "stat"),
               names_pattern = "(.*)_(median|iqr)$",
               values_to = "value") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    Method = recode(Method,
                   shade = "SHADE Hierarchical",
                   flat = "SHADE Flat",
                   gcross = "G-cross",
                   kcross = "K-cross"),
    Method = factor(Method, levels = c("G-cross", "K-cross", "SHADE Flat", "SHADE Hierarchical")),
    compartment_effect = factor(compartment_effect)
  ) %>%
  ggplot(aes(x=compartment_effect, y=median, color=Method, group=Method)) +
  geom_line(linewidth=1) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=pmax(0, median-iqr/2), ymax=pmin(100, median+iqr/2)), width=0.1) +
  geom_hline(yintercept = 5, linetype="dashed", color="red", alpha=0.5) +
  labs(x="Compartment effect strength",
       y="Type I error rate (%)",
       title="Overall False Positive Rate by Compartment Effect",
       subtitle="Aggregated across all density regimes") +
  theme(legend.position = "bottom")
fsave("compartment_overall_type_i_error", height=5, width=7)

cat("\n✓ Figures saved to", figures_folder, "\n")
