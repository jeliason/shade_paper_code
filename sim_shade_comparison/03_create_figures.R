library(tidyverse)

# Load utility functions
source("utils.R")

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
print(paste("Data path:", path))

analysis_file <- paste0(path, "analysis_summary.rds")
if(!file.exists(analysis_file)) {
  stop("Analysis summary not found. Run 03_analyze_results.R first.")
}

analysis_summary <- readRDS(analysis_file)
df <- analysis_summary$df
power_summary <- analysis_summary$power_summary
coverage_summary <- analysis_summary$coverage_summary
type_i_summary <- analysis_summary$type_i_summary
overall_type_i <- analysis_summary$overall_type_i

cat("Loaded analysis summary with", nrow(df), "simulation results (power + null calibration)\n")

# ============================================================================
# COMPREHENSIVE SUMMARY TABLE
# ============================================================================

# Create comprehensive summary stratified by all dimensions
comprehensive_summary <- df %>%
  group_by(tumor_density, t_density, num_images) %>%
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
  pivot_longer(cols = -c(tumor_density, t_density, num_images, n),
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
  select(tumor_density, t_density, num_images, n, method, measure, median, iqr) %>%
  arrange(tumor_density, t_density, num_images, measure, method) %>%
  mutate(across(c(median, iqr), ~round(.x, 3))) %>%
  rename(target_density = tumor_density, source_density = t_density)

cat("\n=== COMPREHENSIVE SUMMARY TABLE ===\n")
cat("Note: 'Source' refers to the conditioning cell type (T cells in simulation)\n")
cat("      'Target' refers to the cell type being modeled (tumor cells in simulation)\n\n")
print(comprehensive_summary, n = Inf)

write_csv(comprehensive_summary, paste0(figures_folder, "comprehensive_summary.csv"))
cat("\n✓ Saved comprehensive summary table\n")

# ============================================================================
# FIGURES
# ============================================================================

# Figure 1: Power by data regime
df %>%
  mutate(num_images = factor(num_images)) %>%
  rename(Gcross=gcross_image_power,Kcross=kcross_image_power,SHADE=shade_image_power,Flat=simple_shade_image_power) %>%
  pivot_longer(c(SHADE,Flat,Gcross,Kcross),names_to="Method") %>%
  mutate(Method=factor(Method,levels=c("Gcross","Kcross","Flat","SHADE"))) %>%
  rename(`Target density`=tumor_density,`Source density`=t_density) %>%
  ggplot(aes(x=num_images,y=value,fill=Method,color=Method)) +
  geom_boxplot() +
  facet_grid(`Source density` ~ `Target density`,labeller=label_both) +
  labs(x="Number of images per patient",y="Detection power (proportion)")
fsave("power_by_regime")

# Figure 2: Coverage by data regime (SHADE only)
df %>%
  select(t_density, tumor_density, num_images, shade_image_coverage, simple_shade_image_coverage) %>%
  pivot_longer(c(shade_image_coverage, simple_shade_image_coverage),
               names_to = "Method", values_to = "coverage") %>%
  mutate(
    Method = recode(Method,
                   shade_image_coverage = "SHADE Hierarchical",
                   simple_shade_image_coverage = "SHADE Flat"),
    num_images = factor(num_images)
  ) %>%
  rename(`Target density`=tumor_density,`Source density`=t_density) %>%
  ggplot(aes(x=num_images,y=coverage,fill=Method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.95, linetype="dashed", color="red") +
  facet_grid(`Source density` ~ `Target density`,labeller=label_both) +
  labs(x="Number of images per patient",y="Coverage (proportion)",
       title="95% Credible Band Coverage (Expected: 0.95)")
fsave("coverage_by_regime")

# Figure 3: Type I error rate by data regime
df %>%
  mutate(num_images = factor(num_images)) %>%
  rename(Gcross=type_i_error_gcross,Kcross=type_i_error_kcross,SHADE=type_i_error_shade,Flat=type_i_error_simple) %>%
  pivot_longer(c(SHADE,Flat,Gcross,Kcross),names_to="Method") %>%
  mutate(Method=factor(Method,levels=c("Gcross","Kcross","Flat","SHADE"))) %>%
  rename(`Target density`=tumor_density,`Source density`=t_density) %>%
  ggplot(aes(x=num_images,y=value,fill=Method,color=Method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05, linetype="dashed", color="red") +
  facet_grid(`Source density` ~ `Target density`,labeller=label_both) +
  labs(x="Number of images per patient",y="Type I error rate (proportion)",
       title="False Positive Rate (Expected: 0.05)")
fsave("type_i_error_by_regime")

cat("\n✓ Figures saved to", figures_folder, "\n")
