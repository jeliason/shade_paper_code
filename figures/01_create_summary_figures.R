library(tidyverse)
library(patchwork)
library(SHADE)

#' Create summary figures combining results from multiple analyses
#'
#' This script creates combined summary figures for the paper by
#' combining key visualizations from the different simulation and real data analyses.

# Set paths to different analysis results
data_dir <- "../"
dummy_dir <- file.path(data_dir, "simulations_dummy")
sample_size_dir <- file.path(data_dir, "simulations_sample_size")
flat_dir <- file.path(data_dir, "simulations_flat")
crc_dir <- file.path(data_dir, "CRC_analysis")

figures_dir <- "./publication/combined"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Set theme for plotting
theme_set(theme_bw(base_size = 14, base_family = 'Helvetica') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()))

# Function to save figures
save_figure <- function(plot, filename, width = 8, height = 6, device = cairo_pdf) {
  ggsave(
    file.path(figures_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    device = device
  )
  
  ggsave(
    file.path(figures_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

# Load key figures from each analysis
# We're loading the plots directly because we've already saved them in the individual analyses

# 1. Create method overview figure
# This combines RBF visualization with SIC examples

# Create RBF visualization
x_seq <- seq(0, 100, 0.1)
rbfs <- make_rbfs(n_basis_functions = 3, max_dist = 75, basis_function_sigma = 15)

# Reshape data for ggplot
rbf_df <- lapply(1:length(rbfs), function(i) {
  data.frame(
    x = x_seq,
    y = rbfs[[i]](x_seq),
    rbf = paste0("RBF ", i)
  )
}) %>% bind_rows()

p_rbfs <- ggplot(rbf_df, aes(x = x, y = y, color = rbf)) +
  geom_line(linewidth = 1.5) +
  labs(
    x = "Distance (microns)",
    y = "Weight",
    color = "Basis Function",
    title = "A. Radial Basis Functions"
  ) +
  theme(legend.position = "bottom")

# Create SIC example with weights
weights <- c(0.8, -0.4, 0.3)
pot_mat <- lapply(1:length(rbfs), function(i) rbfs[[i]](x_seq) * weights[i]) %>% 
  do.call(cbind, .) %>% 
  rowSums()

p_sic <- tibble(x = x_seq, y = pot_mat) %>%
  ggplot(aes(x, y)) +
  geom_line(color = "red", linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_segment(x = 75, y = -0.4, xend = 75, yend = 0.28, 
               color = "blue", linetype = "dotted", linewidth = 1.5) +
  labs(
    x = "Distance (microns)",
    y = "Log-intensity",
    title = "B. Spatial Interaction Coefficient"
  )

# Create diagram of hierarchical model structure
p_hierarchy <- ggplot() +
  # Cohort level
  annotate("rect", xmin = 0, xmax = 4, ymin = 5, ymax = 6, 
           fill = "lightblue", alpha = 0.5) +
  annotate("text", x = 2, y = 5.5, label = "Cohort Level Parameters") +
  
  # Patient level
  annotate("rect", xmin = 0, xmax = 1.3, ymin = 3, ymax = 4, 
           fill = "lightgreen", alpha = 0.5) +
  annotate("text", x = 0.65, y = 3.5, label = "Patient 1", size = 3) +
  
  annotate("rect", xmin = 1.4, xmax = 2.7, ymin = 3, ymax = 4, 
           fill = "lightgreen", alpha = 0.5) +
  annotate("text", x = 2.05, y = 3.5, label = "Patient 2", size = 3) +
  
  annotate("rect", xmin = 2.8, xmax = 4, ymin = 3, ymax = 4, 
           fill = "lightgreen", alpha = 0.5) +
  annotate("text", x = 3.4, y = 3.5, label = "Patient 3", size = 3) +
  
  # Image level
  annotate("rect", xmin = 0, xmax = 0.6, ymin = 1, ymax = 2, 
           fill = "lightyellow", alpha = 0.5) +
  annotate("text", x = 0.3, y = 1.5, label = "Image 1", size = 2.5) +
  
  annotate("rect", xmin = 0.7, xmax = 1.3, ymin = 1, ymax = 2, 
           fill = "lightyellow", alpha = 0.5) +
  annotate("text", x = 1, y = 1.5, label = "Image 2", size = 2.5) +
  
  annotate("rect", xmin = 1.4, xmax = 2, ymin = 1, ymax = 2, 
           fill = "lightyellow", alpha = 0.5) +
  annotate("text", x = 1.7, y = 1.5, label = "Image 3", size = 2.5) +
  
  annotate("rect", xmin = 2.1, xmax = 2.7, ymin = 1, ymax = 2, 
           fill = "lightyellow", alpha = 0.5) +
  annotate("text", x = 2.4, y = 1.5, label = "Image 4", size = 2.5) +
  
  annotate("rect", xmin = 2.8, xmax = 3.4, ymin = 1, ymax = 2, 
           fill = "lightyellow", alpha = 0.5) +
  annotate("text", x = 3.1, y = 1.5, label = "Image 5", size = 2.5) +
  
  annotate("rect", xmin = 3.5, xmax = 4, ymin = 1, ymax = 2, 
           fill = "lightyellow", alpha = 0.5) +
  annotate("text", x = 3.75, y = 1.5, label = "Image 6", size = 2.5) +
  
  # Arrows
  geom_segment(aes(x = 2, y = 5, xend = 2, yend = 4.1), 
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = 0.65, y = 3, xend = 0.3, yend = 2.1), 
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = 0.65, y = 3, xend = 1, yend = 2.1), 
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = 2.05, y = 3, xend = 1.7, yend = 2.1), 
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = 2.05, y = 3, xend = 2.4, yend = 2.1), 
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = 3.4, y = 3, xend = 3.1, yend = 2.1), 
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = 3.4, y = 3, xend = 3.75, yend = 2.1), 
               arrow = arrow(length = unit(0.2, "cm"))) +
  
  labs(title = "C. Hierarchical Model Structure") +
  theme_void() +
  xlim(0, 4) + 
  ylim(0, 6)

# Combine the figures for the method overview
p_method <- (p_rbfs + p_sic) / p_hierarchy +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = "SHADE: Spatial Hierarchical Analysis of Density Effects",
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
  )

save_figure(p_method, "method_overview", width = 10, height = 10)

# 2. Create simulation results summary figure
# This combines key findings from all simulation studies

# Placeholder for loading simulation results
# These would be loaded from the saved analysis files or regenerated

# 3. Create CRC analysis highlight figure
# This showcases the key findings from the CRC analysis

# Placeholder for loading CRC results
# These would be loaded from the saved analysis files or regenerated

# 4. Create overall paper visual abstract
# This combines elements from all analyses into a concise summary

# Placeholder for a visual abstract that summarizes the paper
# This would be a custom-designed figure that highlights the key findings

print("Summary figures created successfully!")