library(tidyverse)
library(patchwork)

# set paths
data_path <- "./sim_timing/data/"
fig_path <- "./manuscript/images/sim_timing/"

# Ensure figure directory exists
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

# Load standardized analysis summary
analysis_summary <- readRDS(paste0(data_path, "analysis_summary.rds"))

timing_summary <- analysis_summary$timing_summary
scaling_feature <- analysis_summary$scaling_exponents$feature
scaling_total <- analysis_summary$scaling_exponents$total

# Create reference lines for O(n) and O(n^2) scaling
n_range <- seq(min(timing_summary$total_cells), max(timing_summary$total_cells), length.out = 100)

# Normalize reference lines to pass through first data point for visual comparison
n0 <- timing_summary$total_cells[1]
t_feature_0 <- timing_summary$mean_feature_sec[1]
t_total_0 <- timing_summary$mean_total_sec[1]

ref_linear_feature <- t_feature_0 * (n_range / n0)
ref_quadratic_feature <- t_feature_0 * (n_range / n0)^2

ref_linear_total <- t_total_0 * (n_range / n0)
ref_quadratic_total <- t_total_0 * (n_range / n0)^2

# Prepare reference data for plotting
ref_data_feature <- tibble(
  n = rep(n_range, 2),
  time = c(ref_linear_feature, ref_quadratic_feature),
  scaling = rep(c("O(n)", "O(n²)"), each = length(n_range))
)

ref_data_total <- tibble(
  n = rep(n_range, 2),
  time = c(ref_linear_total, ref_quadratic_total),
  scaling = rep(c("O(n)", "O(n²)"), each = length(n_range))
)

# Plot 1: Feature construction time
p1 <- ggplot() +
  geom_line(data = ref_data_feature, aes(x = n, y = time, linetype = scaling, color = scaling),
            alpha = 0.4, linewidth = 0.6) +
  geom_pointrange(data = timing_summary,
                  aes(x = total_cells, y = mean_feature_sec,
                      ymin = mean_feature_sec - sd_feature_sec,
                      ymax = mean_feature_sec + sd_feature_sec),
                  color = "#2166AC", size = 0.6, linewidth = 0.8) +
  scale_x_log10(
    breaks = c(1000, 2500, 5000, 10000, 25000, 50000, 100000),
    labels = c("1K", "2.5K", "5K", "10K", "25K", "50K", "100K")
  ) +
  scale_y_log10() +
  scale_color_manual(values = c("O(n)" = "gray50", "O(n²)" = "gray30")) +
  scale_linetype_manual(values = c("O(n)" = "dotted", "O(n²)" = "dashed")) +
  labs(
    x = "Number of cells",
    y = "Time (seconds)",
    title = "Feature Construction Time",
    subtitle = bquote("Empirical scaling:" ~ O(n^.(round(scaling_feature, 2))))
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.key.size = unit(0.8, "lines"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9)
  )

# Plot 2: Total model fitting time
p2 <- ggplot() +
  geom_line(data = ref_data_total, aes(x = n, y = time, linetype = scaling, color = scaling),
            alpha = 0.4, linewidth = 0.6) +
  geom_pointrange(data = timing_summary,
                  aes(x = total_cells, y = mean_total_sec,
                      ymin = mean_total_sec - sd_total_sec,
                      ymax = mean_total_sec + sd_total_sec),
                  color = "#B2182B", size = 0.6, linewidth = 0.8) +
  scale_x_log10(
    breaks = c(1000, 2500, 5000, 10000, 25000, 50000, 100000),
    labels = c("1K", "2.5K", "5K", "10K", "25K", "50K", "100K")
  ) +
  scale_y_log10() +
  scale_color_manual(values = c("O(n)" = "gray50", "O(n²)" = "gray30")) +
  scale_linetype_manual(values = c("O(n)" = "dotted", "O(n²)" = "dashed")) +
  labs(
    x = "Number of cells",
    y = "Time (seconds)",
    title = "Total Model Fitting Time",
    subtitle = bquote("Empirical scaling:" ~ O(n^.(round(scaling_total, 2))))
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.key.size = unit(0.8, "lines"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9)
  )

# Combine plots
combined <- p1 + p2 +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(plot.tag = element_text(face = 'bold', size = 12))
  )
combined
# Save figure
ggsave(
  paste0(fig_path, "scaling_results.pdf"),
  combined,
  width = 10,
  height = 4.5,
  units = "in"
)

print(paste("Figure saved to", paste0(fig_path, "scaling_results.pdf")))

# Also save individual panels for flexibility
ggsave(paste0(fig_path, "feature_scaling.pdf"), p1, width = 5, height = 4.5)
ggsave(paste0(fig_path, "total_scaling.pdf"), p2, width = 5, height = 4.5)

print("Individual panels saved.")
