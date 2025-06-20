library(ggplot2)
library(dplyr)
library(tibble)
library(cowplot)
set.seed(42)

# Rare, tightly clustered A
cluster_center <- c(0.5, 0.5)
df_A_clustered <- data.frame(
  x = rnorm(30, mean = cluster_center[1], sd = 0.05),
  y = rnorm(30, mean = cluster_center[2], sd = 0.05),
  type = "A"
)

# Diffuse B
df_B_diffuse <- data.frame(
  x = runif(100, 0.1, 0.9),
  y = runif(100, 0.1, 0.9),
  type = "B"
)

df_structured <- bind_rows(df_A_clustered, df_B_diffuse)

plot_structured_points <- ggplot(df_structured, aes(x, y, color = type)) +
  geom_point(size = 1.5, alpha = 0.7) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "A is clustered, B is diffuse") +
  scale_color_manual(values = c("A" = "firebrick", "B" = "steelblue"))

distances <- seq(0, 0.2, length.out = 100)
sic_AtoB <- 1 + 4 * exp(-distances / 0.03)
sic_BtoA <- 1 + 1.2 * exp(-distances / 0.07)

df_sic_structured <- tibble(
  r = rep(distances, 2),
  value = c(sic_AtoB, sic_BtoA),
  direction = rep(c("A → B", "B → A"), each = length(distances))
)

plot_structured_sic <- ggplot(df_sic_structured, aes(x = r * 100, y = value, color = direction)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  labs(title = "SICs: A more predictive due to clustering", x = "Distance (µm)", y = "SIC(r)") +
  scale_color_manual(values = c("A → B" = "firebrick", "B → A" = "steelblue"))

panel_structured <- plot_grid(plot_structured_points, plot_structured_sic, ncol = 2, align = "h", rel_widths = c(1, 1.2))
print(panel_structured)
