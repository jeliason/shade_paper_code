# Identify Binary Tumor/Stromal Compartments
#
# Purpose: Delineate tumor vs stromal compartments using tumor cell density
# Approach: Kernel density estimation on tumor cells, median threshold split

library(tidyverse)
library(spatstat)
library(Matrix)

# ============================================================================
# 1. LOAD CRC DATA
# ============================================================================

cat("Loading CRC data...\n")
df_raw <- read_csv("crc_analysis/data/CRC_cleaned.csv")
pt_data <- read_csv("crc_analysis/data/CRC_pt_metadata.csv")

# Recode cell types to match analysis convention
df_raw <- df_raw %>%
  mutate(type = as.factor(type)) %>%
  mutate(type = fct_recode(
    type,
    "CAFs" = "smooth muscle",
    "hybrid E/M" = "stroma",
    "TAMs" = "CD163+ macros",
    "CTLs" = "CD8+ T cells"
  ))

spots <- pt_data %>% pull(Spot)

# ============================================================================
# 2. BINARY COMPARTMENT DETECTION FUNCTION
# ============================================================================

#' Create binary tumor/stromal compartment field
#'
#' @param image_data Data frame with columns X, Y, type for a single image
#' @param sigma Bandwidth for kernel density estimation (in microns)
#' @param threshold Absolute density threshold (cells/μm²). If NULL, uses 75th percentile
#' @return List containing:
#'   - density_im: tumor density spatial field
#'   - compartment_im: binary compartment field (0 = stromal, 1 = tumor)
#'   - threshold: density threshold used
create_binary_compartment_field <- function(image_data,
                                           sigma = 100,
                                           threshold = NULL) {

  # Get tumor cells
  tumor_cells <- image_data %>% filter(type == "tumor cells")

  if (nrow(tumor_cells) < 5) {
    cat("  Warning: <5 tumor cells, skipping\n")
    return(NULL)
  }

  # Create spatial window
  W <- owin(
    xrange = c(min(image_data$X), max(image_data$X)),
    yrange = c(min(image_data$Y), max(image_data$Y))
  )

  # Tumor point pattern
  tumor_ppp <- ppp(tumor_cells$X, tumor_cells$Y, window = W)

  # Estimate density - continuous spatial field
  density_im <- density(tumor_ppp, sigma = sigma)

  # Set threshold
  if (is.null(threshold)) {
    # Fallback: use 75th percentile (image-specific)
    threshold_used <- quantile(density_im, probs = 0.75, na.rm = TRUE)
  } else {
    # Use absolute threshold (same across all images)
    threshold_used <- threshold
  }

  compartment_im <- eval.im(as.integer(density_im > threshold_used))

  return(list(
    density_im = density_im,
    compartment_im = compartment_im,
    threshold = threshold_used
  ))
}

# ============================================================================
# 3. CREATE COMPARTMENT FIELDS FOR ALL IMAGES
# ============================================================================

cat("\nCreating binary tumor/stromal compartment fields...\n")

# SET THRESHOLD HERE:
# Option 1: Use NULL for adaptive (75th percentile per image)
# Option 2: Set absolute value after reviewing density distribution
# Example values to try: 0.001, 0.002, 0.005 cells/μm²
DENSITY_THRESHOLD <- 3.4e-5  # Change this after reviewing density_distribution.png

if (is.null(DENSITY_THRESHOLD)) {
  cat("Using adaptive threshold (75th percentile per image)\n")
} else {
  cat(sprintf("Using absolute threshold: %.6f cells/μm²\n", DENSITY_THRESHOLD))
}

compartment_fields <- lapply(spots, function(spot) {
  img_data <- df_raw %>% filter(Spot == spot)
  result <- create_binary_compartment_field(
    img_data,
    sigma = 100,
    threshold = DENSITY_THRESHOLD
  )
  if (!is.null(result)) {
    result$spot <- spot
  }
  result
})
names(compartment_fields) <- spots

# Remove NULL entries
is_valid <- !sapply(compartment_fields, is.null)
compartment_fields <- compartment_fields[is_valid]

cat(sprintf("Created compartment fields for %d/%d images\n",
            length(compartment_fields), length(spots)))

# ============================================================================
# 3B. EXPLORE DENSITY DISTRIBUTIONS TO INFORM THRESHOLD CHOICE
# ============================================================================

cat("\n--- Density Distribution Analysis ---\n")

# Extract all density values across all images
all_densities <- lapply(compartment_fields, function(field) {
  as.vector(field$density_im$v)
}) %>% unlist()

cat(sprintf("Total density measurements: %d\n", length(all_densities)))
cat("\nDensity quantiles (cells/μm²):\n")
print(quantile(all_densities, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1), na.rm = TRUE))

cat("\nSummary statistics:\n")
cat(sprintf("  Mean: %.6f cells/μm²\n", mean(all_densities, na.rm = TRUE)))
cat(sprintf("  Median: %.6f cells/μm²\n", median(all_densities, na.rm = TRUE)))
cat(sprintf("  SD: %.6f cells/μm²\n", sd(all_densities, na.rm = TRUE)))

# Plot histogram to visualize distribution
p_density_hist <- ggplot(data.frame(density = log(all_densities)), aes(x = density)) +
  geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = median(all_densities, na.rm = TRUE),
             color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = quantile(all_densities, 0.75, na.rm = TRUE),
             color = "orange", linetype = "dashed", linewidth = 1) +
  # scale_x_continuous(trans = "log10") +
  labs(
    title = "Distribution of tumor cell density across all images",
    subtitle = "Red = median, Orange = 75th percentile",
    x = "Tumor cell density (cells/μm², log scale)",
    y = "Count"
  ) +
  theme_minimal()

p_density_hist

ggsave(
  "./compartment_analysis/compartment_analysis_figures/density_distribution.png",
  p_density_hist,
  width = 8,
  height = 5,
  dpi = 300
)

cat("\nSaved density distribution plot\n")
cat("RECOMMENDATION: Review density_distribution.png to choose an appropriate threshold\n")
cat("------------------------------------\n\n")

# ============================================================================
# 4. EVALUATE COMPARTMENT AT CELL LOCATIONS
# ============================================================================

cat("\nEvaluating compartments at cell locations...\n")

df_with_compartments <- df_raw %>%
  filter(Spot %in% names(compartment_fields)) %>%
  group_by(Spot) %>%
  group_modify(function(img_data, keys) {
    spot <- as.character(keys$Spot[1])
    comp_field <- compartment_fields[[spot]]$compartment_im

    # Look up compartment at each cell location
    cell_ppp <- ppp(img_data$X, img_data$Y, window = Window(comp_field))
    img_data$compartment <- safelookup(comp_field, cell_ppp)
    img_data
  }) %>%
  ungroup()

# ============================================================================
# 5. SUMMARY STATISTICS
# ============================================================================

cat("\nCompartment summary:\n")
compartment_summary <- df_with_compartments %>%
  count(Spot, compartment) %>%
  group_by(compartment) %>%
  summarize(
    n_images = n(),
    mean_cells = mean(n),
    sd_cells = sd(n)
  )

print(compartment_summary)

# Cell type composition by compartment
cat("\nCell type composition:\n")
cell_type_by_compartment <- df_with_compartments %>%
  count(compartment, type) %>%
  group_by(compartment) %>%
  mutate(proportion = n / sum(n)) %>%
  arrange(compartment, desc(proportion))

print(cell_type_by_compartment, n = 20)

# ============================================================================
# 6. VISUALIZATIONS
# ============================================================================

cat("\nCreating visualizations...\n")

# Detailed view with density
detailed_spot <- example_spots[1]
detailed_data <- df_with_compartments %>% filter(Spot == detailed_spot)
detailed_density <- as.data.frame(compartment_fields[[detailed_spot]]$density_im)
detailed_compartment <- as.data.frame(compartment_fields[[detailed_spot]]$compartment_im)

p_detailed <- ggplot() +
  geom_tile(
    data = detailed_density,
    aes(x = x, y = y, fill = value)
  ) +
  scale_fill_viridis_c(
    name = "Tumor\ndensity",
    option = "magma"
  ) +
  geom_contour(
    data = detailed_compartment,
    aes(x = x, y = y, z = value),
    breaks = 0.5,
    color = "white",
    linewidth = 1.5,
    linetype = "dashed"
  ) +
  geom_point(
    data = detailed_data,
    aes(x = X, y = Y, color = type),
    size = 1.5, alpha = 0.7
  ) +
  scale_color_manual(
    values = as.vector(pals::glasbey()),
    name = "Cell type"
  ) +
  labs(
    title = paste("Detailed view:", detailed_spot),
    subtitle = "Tumor density field with compartment boundary (dashed line)",
    x = "X (μm)",
    y = "Y (μm)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  coord_fixed()

p_detailed

# ============================================================================
# 7. SAVE DATA
# ============================================================================

cat("\nSaving compartment data...\n")

# Save compartment fields
saveRDS(
  compartment_fields,
  "./compartment_analysis/data/binary_compartment_fields.rds"
)

# Save cell data with compartments
saveRDS(
  df_with_compartments,
  "./compartment_analysis/data/CRC_with_binary_compartments.rds"
)
