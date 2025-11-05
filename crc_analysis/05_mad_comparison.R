library(tidyverse)
library(patchwork)
library(SHADE)

# Load data
source("utils.R")
source("crc_analysis/utils.R")

# Assuming local_sic_df and pt_data are already loaded or created
# from 03_analyze_results.R
path <- "./crc_analysis/data/"
pt_data <- read_csv("crc_analysis/data/CRC_pt_metadata.csv")

figures_folder <- "./manuscript/images/CRC_analysis_paper/"
theme_set(theme_bw(base_size = 14, base_family = 'Helvetica') +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

fsave <- \(fname, height = 5, width = 8) {
  ggsave(paste0(figures_folder, fname, ".pdf"),
         device = cairo_pdf,
         height = height,
         width = width,
         units = "in")
}

# Load local SIC data (assuming this has been created in 03_analyze_results.R)
# If not, you'll need to generate it first
if (!exists("local_sic_df")) {
  stop("local_sic_df not found. Please run 03_analyze_results.R first to generate SIC data.")
}

# Merge patient metadata with SIC data
local_sic_with_group <- local_sic_df %>%
  left_join(pt_data %>% select(Patient, Group), by = "Patient")

# ============================================================================
# Patient-level MAD comparison between CLR and DII
# ============================================================================

# Step 1: Compute cohort means per group (CLR vs DII)
cohort_means_by_group <- local_sic_with_group %>%
  filter(level == "Patient") %>%
  group_by(x, Group, Target, Source) %>%
  summarise(mn_cohort = median(mn, na.rm = TRUE), .groups = "drop")

# Step 2: Merge with patient-level data
patient_merged_by_group <- local_sic_with_group %>%
  filter(level == "Patient") %>%
  left_join(cohort_means_by_group, by = c("x", "Group", "Target", "Source"))

# Step 3: Compute MAD of patient deviation from group cohort mean, per distance
var_within_cohort_x_by_group <- patient_merged_by_group %>%
  group_by(x, Group, Target, Source) %>%
  summarise(
    var_within_cohort = mad(mn - mn_cohort, na.rm = TRUE),
    .groups = "drop"
  )

# Step 4: Average across distances
var_within_cohort_by_group <- var_within_cohort_x_by_group %>%
  group_by(Group, Target, Source) %>%
  summarise(var_within_cohort = median(var_within_cohort, na.rm = TRUE),
            .groups = "drop")

# ============================================================================
# Image-level MAD comparison between CLR and DII
# ============================================================================

# Step 1: Compute patient means per group
patient_means_by_group <- local_sic_with_group %>%
  filter(level == "Patient") %>%
  group_by(x, Patient, Group, Target, Source) %>%
  summarise(mn_patient = median(mn, na.rm = TRUE), .groups = "drop")

# Step 2: Merge with image-level data
image_merged_by_group <- local_sic_with_group %>%
  filter(level == "Image") %>%
  left_join(patient_means_by_group, by = c("x", "Patient", "Group", "Target", "Source"))

# Step 3: Compute MAD of image deviation from patient mean, per distance
var_within_patient_x_by_group <- image_merged_by_group %>%
  group_by(x, Group, Target, Source) %>%
  summarise(
    var_within_patient = mad(mn - mn_patient, na.rm = TRUE),
    .groups = "drop"
  )

# Step 4: Average across distances
var_within_patient_by_group <- var_within_patient_x_by_group %>%
  group_by(Group, Target, Source) %>%
  summarise(var_within_patient = median(var_within_patient, na.rm = TRUE),
            .groups = "drop")

# ============================================================================
# Visualization: Heatmaps comparing CLR vs DII
# ============================================================================

# Patient-level MAD heatmaps (between-patient heterogeneity)
p_patient_CLR <- var_within_cohort_by_group %>%
  filter(Group == "CLR") %>%
  ggplot(aes(x = Target, y = Source, fill = var_within_cohort)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Patient-level MAD",
                       limits = range(var_within_cohort_by_group$var_within_cohort)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(title = "CLR: Between-patient heterogeneity",
       x = "Target cell type",
       y = "Source cell type")

p_patient_DII <- var_within_cohort_by_group %>%
  filter(Group == "DII") %>%
  ggplot(aes(x = Target, y = Source, fill = var_within_cohort)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Patient-level MAD",
                       limits = range(var_within_cohort_by_group$var_within_cohort)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top") +
  labs(title = "DII: Between-patient heterogeneity",
       x = "Target cell type")

# Image-level MAD heatmaps (within-patient heterogeneity)
p_image_CLR <- var_within_patient_by_group %>%
  filter(Group == "CLR") %>%
  ggplot(aes(x = Target, y = Source, fill = var_within_patient)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Image-level MAD",
                       limits = range(var_within_patient_by_group$var_within_patient)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(title = "CLR: Within-patient heterogeneity",
       x = "Target cell type",
       y = "Source cell type")

p_image_DII <- var_within_patient_by_group %>%
  filter(Group == "DII") %>%
  ggplot(aes(x = Target, y = Source, fill = var_within_patient)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Image-level MAD",
                       limits = range(var_within_patient_by_group$var_within_patient)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top") +
  labs(title = "DII: Within-patient heterogeneity",
       x = "Target cell type")

# Combine plots
(p_patient_CLR + p_patient_DII) / (p_image_CLR + p_image_DII) +
  plot_annotation(tag_levels = 'A')

fsave("mad_comparison_clr_dii", width = 12, height = 10)

# ============================================================================
# Compute differences: DII - CLR
# ============================================================================

# Patient-level MAD difference
mad_diff_patient <- var_within_cohort_by_group %>%
  pivot_wider(names_from = Group, values_from = var_within_cohort,
              names_prefix = "MAD_") %>%
  mutate(MAD_diff = MAD_DII - MAD_CLR)

# Image-level MAD difference
mad_diff_image <- var_within_patient_by_group %>%
  pivot_wider(names_from = Group, values_from = var_within_patient,
              names_prefix = "MAD_") %>%
  mutate(MAD_diff = MAD_DII - MAD_CLR)

# Visualize differences as diverging heatmaps
p_diff_patient <- mad_diff_patient %>%
  ggplot(aes(x = Target, y = Source, fill = MAD_diff)) +
  geom_tile() +
  scale_fill_gradient2(low = "purple", mid = "white", high = "green",
                       midpoint = 0, name = "DII - CLR\n(Patient MAD)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(title = "Patient-level heterogeneity difference",
       x = "Target cell type",
       y = "Source cell type")

p_diff_image <- mad_diff_image %>%
  ggplot(aes(x = Target, y = Source, fill = MAD_diff)) +
  geom_tile() +
  scale_fill_gradient2(low = "purple", mid = "white", high = "green",
                       midpoint = 0, name = "DII - CLR\n(Image MAD)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(title = "Image-level heterogeneity difference",
       x = "Target cell type",
       y = "Source cell type")

p_diff_patient + p_diff_image +
  plot_annotation(tag_levels = 'A')

fsave("mad_difference_clr_dii", width = 12, height = 5)

# Save summary statistics
write_csv(mad_diff_patient, paste0(path, "mad_patient_clr_dii_comparison.csv"))
write_csv(mad_diff_image, paste0(path, "mad_image_clr_dii_comparison.csv"))

# Print summary
cat("\n=== Patient-level MAD comparison ===\n")
cat("Interactions with largest DII > CLR heterogeneity:\n")
mad_diff_patient %>%
  arrange(desc(MAD_diff)) %>%
  head(5) %>%
  print()

cat("\nInteractions with largest CLR > DII heterogeneity:\n")
mad_diff_patient %>%
  arrange(MAD_diff) %>%
  head(5) %>%
  print()

cat("\n=== Image-level MAD comparison ===\n")
cat("Interactions with largest DII > CLR heterogeneity:\n")
mad_diff_image %>%
  arrange(desc(MAD_diff)) %>%
  head(5) %>%
  print()

cat("\nInteractions with largest CLR > DII heterogeneity:\n")
mad_diff_image %>%
  arrange(MAD_diff) %>%
  head(5) %>%
  print()
