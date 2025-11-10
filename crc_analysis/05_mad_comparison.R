library(tidyverse)
library(patchwork)
library(SHADE)

# Load data
source("utils.R")
source("crc_analysis/utils.R")

# Assuming local_sic_df and pt_data are already loaded or created
# from 03_analyze_results.R
path <- "./crc_analysis/data/"
pt_data <- read_csv("crc_analysis/data/CRC_pt_metadata.csv") %>%
  mutate(Patient = as.factor(Patient))

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
local_sic_df <- readRDS(paste0(path,"/local_sic_df.rds")) %>%
  mutate(Patient = as.factor(Patient))

# Merge patient metadata with SIC data
local_sic_with_group <- local_sic_df %>%
  left_join(pt_data %>% select(Patient, Group) %>% distinct(), by = "Patient")

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
  filter(level == "Spot") %>%
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

# Combine data for faceting
mad_combined <- bind_rows(
  var_within_cohort_by_group %>%
    mutate(Level = "Between-patient") %>%
    rename(MAD = var_within_cohort),
  var_within_patient_by_group %>%
    mutate(Level = "Within-patient") %>%
    rename(MAD = var_within_patient)
) %>%
  mutate(
    Group = factor(Group, levels = c("CLR", "DII")),
    Level = factor(Level, levels = c("Between-patient",
                                      "Within-patient"))
  )

# Create faceted plot
ggplot(mad_combined, aes(x = Target, y = Source, fill = MAD)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "MAD") +
  facet_grid(Level ~ Group) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        strip.background = element_rect(fill = "white", color = "black")) +
  labs(x = "Target cell type",
       y = "Source cell type")

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
# Combine difference data for faceting
mad_diff_combined <- bind_rows(
  mad_diff_patient %>%
    mutate(Level = "Between-patient"),
  mad_diff_image %>%
    mutate(Level = "Within-patient")
) %>%
  mutate(Level = factor(Level, levels = c("Between-patient", "Within-patient")))

ggplot(mad_diff_combined, aes(x = Target, y = Source, fill = MAD_diff)) +
  geom_tile() +
  scale_fill_gradient2(low = "purple", mid = "white", high = "green",
                       midpoint = 0, name = "DII - CLR") +
  facet_wrap(~ Level, nrow = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        strip.background = element_rect(fill = "white", color = "black")
      ) +
  labs(x = "Target cell type",
       y = "Source cell type")

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
