# ============================================================================
# COMPARISON: Combine Results and Create Summary Tables
# ============================================================================
# Module 3 of method comparison analysis
# Combines SHADE group differences with mFPCA results

cat("\n=== Module 3: Comparison Analysis ===\n")

# ============================================================================
# COMBINE ALL DETECTIONS
# ============================================================================

all_detections <- shade_group_differences %>%
  mutate(
    method = "SHADE",
    level = "Group",
    Group = "Difference"
  ) %>%
  select(method, target, source, level, Group, detected)

# Save results
saveRDS(all_detections, paste0(path, "method_comparison_detections.rds"))
cat("\n✓ Saved SHADE group differences to:", paste0(path, "method_comparison_detections.rds"), "\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n=== SHADE Group Difference Summary ===\n")

summary_stats <- all_detections %>%
  summarise(
    total_pairs = n(),
    n_detected = sum(detected, na.rm = TRUE),
    prop_detected = n_detected / total_pairs
  )

print(summary_stats)

cat("✓ Module 3 complete\n")
