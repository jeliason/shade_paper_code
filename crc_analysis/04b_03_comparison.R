# ============================================================================
# COMPARISON: Combine Results and Create Summary Tables
# ============================================================================
# Module 3 of method comparison analysis
# Combines SHADE group differences with mFPCA/SOFR results

cat("\n=== Module 3: Comparison Analysis ===\n")

# ============================================================================
# COMBINE ALL DETECTIONS
# ============================================================================

all_detections <- bind_rows(
  # SHADE group differences
  shade_group_differences %>%
    mutate(
      method = "SHADE",
      level = "Group",
      Group = "Difference",
      p_value = NA_real_
    ) %>%
    select(method, target, source, level, Group, detected, p_value),

  # G-cross mFPCA/SOFR
  gcross_sofr_results %>%
    mutate(
      level = "Group",
      Group = "Overall"
    ) %>%
    select(method, target, source, level, Group, detected, p_value),

  # L-cross mFPCA/SOFR
  lcross_sofr_results %>%
    mutate(
      level = "Group",
      Group = "Overall"
    ) %>%
    select(method, target, source, level, Group, detected, p_value)
)

# Save results
saveRDS(all_detections, paste0(path, "method_comparison_detections_sofr.rds"))
cat("\n✓ Saved all detections to:", paste0(path, "method_comparison_detections_sofr.rds"), "\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n=== Summary Statistics ===\n")

summary_stats <- all_detections %>%
  group_by(method) %>%
  summarise(
    total_pairs = n(),
    n_detected = sum(detected, na.rm = TRUE),
    prop_detected = n_detected / total_pairs,
    mean_p_value = mean(p_value, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

# ============================================================================
# DETECTION AGREEMENT ANALYSIS
# ============================================================================

cat("\n=== Detection Agreement Across Methods ===\n")

# Create wide format for comparison
detections_wide <- all_detections %>%
  select(target, source, method, detected) %>%
  pivot_wider(
    names_from = method,
    values_from = detected,
    id_cols = c(target, source)
  ) %>%
  mutate(
    pair = paste(source, "->", target),
    n_methods_detected = rowSums(across(starts_with(c("SHADE", "G-cross", "L-cross"))), na.rm = TRUE)
  )

# Count agreement patterns
agreement_summary <- detections_wide %>%
  count(SHADE, `G-cross (mFPCA)`, `L-cross (mFPCA)`, name = "n_pairs") %>%
  arrange(desc(n_pairs))

print(agreement_summary)

# Which pairs detected by all methods?
all_methods_agree <- detections_wide %>%
  filter(n_methods_detected == 3)

cat("\nPairs detected by all methods:\n")
print(all_methods_agree %>% select(pair))

# Which pairs show disagreement?
disagreement <- detections_wide %>%
  filter(n_methods_detected > 0 & n_methods_detected < 3)

cat("\nPairs with partial detection:\n")
print(disagreement %>% select(pair, SHADE, `G-cross (mFPCA)`, `L-cross (mFPCA)`))

cat("✓ Module 3 complete\n")
