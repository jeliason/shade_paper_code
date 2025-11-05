library(tidyverse)

# Load utilities
source("utils.R")

# Set environment and get data path
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  sim_indices <- 1:10  # subset for local testing
} else {
  sim_indices <- NULL  # Process all available results
}

# Get data path (handles HPC vs local automatically)
path <- get_data_path("sim_timing")
print(paste("Data path:", path))

# parameters (must match 01_generate_data.R)
grid <- expand.grid(
  total_cells = c(5000, 10000, 25000, 50000, 100000, 250000),
  sim = 1:20
)

# grid <- grid[sim_indices, ]

# Load all timing results
timing_data <- pmap(grid, \(total_cells, sim) {
  file_timing <- paste0(path, "timing_sim", sim, "_cells_", total_cells, ".rds")

  if(file.exists(file_timing)) {
    timing <- readRDS(file_timing)
    return(tibble(
      total_cells = timing$total_cells,
      sim = timing$sim,
      time_feature_sec = timing$time_feature_sec,
      time_fitting_sec = timing$time_fitting_sec,
      time_total_sec = timing$time_total_sec
    ))
  } else {
    return(NULL)
  }
}) %>%
  bind_rows()

print(paste("Loaded", nrow(timing_data), "timing results"))

# Compute summary statistics by cell count
timing_summary <- timing_data %>%
  group_by(total_cells) %>%
  summarise(
    n_reps = n(),
    mean_feature_sec = mean(time_feature_sec),
    sd_feature_sec = sd(time_feature_sec),
    mean_fitting_sec = mean(time_fitting_sec),
    sd_fitting_sec = sd(time_fitting_sec),
    mean_total_sec = mean(time_total_sec),
    sd_total_sec = sd(time_total_sec),
    .groups = "drop"
  )

print("Timing summary by cell count:")
print(timing_summary)

# Fit log-log regression to estimate scaling exponent
# For feature construction (expected ~quadratic: O(n^2))
fit_feature <- lm(log(mean_feature_sec) ~ log(total_cells), data = timing_summary)
scaling_feature <- coef(fit_feature)[2]
print(paste("Feature construction scaling exponent:", round(scaling_feature, 3)))

# For total fitting time (expected sublinear due to parallelization)
fit_total <- lm(log(mean_total_sec) ~ log(total_cells), data = timing_summary)
scaling_total <- coef(fit_total)[2]
print(paste("Total fitting scaling exponent:", round(scaling_total, 3)))

# Save standardized analysis summary for fetching/plotting
analysis_summary <- list(
  timing_data = timing_data,
  timing_summary = timing_summary,
  scaling_exponents = list(
    feature = scaling_feature,
    total = scaling_total
  ),
  fit_feature = fit_feature,
  fit_total = fit_total,
  metadata = list(
    sim_dir = "sim_timing",
    date_analyzed = Sys.time(),
    n_simulations = nrow(timing_data),
    cell_counts = unique(timing_data$total_cells)
  )
)

saveRDS(analysis_summary, paste0(path, "analysis_summary.rds"))
print(paste("Saved analysis summary to:", paste0(path, "analysis_summary.rds")))

# Print summary table
cat("\n=== Computational Scaling Summary ===\n")
cat(sprintf("Cell Count | Feature (sec) | Fitting (sec) | Total (sec)\n"))
cat(sprintf("---------- | ------------- | ------------- | -----------\n"))
for(i in 1:nrow(timing_summary)) {
  cat(sprintf("%10d | %6.2f ± %5.2f | %6.2f ± %5.2f | %6.2f ± %5.2f\n",
              timing_summary$total_cells[i],
              timing_summary$mean_feature_sec[i],
              timing_summary$sd_feature_sec[i],
              timing_summary$mean_fitting_sec[i],
              timing_summary$sd_fitting_sec[i],
              timing_summary$mean_total_sec[i],
              timing_summary$sd_total_sec[i]))
}
cat("\n")
cat(sprintf("Feature construction scales as O(n^%.2f)\n", scaling_feature))
cat(sprintf("Total fitting time scales as O(n^%.2f)\n", scaling_total))
