# ============================================================================
# METHOD COMPARISON: SHADE vs mFPCA for CRC Data
# ============================================================================
# Main script for visual comparison of SHADE group differences with
# multilevel functional PCA (mFPCA) group-level curves for G-cross
# and L-cross summary statistics
#
# This script sources modular analysis files for better organization:
#   - 04b_00_setup.R: Load data and prepare spatial patterns
#   - 04b_01_shade_analysis.R: SHADE group-level detections and differences
#   - 04b_02_mfpca.R: mFPCA for G-cross and L-cross (by group)
#   - 04b_03_comparison.R: Combine results and create summary tables
#   - 04b_04_visualizations.R: Create SHADE SIC figures
# ============================================================================

cat(paste(rep("=", 80), collapse=""), "\n")
cat("METHOD COMPARISON: SHADE vs mFPCA\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

start_time <- Sys.time()

# ============================================================================
# MODULE 0: SETUP
# ============================================================================

cat("\n*** Running Module 0: Setup ***\n")
source("crc_analysis/04b_00_setup.R")

# ============================================================================
# MODULE 1: SHADE ANALYSIS
# ============================================================================

cat("\n*** Running Module 1: SHADE Analysis ***\n")
source("crc_analysis/04b_01_shade_analysis.R")

# ============================================================================
# MODULE 2: mFPCA ANALYSIS
# ============================================================================

cat("\n*** Running Module 2: mFPCA Analysis ***\n")
source("crc_analysis/04b_02_mfpca.R")

# ============================================================================
# MODULE 3: COMPARISON
# ============================================================================

cat("\n*** Running Module 3: Comparison ***\n")
source("crc_analysis/04b_03_comparison.R")

# ============================================================================
# MODULE 4: VISUALIZATIONS
# ============================================================================

cat("\n*** Running Module 4: Visualizations ***\n")
source("crc_analysis/04b_04_visualizations.R")

# ============================================================================
# COMPLETE
# ============================================================================

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("✓ METHOD COMPARISON COMPLETE\n")
cat("Total runtime:", round(runtime, 2), "minutes\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")
