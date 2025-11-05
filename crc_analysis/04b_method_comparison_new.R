# ============================================================================
# METHOD COMPARISON: SHADE vs mFPCA/SOFR for CRC Data
# ============================================================================
# Main script for comparing SHADE group differences with functional data
# analysis approaches (mFPCA + scalar-on-function regression) on G-cross
# and L-cross summary statistics
#
# This script sources modular analysis files for better organization:
#   - 04b_00_setup.R: Load data and prepare spatial patterns
#   - 04b_01_shade_analysis.R: SHADE group-level detections and differences
#   - 04b_02_mfpca_sofr.R: mFPCA/SOFR for G-cross and L-cross
#   - 04b_03_comparison.R: Combine results and create summary tables
#   - 04b_04_visualizations.R: Create all figures
# ============================================================================

cat("=" %>% rep(80) %>% paste(collapse=""), "\n")
cat("METHOD COMPARISON: SHADE vs mFPCA/SOFR\n")
cat("=" %>% rep(80) %>% paste(collapse=""), "\n\n")

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
# MODULE 2: mFPCA/SOFR ANALYSIS
# ============================================================================

cat("\n*** Running Module 2: mFPCA/SOFR Analysis ***\n")
source("crc_analysis/04b_02_mfpca_sofr.R")

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

cat("\n", "=" %>% rep(80) %>% paste(collapse=""), "\n")
cat("✓ METHOD COMPARISON COMPLETE\n")
cat("Total runtime:", round(runtime, 2), "minutes\n")
cat("=" %>% rep(80) %>% paste(collapse=""), "\n\n")
