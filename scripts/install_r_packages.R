#!/usr/bin/env Rscript
# Install R packages required for SHADE paper code
# Run this on HPC before submitting jobs to ensure all dependencies are available

cat("=== Installing R packages for SHADE paper code ===\n\n")

cat(sprintf("Using R library path: %s\n", .libPaths()[1]))
cat(sprintf("R version: %s\n\n", R.version.string))

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# List of CRAN packages
cran_packages <- c(
  "tidyverse",      # Data manipulation and plotting
  "Matrix",         # Sparse matrices
  "spatstat",       # Spatial statistics
  "cmdstanr",       # Stan interface
  "posterior",      # Posterior analysis
  "ggdist",         # Distribution visualizations
  "survival",       # Survival analysis
  "latex2exp",      # LaTeX expressions in plots
  "patchwork",      # Combine plots
  "splines2",       # Spline functions
  "survminer",      # Survival plot visualization
  "kableExtra",     # Table formatting
  "tictoc"          # Timing utilities
)

cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  - Installing %s\n", pkg))
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

cat("\n")

# cmdstanr requires additional setup
cat("Checking cmdstanr installation...\n")
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  cat("  - Installing cmdstanr from GitHub\n")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("stan-dev/cmdstanr")
}

# Check if CmdStan is installed
library(cmdstanr)
if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  cat("  - CmdStan not found, installing...\n")
  cat("    This may take several minutes...\n")
  cmdstanr::install_cmdstan()
} else {
  cat(sprintf("  ✓ CmdStan %s already installed\n", cmdstanr::cmdstan_version()))
}

cat("\n")

# SHADE package installation
cat("Installing SHADE package...\n")

# Check environment variable for SHADE source
shade_source <- Sys.getenv("SHADE_SOURCE", "")
force_reinstall <- Sys.getenv("SHADE_FORCE_REINSTALL", "false") == "true"

if (shade_source == "") {
  cat("  ⚠️  SHADE_SOURCE environment variable not set\n")
  cat("     Please set one of:\n")
  cat("     - SHADE_SOURCE=github:username/repo  (for GitHub)\n")
  cat("     - SHADE_SOURCE=/path/to/local/SHADE (for local install)\n")
  cat("\n")
  cat("     Example:\n")
  cat("     export SHADE_SOURCE=github:yourusername/SHADE\n")
  cat("     Rscript scripts/install_r_packages.R\n")
  cat("\n")
  cat("     Optional: Set SHADE_FORCE_REINSTALL=true to force reinstall\n")
  cat("\n")
} else if (grepl("^github:", shade_source)) {
  # Install from GitHub
  repo <- sub("^github:", "", shade_source)
  if (force_reinstall) {
    cat(sprintf("  - Force reinstalling SHADE from GitHub: %s\n", repo))
  } else {
    cat(sprintf("  - Installing SHADE from GitHub: %s\n", repo))
  }
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github(repo, force = force_reinstall)
} else if (file.exists(shade_source)) {
  # Install from local path
  if (force_reinstall) {
    cat(sprintf("  - Force reinstalling SHADE from local path: %s\n", shade_source))
    # Remove existing SHADE first
    if (requireNamespace("SHADE", quietly = TRUE)) {
      remove.packages("SHADE")
    }
  } else {
    cat(sprintf("  - Installing SHADE from local path: %s\n", shade_source))
  }
  install.packages(shade_source, repos = NULL, type = "source")
} else {
  cat(sprintf("  ⚠️  Invalid SHADE_SOURCE: %s\n", shade_source))
  cat("     Path does not exist or invalid format\n")
}

# Verify SHADE is available
if (requireNamespace("SHADE", quietly = TRUE)) {
  cat("  ✓ SHADE package installed successfully\n")
} else {
  cat("  ⚠️  SHADE package not available - you'll need to install it manually\n")
}

cat("\n=== Installation complete ===\n")
cat("\nVerifying all packages are loadable...\n")

# Test loading key packages
test_packages <- c("tidyverse", "spatstat", "Matrix", "cmdstanr", "posterior", "SHADE")
all_ok <- TRUE

for (pkg in test_packages) {
  tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    cat(sprintf("  ✓ %s loads successfully\n", pkg))
  }, error = function(e) {
    cat(sprintf("  ✗ %s failed to load: %s\n", pkg, e$message))
    all_ok <<- FALSE
  })
}

cat("\n")
if (all_ok) {
  cat("✓ All packages installed and verified!\n")
} else {
  cat("⚠️  Some packages failed - please check the errors above\n")
  quit(status = 1)
}
