#!/usr/bin/env Rscript
# Reinstall SHADE package from local repository

library(devtools)

# Path to local SHADE repository
shade_path <- "/Users/joeleliason/Projects/SHADE"

cat("Reinstalling SHADE package from:", shade_path, "\n")

# Remove existing installation
try(remove.packages("SHADE"), silent = TRUE)

# Install from local source
install(shade_path, upgrade = "never")

cat("\nSHADE package reinstalled successfully!\n")

# Verify installation
cat("\nVerifying installation...\n")
library(SHADE)
cat("SHADE version:", as.character(packageVersion("SHADE")), "\n")
