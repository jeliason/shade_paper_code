#!/usr/bin/env Rscript
# Setup script for sim_shade_comparison
# Compiles Stan model on HPC before running simulations

library(cmdstanr)

cat("=== Compiling super_simple_shade.stan ===\n")

# Path to Stan model
stan_file <- "sim_shade_comparison/super_simple_shade.stan"

if (!file.exists(stan_file)) {
  stop("Stan file not found: ", stan_file)
}

# Force compilation
cat("Compiling Stan model with HPC toolchain...\n")
mod <- cmdstan_model(stan_file, compile=TRUE, force_recompile = TRUE)

cat("✓ Stan model compiled successfully\n")
cat("Executable: ", mod$exe_file(), "\n")

# Test that the model can be loaded
cat("\nTesting model can be loaded...\n")
test_mod <- cmdstan_model(stan_file)
cat("✓ Model loads successfully\n")

cat("\n=== Setup complete ===\n")
