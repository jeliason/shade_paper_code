library(tidyverse)
library(spatstat)
library(patchwork)

#' Plot diagnostic visualizations for confounder simulations
#'
#' Creates 4-panel plot showing:
#' 1. Source contribution to target density
#' 2. Confounder contribution to target density
#' 3. Combined log-intensity (source + confounder)
#' 4. Actual simulated target points
#'
#' @param pattern_data List containing pattern info (from generate_spatial_patterns_*)
#' @param pattern_index Which pattern to plot (default: 1)
#' @param confounder_type Type of confounder: "gradient", "compartment", "unmeasured"
plot_confounder_diagnostics <- function(pattern_data, pattern_index = 1,
                                        confounder_type = "gradient") {

  pattern_info <- pattern_data[[pattern_index]]
  combined_pattern <- pattern_info$pattern

  W <- Window(combined_pattern)
  size_im <- max(W$xrange)

  # Extract cell types
  all_points <- as.data.frame(combined_pattern)
  t_cells <- all_points %>% filter(marks == 1)
  tumor_cells <- all_points %>% filter(marks == 3)

  # Recreate components for visualization
  potentials <- SHADE::make_rbfs(n_basis_functions = 3, max_dist = 75,
                                 basis_function_sigma = 15)
  image_coeffs <- pattern_info$image_coeffs

  # --- Panel 1: Source contribution ---
  if(nrow(t_cells) > 0) {
    t_ppp <- ppp(t_cells$x, t_cells$y, window = W)

    custom_kernel <- Vectorize(function(x, y) {
      d <- sqrt(x^2 + y^2)
      kernel_value <- sum(sapply(seq_along(image_coeffs), function(k) {
        image_coeffs[k] * potentials[[k]](d)
      }))
      return(kernel_value)
    })

    source_density <- smooth_density_fft(t_ppp, custom_kernel, resolution = 128)
  } else {
    source_density <- as.im(0, W = W)
  }

  # --- Panel 2: Confounder contribution ---
  if(confounder_type == "gradient") {
    gradient_type <- pattern_info$gradient_type
    gradient_strength <- pattern_info$gradient_strength

    gradient_fn <- function(x, y) {
      if(gradient_type == "linear") {
        return(gradient_strength * (x - size_im/2) / (size_im/2))
      } else if(gradient_type == "radial") {
        cx <- size_im / 2
        cy <- size_im / 2
        dist_from_center <- sqrt((x - cx)^2 + (y - cy)^2)
        max_dist <- sqrt(2) * size_im / 2
        return(gradient_strength * (dist_from_center / max_dist - 0.5))
      } else if(gradient_type == "sinusoidal") {
        return(gradient_strength * sin(2 * pi * x / size_im))
      }
    }

    confounder_density <- as.im(function(x, y) gradient_fn(x, y), W = W,
                                 eps = source_density$xstep)
    confounder_label <- paste0(gradient_type, " gradient (strength=",
                               gradient_strength, ")")

  } else if(confounder_type == "compartment") {
    # Recreate compartment structure from stored info
    compartment_image <- pattern_info$compartment_image
    compartment_effects <- pattern_info$compartment_effects

    # Create density image from compartment effects using stored compartment image
    effect_mat <- matrix(log(compartment_effects[compartment_image$v]),
                        nrow = nrow(compartment_image$v),
                        ncol = ncol(compartment_image$v))
    confounder_density <- im(effect_mat,
                             xcol = compartment_image$xcol,
                             yrow = compartment_image$yrow)

    confounder_label <- paste0(pattern_info$n_compartments, " compartments (effect=",
                               pattern_info$compartment_effect, ")")

  } else if(confounder_type == "unmeasured") {
    # Visualize unmeasured source
    unmeasured_pattern <- pattern_info$unmeasured_pattern
    if(npoints(unmeasured_pattern) > 0) {
      unmeasured_strength <- pattern_info$unmeasured_strength

      custom_kernel_unmeasured <- Vectorize(function(x, y) {
        d <- sqrt(x^2 + y^2)
        kernel_value <- unmeasured_strength * sum(sapply(seq_along(image_coeffs), function(k) {
          image_coeffs[k] * potentials[[k]](d)
        }))
        return(kernel_value)
      })

      confounder_density <- smooth_density_fft(unmeasured_pattern,
                                               custom_kernel_unmeasured,
                                               resolution = 128)
      confounder_label <- paste0("Unmeasured source (strength=",
                                 unmeasured_strength, "x)")
    } else {
      confounder_density <- as.im(0, W = W)
      confounder_label <- "Unmeasured source (none)"
    }
  }

  # --- Panel 3: Combined log-intensity ---
  combined_density <- source_density + confounder_density

  # Convert to data frames for ggplot
  source_df <- as.data.frame(source_density)
  confounder_df <- as.data.frame(confounder_density)
  combined_df <- as.data.frame(combined_density)

  # Create plots
  p1 <- ggplot(source_df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    geom_point(data = t_cells, aes(x = x, y = y),
               color = "white", size = 0.5, alpha = 0.6, inherit.aes = FALSE) +
    scale_fill_viridis_c(option = "magma", name = "Log\nintensity") +
    coord_equal() +
    labs(title = "Source contribution (T cells)",
         subtitle = paste(nrow(t_cells), "T cells")) +
    theme_minimal() +
    theme(legend.position = "right")

  p2 <- ggplot(confounder_df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "viridis", name = "Log\nintensity") +
    coord_equal() +
    labs(title = "Confounder contribution",
         subtitle = confounder_label) +
    theme_minimal() +
    theme(legend.position = "right")

  p3 <- ggplot(combined_df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma", name = "Log\nintensity") +
    coord_equal() +
    labs(title = "Combined log-intensity",
         subtitle = "Source + Confounder") +
    theme_minimal() +
    theme(legend.position = "right")

  p4 <- ggplot() +
    geom_point(data = t_cells, aes(x = x, y = y),
               color = "blue", size = 1, alpha = 0.5) +
    geom_point(data = tumor_cells, aes(x = x, y = y),
               color = "red", size = 1, alpha = 0.5) +
    coord_equal() +
    xlim(0, size_im) + ylim(0, size_im) +
    labs(title = "Simulated pattern",
         subtitle = paste(nrow(tumor_cells), "tumor cells (red),",
                         nrow(t_cells), "T cells (blue)")) +
    theme_minimal() +
    theme(panel.grid = element_blank())

  # Combine panels
  combined_plot <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title = paste("Confounder Diagnostic: Pattern", pattern_index),
      subtitle = paste("Confounder type:", confounder_type),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )

  return(combined_plot)
}

# Example usage function
#' Run a quick diagnostic plot for gradient confounder
#'
#' @param t_density Source cell density ("high" or "low")
#' @param tumor_density Target cell density ("high" or "low")
#' @param gradient_type Type of gradient ("linear", "radial", "sinusoidal")
#' @param gradient_strength Strength of gradient
example_gradient_diagnostic <- function(t_density = "high",
                                        tumor_density = "high",
                                        gradient_type = "radial",
                                        gradient_strength = 1.0) {

  source("./utils.R")
  source("sim_shade_comparison/comparison_functions.R")

  # Quick setup
  structure <- create_patient_structure(num_patients = 1, num_images = 1, seed = 2024)

  coefficients <- generate_spatial_coefficients(
    structure = structure,
    num_potentials = 3,
    cohort_mean = c(1.5, 1.0, 0.5),
    seed = 2024
  )

  density_params <- list(
    np_t = ifelse(t_density == "high", 150, 15),
    np_tumor = ifelse(tumor_density == "high", 150, 15),
    np_b = ifelse(t_density == "high", 150, 15)
  )

  # Generate pattern
  patterns <- generate_spatial_patterns_gradient(
    structure = structure,
    coefficients = coefficients,
    density_params = density_params,
    gradient_type = gradient_type,
    gradient_strength = gradient_strength,
    seed = 2024
  )

  # Plot
  p <- plot_confounder_diagnostics(patterns, pattern_index = 1,
                                   confounder_type = "gradient")
  print(p)

  return(invisible(patterns))
}

# Example for unmeasured source
example_unmeasured_diagnostic <- function(t_density = "high",
                                          tumor_density = "high",
                                          unmeasured_strength = 1.0) {

  source("./utils.R")
  source("sim_shade_comparison/comparison_functions.R")

  structure <- create_patient_structure(num_patients = 1, num_images = 1, seed = 2024)

  coefficients <- generate_spatial_coefficients(
    structure = structure,
    num_potentials = 3,
    cohort_mean = c(1.5, 1.0, 0.5),
    seed = 2024
  )

  density_params <- list(
    np_t = ifelse(t_density == "high", 150, 15),
    np_tumor = ifelse(tumor_density == "high", 150, 15),
    np_b = ifelse(t_density == "high", 150, 15)
  )

  patterns <- generate_spatial_patterns_unmeasured_source(
    structure = structure,
    coefficients = coefficients,
    density_params = density_params,
    unmeasured_strength = unmeasured_strength,
    seed = 2024
  )

  p <- plot_confounder_diagnostics(patterns, pattern_index = 1,
                                   confounder_type = "unmeasured")
  print(p)

  return(invisible(patterns))
}

# Example for compartment confounder
example_compartment_diagnostic <- function(t_density = "high",
                                           tumor_density = "high",
                                           n_compartments = 3,
                                           compartment_effect = 0.5) {

  source("./utils.R")
  source("sim_shade_comparison/comparison_functions.R")

  structure <- create_patient_structure(num_patients = 1, num_images = 1, seed = 2024)

  coefficients <- generate_spatial_coefficients(
    structure = structure,
    num_potentials = 3,
    cohort_mean = c(1.5, 1.0, 0.5),
    seed = 2024
  )

  density_params <- list(
    np_t = ifelse(t_density == "high", 150, 15),
    np_tumor = ifelse(tumor_density == "high", 150, 15),
    np_b = ifelse(t_density == "high", 150, 15)
  )

  patterns <- generate_spatial_patterns_compartments(
    structure = structure,
    coefficients = coefficients,
    density_params = density_params,
    n_compartments = n_compartments,
    compartment_effect = compartment_effect,
    seed = 2024
  )

  p <- plot_confounder_diagnostics(patterns, pattern_index = 1,
                                   confounder_type = "compartment")
  print(p)

  return(invisible(patterns))
}
