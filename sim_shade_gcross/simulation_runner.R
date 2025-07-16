library(spatstat)
library(tidyverse)
library(Matrix)
library(SHADE)
library(cmdstanr)

start <- Sys.time()
# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Create patient metadata with proper indexing
create_patient_structure <- function(num_pts_per_group, num_images, seed = 2024) {
  set.seed(seed)
  
  # Generate images per patient based on imbalance level
  total_patients <- num_pts_per_group * 2
  
  # if(imbalance_level == "balanced") {
  images_per_patient <- rep(num_images, total_patients)
  # } else if(imbalance_level == "moderate") {
  #   images_per_patient <- sample(c(1,2,3,4,5), total_patients, replace = TRUE,
  #                                prob = c(0.2, 0.15, 0.3, 0.15, 0.2))
  # } else { # severe
  #   images_per_patient <- sample(c(1,2,3,6), total_patients, replace = TRUE,
  #                                prob = c(0.5, 0.3, 0.15, 0.05))
  # }
  
  # Create patient metadata
  patients_df <- tibble(
    patient_id = 1:total_patients,
    group = c(rep("responder", num_pts_per_group), rep("nonresponder", num_pts_per_group)),
    group_idx = c(1:num_pts_per_group, 1:num_pts_per_group),
    images_count = images_per_patient
  )
  
  # Create image metadata
  total_images <- sum(images_per_patient)
  images_df <- tibble(
    image_id = 1:total_images,
    patient_id = rep(patients_df$patient_id, patients_df$images_count),
    group = rep(patients_df$group, patients_df$images_count),
    image_within_patient = sequence(patients_df$images_count)
  )
  
  return(list(
    patients_df = patients_df,
    images_df = images_df,
    total_patients = total_patients,
    total_images = total_images
  ))
}

#' Generate hierarchical spatial coefficients
generate_spatial_coefficients <- function(structure, num_potentials, seed = 2024) {
  set.seed(seed)
  
  num_responders <- sum(structure$patients_df$group == "responder")
  num_nonresponders <- sum(structure$patients_df$group == "nonresponder")
  
  # Group-level means (fixed effects)
  group_means <- matrix(0, nrow = 2, ncol = num_potentials)
  rownames(group_means) <- c("responder", "nonresponder")
  colnames(group_means) <- paste0("potential_", 1:num_potentials)
  
  # Responders: strong clustering (negative coefficients, decreasing with distance)
  group_means["responder", ] <- c(-1.5, -1, -0.5)  # Short, medium, long range
  
  # Non-responders: weak repulsion (positive coefficients, decreasing with distance)  
  group_means["nonresponder", ] <- c(1.5, 1.0, 0.5)  # Short, medium, long range
  
  # Individual-level effects (random effects around group means)
  sigma_individual <- 0.1
  
  individual_effects <- array(0, dim = c(structure$total_patients, num_potentials))
  rownames(individual_effects) <- paste0("patient_", 1:structure$total_patients)
  colnames(individual_effects) <- paste0("potential_", 1:num_potentials)
  
  for(i in 1:structure$total_patients) {
    patient_group <- structure$patients_df$group[i]
    group_mean_vec <- group_means[patient_group, ]
    individual_effects[i, ] <- rnorm(num_potentials, mean = group_mean_vec, sd = sigma_individual)
  }
  
  # Image-level effects (random effects around individual means)
  sigma_image <- 0.1
  
  image_effects <- array(0, dim = c(structure$total_images, num_potentials))
  rownames(image_effects) <- paste0("image_", 1:structure$total_images)
  colnames(image_effects) <- paste0("potential_", 1:num_potentials)
  
  for(i in 1:structure$total_images) {
    patient_id <- structure$images_df$patient_id[i]
    individual_mean_vec <- individual_effects[patient_id, ]
    image_effects[i, ] <- rnorm(num_potentials, mean = individual_mean_vec, sd = sigma_image)
  }
  
  return(list(
    group_means = group_means,
    individual_effects = individual_effects,
    image_effects = image_effects
  ))
}

#' Generate spatial patterns for all images
generate_spatial_patterns <- function(structure, coefficients, density_params, size_im = 1500, seed = 2024) {
  set.seed(seed)
  
  W <- owin(c(0, size_im), c(0, size_im))
  area <- size_im^2
  
  # Create RBF potentials
  potentials <- make_rbfs(n_basis_functions = 3, max_dist = 75, basis_function_sigma = 15)
  
  patterns <- vector("list", structure$total_images)
  
  for(i in 1:structure$total_images) {
    if(i %% 5 == 0) cat("Generating pattern", i, "of", structure$total_images, "\n")
    
    # Get image-specific parameters
    image_coeffs <- coefficients$image_effects[i, ]
    
    # Generate source cells (T cells and B cells)
    pat_sources <- rmpoispp(lambda = c(density_params$np_t/area, density_params$np_b/area), win = W)
    
    # Generate tumor cell density based on T cells only
    if(sum(pat_sources$marks == 1) > 0) {
      t_cells <- unmark(subset(pat_sources, marks == 1))
      
      # Create custom kernel using image-specific coefficients
      custom_kernel <- Vectorize(function(x, y) {
        d <- sqrt(x^2 + y^2)
        kernel_value <- sum(sapply(seq_along(image_coeffs), function(k) {
          image_coeffs[k] * potentials[[k]](d)
        }))
        return(kernel_value)
      })
      
      # Generate density surface
      dens <- smooth_density_fft(t_cells, custom_kernel, resolution = 128)
      
      # Generate tumor cells
      lambda_integral <- sum(exp(dens$v)) * (dens$xstep * dens$ystep)
      
      if(lambda_integral > 0) {
        beta0 <- log(density_params$np_tumor / lambda_integral)
        tumor_cells <- rpoispp(lambda = exp(dens + beta0))
      } else {
        tumor_cells <- rpoispp(lambda = density_params$np_tumor/area, win = W)
        beta0 <- log(density_params$np_tumor/area)
      }
    } else {
      # No T cells - uniform tumor distribution
      tumor_cells <- rpoispp(lambda = density_params$np_tumor/area, win = W)
      beta0 <- log(density_params$np_tumor/area)
    }
    
    # Set tumor cell marks
    marks(tumor_cells) <- factor(3)
    # Combine all cell types
    combined_pattern <- superimpose(pat_sources, tumor_cells)
    
    patterns[[i]] <- list(
      pattern = combined_pattern,
      beta0 = beta0,
      image_coeffs = image_coeffs,
      patient_id = structure$images_df$patient_id[i],
      group = structure$images_df$group[i]
    )
  }
  
  return(patterns)
}

#' Run SHADE analysis
run_shade_analysis <- function(patterns, structure, num_potentials) {
  cat("Preparing data for SHADE...\n")
  
  # Prepare coordinates data
  coords <- do.call(rbind, lapply(seq_along(patterns), function(i) {
    pat <- patterns[[i]]$pattern
    
    # pat <- subset(pat,marks %in% c(1,3))
    
    tibble(
      x = pat$x,
      y = pat$y,
      type = marks(pat),
      image_id = paste0("img_", i)
    )
  }))
  
  coords %>%
    filter(image_id == "img_3") %>%
    ggplot(aes(x,y)) +
    geom_point() +
    facet_wrap(~type)
  
  # Create patient metadata for SHADE
  image_ids <- unique(coords$image_id)
  coords$image_id <- factor(coords$image_id,level = image_ids)
  patient_metadata <- tibble(
    Spot = image_ids,
    Patient = paste0("pt_", structure$images_df$patient_id),
    Group = structure$images_df$group
  )
  # undebug(prepare_spatial_model_data)
  # undebug(make_data)
  # undebug(make_dispersions)
  # Prepare model data
  prep <- prepare_spatial_model_data(
    x = coords$x,
    y = coords$y,
    cell_type = coords$type,
    image_id = coords$image_id,
    patient_metadata = patient_metadata,
    type_idx = 3,  # Tumor cells are type 3
    n_dummy = 1e3,
    n_basis_functions = num_potentials
  )
  stan_data <- prep$stan_data
  
  # id <- stan_data$sample_id
  # x_cells <- as.matrix(prep$sparse_matrix)[id == i,]
  # is_cell <- stan_data$is_cell[id == i]
  # oset <- stan_data$oset[id == i]
  # summary(glm(y ~ . - 1, data = as.data.frame(cbind(y=is_cell,x_cells)),offset = oset,family = "binomial"))
  stan_data$scale_sigma_betas <- c(10,10,10)
  stan_data$scale_sigmas <- 10
  cat("Fitting SHADE model...\n")
  shade_fit <- run_SHADE_model(stan_data, 
                               method = "variational", 
                               draws = 1e3,
                               # iter_warmup = 500,
                               # iter_sampling = 500,
                               refresh = 20,
                               # chains = 1,
                               threads = 2)
  
  return(list(fit = shade_fit, prep = prep))
}

simple_shade <- function(pattern) {
  tryCatch({
    # SHADE model parameters (adjust as needed)
    n_dummy <- 1000  # number of dummy points
    type <- "3"
    
    Q <- make_quadrature(pattern,n_dummy = n_dummy)
    # plot(Q)
    
    offset <- log(spatstat.geom::intensity(Q$dummy)) |>
      tibble::enframe() |>
      dplyr::right_join(tibble::tibble(name = spatstat.geom::marks(Q)), by = "name") |>
      dplyr::filter(name == type) |>
      dplyr::pull(value)
      
    data_list <- make_data(Q, potentials, type, verbose = FALSE)
      
    x_cells <- as.matrix(data_list$data)[,-1]
    is_cell <- data_list$response
      
    Q_grid <- make_quadrature(pattern,n_dummy = 50,dist="grid")
    plot(Q_grid)
    data_list_grid <- make_data(Q_grid, potentials, type, verbose = FALSE)
      
    X_new <- as.matrix(data_list_grid$data)[,-1]
    resp_new <- data_list_grid$response
    X_new <- X_new[resp_new == 0,]
    N_new <- nrow(X_new)
    offset_new <- rep(unique(offset),N_new)
      
    data_stan <- list(
      N = nrow(x_cells),
      K = ncol(x_cells),
      y = is_cell,
      X = x_cells,
      oset = -offset,
      N_new = N_new,
      X_new = X_new,
      offset_new = offset_new
    )
      
    mod <- cmdstan_model("sim_shade_gcross/simple_shade.stan")
    fit <- mod$variational(data_stan,draws=1e3,refresh=100)
    
    list(fit=fit,data_stan=data_stan)
  }, error = function(e) {
    print(e)
  })
}

#' Run G-cross analysis
run_gcross_analysis <- function(patterns, structure) {
  cat("Running G-cross analysis...\n")
  
  # Image-level G-cross results
  image_results <- lapply(seq_along(patterns), function(i) {
    if(i %% 10 == 0) {
      print(i)
    }

    out <- tryCatch({

      pat <- patterns[[i]]$pattern
      true_group <- patterns[[i]]$group
      plot(pat)
      # Check cell counts
      n_t <- sum(pat$marks == "1")
      n_tumor <- sum(pat$marks == "3")
      
      if(n_t < 3 || n_tumor < 3) {
        stop("not enough cells!")
      }
      num_env <- 10
      # maybe do this a bunch of times?
      detects <- lapply(1:num_env,\(j) {
        
      
      # Calculate G-cross with envelope
      env <- envelope(pat, Gcross, i = "1", j = "3", nsim = 39, 
                      correction = "border", verbose = FALSE)
      # plot(env)
      # Test for significant deviations at short range
      # Define target distances
      target_r <- c(20, 40, 60)
      
      # Find nearest available r indices
      target_indices <- sapply(target_r, function(rval) which.min(abs(env$r - rval)))
      
      # Check for clustering or repulsion at those distances
      clustering_detected <- sum(env$obs[target_indices] > env$hi[target_indices], na.rm = TRUE)
      repulsion_detected  <- sum(env$obs[target_indices] < env$lo[target_indices], na.rm = TRUE)
      
      # Check if detection matches expected pattern
      if(true_group == "nonresponder") {
        return(clustering_detected)  # Should detect clustering
      } else {
        return(repulsion_detected)   # Should detect repulsion
      }
      })
      
      target_r <- c(20, 40, 60)
      
      # Find nearest available r indices
      mean_detect <- mean(unlist(detects))
      gx <- Gcross(pat,i = "1",j = "3")
      target_indices <- sapply(target_r, function(rval) which.min(abs(gx$r - rval)))
      
      out <- list(mean_detect=mean_detect,gx=gx$km[target_indices])
    }, error = function(e) {
      out <- list(mean_detect=NA,gx=NA)
    })
    return(out)
  })
  
  return(image_results)
}

# ============================================================================
# MAIN SIMULATION SCRIPT
# ============================================================================

# Set environment and get parameters
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  path <- "./sim_shade_gcross/data/"
  sim_idx <- 1
} else {
  path <- "./sim_shade_gcross/data/"
  args <- commandArgs(trailingOnly=TRUE)
  sim_idx <- as.numeric(args[1])
  cat("Running simulation", sim_idx, "\n")
}

file_out <- paste0(path,"sim_",sim_idx,".rds")


if (file.exists(file_out) & SYSTEM_ENV == "HPC") {
  cat("Simulation", sim_idx, "already completed. Skipping.\n")
  quit(save="no")
}

grid <- expand.grid(
  # imbalance_level = c("balanced", "moderate", "severe"),
  t_density = c("high", "low"), 
  tumor_density = c("high", "low"),
  num_images = c(1,2,3),
  sim_rep = 1:30
)

# Extract current condition parameters
condition <- list(
  # imbalance_level = as.character(grid$imbalance_level[sim_idx]),
  num_images = as.character(grid$num_images[sim_idx]),
  t_density = as.character(grid$t_density[sim_idx]),
  tumor_density = as.character(grid$tumor_density[sim_idx]),
  sim_rep = grid$sim_rep[sim_idx]
)

cat("Condition:", condition$t_density, condition$tumor_density, condition$num_images,
    "Rep:", condition$sim_rep, "\n")

# Fixed parameters
num_pts_per_group <- 20
num_potentials <- 3
size_im <- 1500

# Density parameters
density_params <- list(
  np_t = ifelse(condition$t_density == "high", 150, 15),
  np_tumor = ifelse(condition$tumor_density == "high", 150, 15),
  np_b = ifelse(condition$t_density == "high", 150, 15)  # B cells same as T cells
)

# ============================================================================
# RUN SIMULATION
# ============================================================================

cat("Creating patient structure...\n")
structure <- create_patient_structure(
  num_pts_per_group = num_pts_per_group,
  num_images = grid$num_images[sim_idx],
  # imbalance_level = condition$imbalance_level,
  seed = 2024 + sim_idx
)

cat("Generating spatial coefficients...\n")
coefficients <- generate_spatial_coefficients(
  structure = structure,
  num_potentials = num_potentials,
  seed = 2024 + condition$sim_rep
)

cat("Generating spatial patterns...\n")
patterns <- generate_spatial_patterns(
  structure = structure,
  coefficients = coefficients,
  density_params = density_params,
  size_im = size_im,
  seed = 2024 + sim_idx + 1000
)

cat("Running SHADE analysis...\n")
shade_results <- run_shade_analysis(patterns, structure, num_potentials)

cat("Running G-cross analysis...\n")
gcross_results <- run_gcross_analysis(patterns, structure)

# ============================================================================
# EXTRACT RESULTS AND CALCULATE METRICS
# ============================================================================

# Extract SHADE results
library(posterior)
shade_draws <- as_draws_rvars(shade_results$fit$draws())
potentials <- make_rbfs(max_dist = 75, n_basis_functions = 3, basis_function_sigma = 15)
# Calculate SHADE image-level power
simple_shade_image_results <- sapply(1:structure$total_images, function(i) {
  true_group <- structure$images_df$group[i]
  pattern <- patterns[[i]]$pattern
  # plot(pattern)
  # plot(split(pattern))
  
  tryCatch({
    
  
  results <- simple_shade(pattern)
  shade_draws <- as_draws_rvars(results$fit$draws())
  
  # Get image-level coefficients for all 3 potentials
  beta <- as.vector(shade_draws$beta[1:3])
  beta_true <- coefficients$image_effects[i,]
  
  # coefficients$individual_effects[22,]
  
  x_seq <- seq(0,80,1)
  x_des <- lapply(potentials,\(pot) pot(x_seq)) %>% do.call(cbind,.)

  lp <- as.vector(x_des %*% beta)
  lp_true <- as.vector(x_des %*% beta_true)
  
  alpha <- 0.05  # for 95% simultaneous bands
  
  # Assume lp is an rvar data.frame where each column (except x) is an rvar
  
  # Step 1: Extract mean and SD per x
  lp %>%
    as.matrix() %>%
    as.data.frame()  %>%
    mutate(x = x_seq) %>%
    mutate(
      mn = as.vector(E(V1)),
      sd = as.vector(sd(V1))
    ) -> df_summary
  
  # Step 2: Standardize residuals across samples
  lp_std <- lp %>%
    as.matrix() %>%
    as.data.frame() %>%
    mutate(across(
      everything(),
      ~ (. - mean(.)) / sd(.)
    ))
  
  # Step 3: Find maximum absolute deviation across x for each posterior sample (column-wise)
  max_dev <- max(abs(lp_std$V1)) %>%
    as.matrix() %>%
    as.data.frame()
  
  # Step 4: Find critical value (z_score_band) for simultaneous coverage
  z_score_band <- max_dev %>%
    pivot_longer(everything(), names_to = "sample", values_to = "max_dev") %>%
    pull(max_dev) %>%
    quantile(probs = 1 - alpha)
  
  # Step 5: Construct simultaneous lower and upper bands
  df_summ <- df_summary %>%
    mutate(
      lo_simul = mn - z_score_band * sd,
      hi_simul = mn + z_score_band * sd,
      lo_pw = as.vector(quantile(V1,probs=alpha/2)),
      hi_pw = as.vector(quantile(V1,probs=1-(alpha/2))),
      lp_true = lp_true
    )
  
  # df_summ %>%
  #   head()
  df_summ %>%
    ggplot(aes(x)) +
    geom_line(aes(y=mn)) +
    geom_line(aes(y=lp_true),color="steelblue") +
    geom_ribbon(aes(ymin=lo_simul,ymax=hi_simul),alpha=0.2) +
    geom_ribbon(aes(ymin=lo_pw,ymax=hi_pw),alpha=0.2,fill="red") +
    geom_hline(yintercept=0,color="red",linetype="dashed")
  
  # Specify target x values
  target_x <- c(20, 40, 60)
  
  # Find indices of x values in df_summ closest to each target
  target_indices <- sapply(target_x, function(x_val) which.min(abs(df_summ$x - x_val)))
  
  # Subset the relevant rows
  target_rows <- df_summ[target_indices, ]
  
  # Determine direction
  is_negative <- with(target_rows, hi_pw < 0)
  is_positive <- with(target_rows, lo_pw > 0)
  
  # Assume `responder_status` is a variable: TRUE for responder, FALSE for non-responder
  # (you'll need to define this in context)
  if (true_group == "responder") {
    sum_pass <- sum(is_negative)
  } else {
    sum_pass <- sum(is_positive)
  }
  
  return(sum_pass)
  
  },
  error=\(e) {
    print(i)
    print(e)
    return(NULL)
  })
})

# Calculate SHADE image-level power
shade_image_results <- sapply(1:structure$total_images, function(i) {
  true_group <- structure$images_df$group[i]
  
  # Get image-level coefficients for all 3 potentials
  beta <- as.vector(shade_draws$beta_local[i,2:4])
  beta_true <- coefficients$image_effects[i,]
  
  coefficients$individual_effects[22,]
  
  x_seq <- seq(0,80,1)
  x_des <- lapply(potentials,\(pot) pot(x_seq)) %>% do.call(cbind,.)
  
  lp <- as.vector(x_des %*% beta)
  lp_true <- as.vector(x_des %*% beta_true)
  
  alpha <- 0.05  # for 95% simultaneous bands
  
  # Assume lp is an rvar data.frame where each column (except x) is an rvar
  
  # Step 1: Extract mean and SD per x
  lp %>%
    as.matrix() %>%
    as.data.frame()  %>%
    mutate(x = x_seq) %>%
    mutate(
      mn = as.vector(E(V1)),
      sd = as.vector(sd(V1))
    ) -> df_summary
  
  # Step 2: Standardize residuals across samples
  lp_std <- lp %>%
    as.matrix() %>%
    as.data.frame() %>%
    mutate(across(
      everything(),
      ~ (. - mean(.)) / sd(.)
    ))
  
  # Step 3: Find maximum absolute deviation across x for each posterior sample (column-wise)
  max_dev <- max(abs(lp_std$V1)) %>%
    as.matrix() %>%
    as.data.frame()
  
  # Step 4: Find critical value (z_score_band) for simultaneous coverage
  z_score_band <- max_dev %>%
    pivot_longer(everything(), names_to = "sample", values_to = "max_dev") %>%
    pull(max_dev) %>%
    quantile(probs = 1 - alpha)
  
  # Step 5: Construct simultaneous lower and upper bands
  df_summ <- df_summary %>%
    mutate(
      lo_simul = mn - z_score_band * sd,
      hi_simul = mn + z_score_band * sd,
      lo_pw = as.vector(quantile(V1,probs=alpha/2)),
      hi_pw = as.vector(quantile(V1,probs=1-(alpha/2))),
      lp_true = lp_true
    )
  
  # df_summ %>%
  #   head()
  df_summ %>%
    ggplot(aes(x)) +
    geom_line(aes(y=mn)) +
    geom_line(aes(y=lp_true),color="steelblue") +
    geom_ribbon(aes(ymin=lo_simul,ymax=hi_simul),alpha=0.2) +
    geom_ribbon(aes(ymin=lo_pw,ymax=hi_pw),alpha=0.2,fill="red") +
    geom_hline(yintercept=0,color="red",linetype="dashed")
  
  # Specify target x values
  target_x <- c(20, 40, 60)
  
  # Find indices of x values in df_summ closest to each target
  target_indices <- sapply(target_x, function(x_val) which.min(abs(df_summ$x - x_val)))
  
  # Subset the relevant rows
  target_rows <- df_summ[target_indices, ]
  
  # Determine direction
  is_negative <- with(target_rows, hi_pw < 0)
  is_positive <- with(target_rows, lo_pw > 0)
  
  # Assume `responder_status` is a variable: TRUE for responder, FALSE for non-responder
  # (you'll need to define this in context)
  if (true_group == "responder") {
    sum_pass <- sum(is_negative)
  } else {
    sum_pass <- sum(is_positive)
  }
  
  return(sum_pass)
})

# Calculate SHADE group-level significance
# group_effect <- shade_draws$beta_global[2:4, 2] - shade_draws$beta_global[2:4, 1]
# group_ci_lower <- as.vector(quantile(group_effect, 0.025))
# group_ci_upper <- as.vector(quantile(group_effect, 0.975))
# 
# all(group_ci_lower > 0)

# Calculate performance metrics
shade_image_power <- mean(shade_image_results, na.rm = TRUE)
gx_detect <- unlist(lapply(gcross_results,\(o) o$mean_detect))
gcross_image_power <- mean(gx_detect, na.rm = TRUE)
simple_shade_image_power <- mean(unlist(compact(simple_shade_image_results)))

out <- list(shade_image_power=shade_image_power,
            simple_shade_image_power=simple_shade_image_power,
            gcross_image_power=gcross_image_power)

saveRDS(out,file_out)

end <- Sys.time()

end - start
# gx <- lapply(gcross_results,\(o) o$gx)
# gx <- do.call(rbind,gx) %>%
#   as.data.frame() %>%
#   mutate(image_id=structure$images_df$image_id,
#          patient_id=structure$images_df$patient_id,
#          group = structure$images_df$group)

#### extra plots
theme_set(theme_bw(base_size=14, base_family='Helvetica')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

fsave <- \(fname,height=5,width=8) {
  ggsave(paste0(figures_folder,fname,".pdf"),device=cairo_pdf, height=height, width=width, units="in")
}
figures_folder <- "./sim_shade_gcross/gx_figures/"
i <- 1
pat <- patterns[[i]]$pattern
marks(pat) <- fct_recode(marks(pat),"Tumor cells"="1","B cells"="2","T cells"="3")
plot(pat)

p1 <- pat %>%
  as.data.frame() %>%
  rename(Type=marks) %>%
  mutate(Type = fct_recode(Type,"Tumor cells"="1","B cells"="2","T cells"="3")) %>%
  filter(Type != "B cells") %>%
  # mutate(is_B = Type != "B cells") %>%
  ggplot(aes(x,y,color=Type,shape=Type)) +
  geom_point(size=1.5) +
  coord_fixed() +
  labs(x="X (microns)",y="Y (microns)")
p1
fsave("pat_example")

pat <- patterns[[i]]$pattern
marks(pat) <- fct_recode(marks(pat),"Tumor"="1","B"="2","T"="3")
env <- envelope(pat, Gcross, i = "Tumor", j = "T", nsim = 39, 
                correction = "border", verbose = FALSE)
plot(env)

p2 <- env %>%
  as.data.frame() %>%
  filter(r < 80) %>%
  ggplot(aes(x=r)) +
  geom_line(aes(y=obs)) +
  geom_line(aes(y=theo),color="red",linetype="dashed") +
  geom_ribbon(aes(ymin=lo,ymax=hi),color="gray50",alpha=0.2,linetype=0) +
  labs(x="Distance (microns)",y="Gcross")
p2
fsave("gx_example")

# Get image-level coefficients for all 3 potentials
beta <- as.vector(shade_draws$beta_local[i,2:4])
beta_true <- coefficients$image_effects[i,]

# coefficients$individual_effects[22,]

x_seq <- seq(0,80,1)
x_des <- lapply(potentials,\(pot) pot(x_seq)) %>% do.call(cbind,.)

lp <- as.vector(x_des %*% beta)
lp_true <- as.vector(x_des %*% beta_true)

alpha <- 0.05  # for 95% simultaneous bands

# Step 1: Extract mean and SD per x
lp %>%
  as.matrix() %>%
  as.data.frame()  %>%
  mutate(x = x_seq) %>%
  mutate(
    mn = as.vector(E(V1)),
    sd = as.vector(sd(V1))
  ) -> df_summary

# Step 5: Construct simultaneous lower and upper bands
df_summ <- df_summary %>%
  mutate(
    lo_pw = as.vector(quantile(V1,probs=alpha/2)),
    hi_pw = as.vector(quantile(V1,probs=1-(alpha/2))),
    lp_true = lp_true
  )

pat <- patterns[[i]]$pattern
results <- simple_shade(pat)
simple_shade_draws <- as_draws_rvars(results$fit$draws())

beta_simple <- as.vector(simple_shade_draws$beta[1:3])

lp_simple <- as.vector(x_des %*% beta_simple)

lp %>%
  as.matrix() %>%
  as.data.frame()  %>%
  mutate(x = x_seq) %>%
  mutate(
    mn = as.vector(E(V1)),
    sd = as.vector(sd(V1))
  ) -> df_summary


df_summ <- df_summ %>%
  mutate(
    mn_simple = as.vector(E(lp_simple)),
    lo_simple = as.vector(quantile(lp_simple,probs=alpha/2)),
    hi_simple = as.vector(quantile(lp_simple,probs=1-(alpha/2))),
  )

# Create a long-format data frame for the lines
df_lines <- df_summ %>%
  select(x, mn, mn_simple, lp_true) %>%
  pivot_longer(cols = c(mn, mn_simple, lp_true), names_to = "type", values_to = "value")

# Create a long-format data frame for the ribbons
df_ribbons <- df_summ %>%
  transmute(x,
            type = "mn",
            ymin = lo_pw,
            ymax = hi_pw) %>%
  bind_rows(
    df_summ %>%
      transmute(x,
                type = "mn_simple",
                ymin = lo_simple,
                ymax = hi_simple)
  )

# Plot
p3 <- ggplot() +
  # Ribbons
  geom_ribbon(data = df_ribbons, aes(x = x, ymin = ymin, ymax = ymax, fill = type), alpha = 0.2) +
  # Lines
  geom_line(data = df_lines, aes(x = x, y = value, color = type), size = 1) +
  geom_hline(yintercept=0,color="black",linetype="dashed") +
  scale_color_manual(values = c(
    mn = "#e6194b",
    mn_simple = "#3cb44b",
    lp_true = "#4363d8"
  ),
  labels = c(
    mn = "SHADE",
    mn_simple = "Flat",
    lp_true = "True"
  )) +
  scale_fill_manual(values = c(
    mn = "#e6194b",
    mn_simple = "#3cb44b"
  )) +
  guides(fill="none") +
  labs(x = "Distance (microns)", y = "Estimated SIC", color = "Method")
p3
fsave("sic_comparison")

library(patchwork)

(p1 + p2) / (p3 + plot_spacer()) + plot_annotation(tag_levels='a')
fsave("three_plots")
