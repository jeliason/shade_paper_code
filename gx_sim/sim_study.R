# Bayesian Hierarchical Spatial Model vs G-cross Function Comparison
# Simulation Study for Statistical Power and Type I Error

library(spatstat)
library(tidyverse)
library(parallel)
library(SHADE)
library(posterior)
library(cmdstanr)

# simulate with SHADE, two separate groups. Unbalanced data, should
# Want to see if we can identify that the SICs and Gcross between the two groups are in fact different
# Also want to see how well SICs and Gcross sparse, unbalanced images are estimated (pooling should help?)
# don't do patient-level comparison, too many layers

# Set simulation parameters
set.seed(12345)
n_simulations <- 100
lambda_A <- 20  # intensity for type A points
lambda_B <- 20  # intensity for type B points

# Create observation window
W <- owin(c(0,1000),c(0,1000))

num_pot <- 3
potentials <- make_rbfs(
  n_basis_functions = num_pot,
  max_dist = 75,
  basis_function_sigma = 15
)

beta <- c(0.2,-0.3,0.1)
create_shade_pattern <- function(beta) {
  # predict the last type from all of the others
  pat <- rpoispp(lambda=lambda_A/area(W),win = W)
  marks(pat) <- factor("A")
  # plot(pat)
  custom_kernel <- Vectorize(function(x, y) {
    d <- sqrt(x^2 + y^2)  # Compute distance
    
    # Compute weighted sum of basis functions
    kernel_value <- sum(sapply(seq_along(beta), function(i) {
      beta[i] * potentials[[i]](d)
    }))
    
    return(kernel_value)
  })
  
  subs <- unmark(pat)
  dens <- smooth_density_fft(subs, custom_kernel, resolution = 128)
  # plot(dens)
  
  lambda_integral <- sum(exp(dens$v)) * (dens$xstep * dens$ystep)  # Approximate integral using grid summation
  
  # Compute the necessary beta0[i] to achieve expected np points
  beta0 <- log(lambda_B / lambda_integral)
  pat2 <- rpoispp(lambda = exp(dens + beta0))
  marks(pat2) <- factor("B")

  pat <- superimpose(pat,pat2)
  # list(pat=pat,beta0=beta0)
}

# plot_rbfs <- function(potentials, x_range = c(0, 1), n_points = 100) {
#   # Create sequence of x values
#   x_vals <- seq(from = x_range[1], to = x_range[2], length.out = n_points)
#   
#   # Evaluate each RBF function over the x values
#   rbf_data <- map_dfr(names(potentials), function(rbf_name) {
#     y_vals <- map_dbl(x_vals, potentials[[rbf_name]])
#     
#     tibble(
#       x = x_vals,
#       y = y_vals,
#       rbf = rbf_name
#     )
#   })
#   
#   # Create the plot
#   ggplot(rbf_data, aes(x = x, y = y, color = rbf)) +
#     geom_line(linewidth = 1) +
#     labs(
#       title = "Radial Basis Functions",
#       x = "Distance",
#       y = "Potential",
#       color = "RBF"
#     ) +
#     theme_minimal() +
#     theme(
#       legend.position = "right",
#       panel.grid.minor = element_blank()
#     )
# }

# Function to simulate null model (independent processes)
# simulate_null <- function() {
#   A_points <- rpoispp(lambda_A, win = W)
#   B_points <- rpoispp(lambda_B, win = W)
#   marks_A <- rep("A", A_points$n)
#   marks_B <- rep("B", B_points$n)
#   
#   all_x <- c(A_points$x, B_points$x)
#   all_y <- c(A_points$y, B_points$y)
#   all_marks <- c(marks_A, marks_B)
#   
#   ppp(all_x, all_y, window = W, marks = factor(all_marks))
# }

# pp <- simulate_null()
# plot(pp)
# Function to simulate clustering model (Thomas cluster process with saveparents)
simulate_clustering <- function() {
  # Generate Thomas cluster process with parents saved
  cluster_pattern <- rThomas(kappa = lambda_A/3/area(W), scale = 80, mu = 3,
                             win = W, saveparents = TRUE)

  # Extract parent points (type A) and offspring points (type B)
  parents <- as.data.frame(attr(cluster_pattern, "parents"))

  if(nrow(parents) > 0) {
    # Parent points are type A
    A_x <- parents[,1]
    A_y <- parents[,2]
    marks_A <- rep("A", length(A_x))

    # Offspring points are type B
    B_x <- cluster_pattern$x
    B_y <- cluster_pattern$y
    marks_B <- rep("B", length(B_x))

    all_x <- c(A_x, B_x)
    all_y <- c(A_y, B_y)
    all_marks <- c(marks_A, marks_B)

    pat <- ppp(all_x, all_y, window = W, marks = factor(all_marks))
  } else {
    # Fallback to independent processes if no parents generated
    simulate_null()
  }
}

pp <- simulate_clustering()
plot(pp)
# Function to simulate repulsion model using MultiHard
simulate_repulsion <- function() {
  
  # Multitype Strauss:
  beta <- c(0.027,0.008)
  gmma <- matrix(c(0.98,0.36,0.36,0.98),2,2)
  r    <- matrix(c(45,45,45,45),2,2)
  mod08 <- list(cif="straussm",par=list(beta=beta,gamma=gmma,radii=r),
                w=c(0,1000,0,1000))
  nr   <- 1e5
  nv  <- 5000
  ns <- 200
  pattern <- rmh(model=mod08,start=list(n.start=c(lambda_A,lambda_B)),
                 control=list(fixall=TRUE,p=1,ptypes=c(0.75,0.25),
                              nrep=nr,nverb=nv))
  marks(pattern) <- fct_recode(marks(pattern),"A"="1","B"="2")
  pattern
}

# beta <- c(1.0,-0.5,0.5)
# pattern <- create_shade_pattern(beta)
pattern <- simulate_clustering()
# pattern <- simulate_repulsion()
plot(pattern)
r <- seq(0,80,length.out=513)
# Function to perform G-cross analysis
analyze_gcross <- function(pattern) {
  tryCatch({
    # Calculate G-cross function
    g_AB <- Gcross(pattern, i = "A", j = "B", correction = "border",r=r)
    g_BA <- Gcross(pattern, i = "B", j = "A", correction = "border",r=r)
    
    plot(g_AB)
    plot(g_BA)
    # Monte Carlo envelopes
    env_AB <- envelope(pattern, Gcross, i = "A", j = "B", nsim = 99, 
                       correction = "border", verbose = FALSE,global=TRUE,nrank=5,r=r)
    env_BA <- envelope(pattern, Gcross, i = "B", j = "A", nsim = 99, 
                       correction = "border", verbose = FALSE,global=TRUE,nrank=5,r=r)
    plot(env_AB)
    plot(env_BA)
    
    # Test for significant deviation
    sig_AB_clustering <- any(env_AB$obs > env_AB$hi, na.rm = TRUE)
    sig_AB_repulsion <- any(env_AB$obs < env_AB$lo, na.rm = TRUE)
    sig_BA_clustering <- any(env_BA$obs > env_BA$hi, na.rm = TRUE)
    sig_BA_repulsion <- any(env_BA$obs < env_BA$lo, na.rm = TRUE)
    
    list(
      sig_AB_clustering = sig_AB_clustering,
      sig_AB_repulsion = sig_AB_repulsion,
      sig_BA_clustering = sig_BA_clustering,
      sig_BA_repulsion = sig_BA_repulsion,
      g_AB = g_AB,
      g_BA = g_BA,
      env_AB = env_AB,
      env_BA = env_BA
    )
  }, error = function(e) {
    print(e)
    list(
      sig_AB_clustering = FALSE,
      sig_AB_repulsion = FALSE,
      sig_BA_clustering = FALSE,
      sig_BA_repulsion = FALSE,
      error = TRUE
    )
  })
}
pattern <- simulate_repulsion()
plot(pattern)
# SHADE model analysis
analyze_shade <- function(pattern) {
  tryCatch({
    # SHADE model parameters (adjust as needed)
    n_dummy <- 1000  # number of dummy points
    warmup <- 500
    iters <- 500

    type <- "B"
    
    # potentials
    # plot_rbfs(potentials)
    
    Q <- make_quadrature(pattern,n_dummy = n_dummy)
    plot(Q)
    
    get_significance_SHADE <- function(type) {
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
      
      mod <- cmdstan_model("gx_sim/simple_shade.stan")
      fit <- mod$sample(data_stan,iter_warmup = warmup,iter_sampling = iters,chains=1,parallel_chains = 1,refresh=100)
      
      x_seq <- seq(0,80,1)
      x_des <- lapply(potentials,\(pot) pot(x_seq)) %>% do.call(cbind,.)
      
      draws <- as_draws_rvars(fit$draws())
      draws
      beta <- as.vector(draws$beta)
      lp <- as.vector(x_des %*% beta)
  
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
          hi_simul = mn + z_score_band * sd
        )
      
      df_summ %>%
        ggplot(aes(x)) +
        geom_line(aes(y=mn)) +
        geom_ribbon(aes(ymin=lo_simul,ymax=hi_simul),alpha=0.2) +
        geom_hline(yintercept=0,color="red",linetype="dashed")
      
      sig_clustering <- any(df_summ$lo_simul > 0, na.rm = TRUE)
      sig_repulsion <- any(df_summ$hi_simul < 0, na.rm = TRUE)
      
      dummy <- Q_grid$dummy
      dummy <- subset(dummy,marks == type)
      coords <- cbind(x=dummy$x,y=dummy$y)
      lp <- draws$linear_pred_new
      
     coords %>%
        as.data.frame() %>%
       mutate(lp=E(lp)) %>%
       filter(x < 1000,y < 1000,x>0,y>0) %>%
       ggplot() +
       geom_raster(aes(x,y,fill=lp)) +
       geom_point(data=as.data.frame(pattern),aes(x,y,shape=marks,color=marks))
      
      list(sig_clustering = sig_clustering,
           sig_repulsion = sig_repulsion)
    }
    
    AB_sig <- get_significance_SHADE("B")
    BA_sig <- get_significance_SHADE("A")

    
    list(
      sig_AB_clustering = AB_sig$sig_clustering,  
      sig_AB_repulsion = AB_sig$sig_repulsion,  
      sig_BA_clustering = BA_sig$sig_clustering,  
      sig_BA_repulsion = BA_sig$sig_repulsion
    )
    
  }, error = function(e) {
    print(e)
    # Return failed analysis if SHADE fitting fails
    list(
      sig_AB_clustering = FALSE,
      sig_AB_repulsion = FALSE,
      sig_BA_clustering = FALSE,
      sig_BA_repulsion = FALSE,
      error = TRUE,
      error_message = e$message
    )
  })
}
n_simulations <- n_sim <- 10
# Run simulation study
run_simulation_study <- function(scenario, n_sim = n_simulations) {
  
  sim_function <- switch(scenario,
                         "null" = simulate_null,
                         "clustering" = simulate_clustering,
                         "repulsion" = simulate_repulsion)
  
  results_gcross <- list()
  results_shade <- list()
  
  cat("Running", scenario, "scenario...\n")
  
  for(i in 1:n_sim) {
    if(i %% 10 == 0) cat("  Simulation", i, "of", n_sim, "\n")
    
    # Generate pattern
    pattern <- sim_function()
    
    # Analyze with G-cross
    gcross_result <- analyze_gcross(pattern)
    results_gcross[[i]] <- gcross_result
    
    # Analyze with SHADE (placeholder)
    shade_result <- analyze_shade(pattern)
    results_shade[[i]] <- shade_result
  }
  
  list(
    scenario = scenario,
    gcross = results_gcross,
    shade = results_shade
  )
}

# Run all scenarios
scenarios <- c("null", "clustering", "repulsion")
all_results <- list()

for(scenario in scenarios) {
  all_results[[scenario]] <- run_simulation_study(scenario)
}

# Function to calculate performance metrics
calculate_metrics <- function(results, scenario) {
  
  n_valid <- length(results)
  
  # Extract detection results
  gcross_AB_clustering <- sapply(results$gcross, function(x) x$sig_AB_clustering)
  gcross_AB_repulsion <- sapply(results$gcross, function(x) x$sig_AB_repulsion)
  gcross_BA_clustering <- sapply(results$gcross, function(x) x$sig_BA_clustering)
  gcross_BA_repulsion <- sapply(results$gcross, function(x) x$sig_BA_repulsion)
  
  shade_AB_clustering <- sapply(results$shade, function(x) x$sig_AB_clustering)
  shade_AB_repulsion <- sapply(results$shade, function(x) x$sig_AB_repulsion)
  shade_BA_clustering <- sapply(results$shade, function(x) x$sig_BA_clustering)
  shade_BA_repulsion <- sapply(results$shade, function(x) x$sig_BA_repulsion)
  
  # Calculate performance based on scenario
  if(scenario == "null") {
    # Type I error (false positives)
    gcross_type1 <- mean(gcross_AB_clustering | gcross_AB_repulsion | 
                           gcross_BA_clustering | gcross_BA_repulsion, na.rm = TRUE)
    shade_type1 <- mean(shade_AB_clustering | shade_AB_repulsion | 
                          shade_BA_clustering | shade_BA_repulsion, na.rm = TRUE)
    
    data.frame(
      Method = c("G-cross", "SHADE"),
      TypeI_Error = c(gcross_type1, shade_type1),
      Power_Clustering = c(NA, NA),
      Power_Repulsion = c(NA, NA)
    )
    
  } else if(scenario == "clustering") {
    # Statistical power for detecting clustering
    gcross_power <- mean(gcross_AB_clustering | gcross_BA_clustering, na.rm = TRUE)
    shade_power <- mean(shade_AB_clustering | shade_BA_clustering, na.rm = TRUE)
    
    # Type I error for repulsion detection
    gcross_type1 <- mean(gcross_AB_repulsion | gcross_BA_repulsion, na.rm = TRUE)
    shade_type1 <- mean(shade_AB_repulsion | shade_BA_repulsion, na.rm = TRUE)
    
    data.frame(
      Method = c("G-cross", "SHADE"),
      TypeI_Error = c(gcross_type1, shade_type1),
      Power_Clustering = c(gcross_power, shade_power),
      Power_Repulsion = c(NA, NA)
    )
    
  } else if(scenario == "repulsion") {
    # Statistical power for detecting repulsion
    gcross_power <- mean(gcross_AB_repulsion | gcross_BA_repulsion, na.rm = TRUE)
    shade_power <- mean(shade_AB_repulsion | shade_BA_repulsion, na.rm = TRUE)
    
    # Type I error for clustering detection
    gcross_type1 <- mean(gcross_AB_clustering | gcross_BA_clustering, na.rm = TRUE)
    shade_type1 <- mean(shade_AB_clustering | shade_BA_clustering, na.rm = TRUE)
    
    data.frame(
      Method = c("G-cross", "SHADE"),
      TypeI_Error = c(gcross_type1, shade_type1),
      Power_Clustering = c(NA, NA),
      Power_Repulsion = c(gcross_power, shade_power)
    )
  }
}

# Calculate and display results
performance_results <- list()

for(scenario in scenarios) {
  cat("\n=== Performance Metrics for", toupper(scenario), "Scenario ===\n")
  perf <- calculate_metrics(all_results[[scenario]], scenario)
  performance_results[[scenario]] <- perf
  print(perf)
}

# Summary table
summary_table <- do.call(rbind, lapply(names(performance_results), function(s) {
  df <- performance_results[[s]]
  df$Scenario <- s
  df
}))

write.csv(summary_table, "simulation_results.csv", row.names = FALSE)

# Create visualization of G-cross functions for example patterns
par(mfrow = c(2, 3))

for(scenario in scenarios) {
  # Get first valid result for visualization
  example_idx <- 1
  while(is.null(all_results[[scenario]]$gcross[[example_idx]]$g_AB) && 
        example_idx <= length(all_results[[scenario]]$gcross)) {
    example_idx <- example_idx + 1
  }
  
  if(example_idx <= length(all_results[[scenario]]$gcross)) {
    gcross_result <- all_results[[scenario]]$gcross[[example_idx]]
    
    # Plot G_AB
    plot(gcross_result$env_AB, main = paste("G-cross A→B:", scenario))
    
    # Plot G_BA  
    plot(gcross_result$env_BA, main = paste("G-cross B→A:", scenario))
  }
}

cat("\n=== Simulation Study Complete ===\n")
cat("Results saved to simulation_results.csv\n")
cat("Performance metrics calculated for", n_simulations, "simulations per scenario\n")
cat("SHADE analysis placeholder implemented - replace with actual SHADE model\n")