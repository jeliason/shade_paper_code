library(tidyverse)
library(latex2exp)
library(spatstat)
library(SHADE)

theme_set(theme_bw(base_size=14, base_family='Helvetica')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

fsave <- \(fname,height=5,width=8) {
  ggsave(paste0(figures_folder,fname,".pdf"),device=cairo_pdf, height=height, width=width, units="in")
}
figures_folder <- "./figures/demo_plots/"



x_seq <- seq(0,100,0.1)
weights <- c(0.8,-0.4,0.3)
rbfs <- make_rbfs(n_basis_functions = 3, max_dist = 75, basis_function_sigma = 15)
pot_mat <- lapply(1:length(rbfs),\(i) rbfs[[i]](x_seq)*weights[i]) %>% do.call(cbind,.) %>% rowSums()
tibble(x=x_seq,y=pot_mat) %>%
  ggplot(aes(x,y)) +
  geom_line(color="red",linewidth=1.5) +
  geom_hline(yintercept=0,linetype="dashed") +
  geom_segment(x=75,y=-0.4,xend=75,yend=0.28,color="blue",linetype="dotted",linewidth=1.5) +
  annotate(
    'segment',
    x = 50, # Play around with the coordinates until you're satisfied
    y = 0.48,
    yend = 0.28,
    xend = 75,
    linewidth = 1,
    arrow = arrow(length = unit(0.5, 'cm'))
  ) +
  annotate(
    'text',
    x = 60, # Play around with the coordinates until you're satisfied
    y = 0.5,
    size=4,
    fontface="bold",
    label = expression(paste("Type ", italic("i"), " cell presence ", "\u2192", "  increase of 0.35 in log-intensity of type ", italic("j"), " cells at 75 microns")),
  ) +
  labs(x="Distance (microns)",
       y=expression(paste("Log-intensity of type ",italic("j"))))
fsave("sic_explainer")

# svc_weights <- c(0.3,0.5,-0.1)
# svc <- lapply(1:length(rbfs),\(i) rbfs[[i]](x_seq)*svc_weights[i]) %>% do.call(cbind,.) %>% rowSums()
# p_tb <- tibble(x=x_seq,SIC=pot_mat,SVC=svc) %>%
#   mutate(`Linear Predictor`=SIC*SVC) %>%
#   pivot_longer(-x)
# 
# p_tb %>%
#   ggplot(aes(x,value,color=name)) +
#   geom_line() +
#   geom_ribbon(data=subset(p_tb,name == "Linear Predictor"),aes(ymin = 0, ymax = value),fill = "lightblue", alpha = 0.5,show.legend = FALSE) + # Fill under the line
#   geom_hline(yintercept=0,linetype="dashed") +
#   scale_color_manual(values = c("Linear Predictor" = "blue", "SVC" = "green", "SIC" = "red")) +
#   labs(x="Distance",
#        y="Log-intensity of type j",
#        color="Curve") +
#   guides(fill="none")
# fsave("svc_sic_explainer.png")


### Image of densities
sim_idx <- 125
num_pts <- 20
images_per_pt <- 1
# parameters to adjust
grid <- expand.grid(ratio=c(0.5,1,2,5,10),
                    num_points_per_type=c(20,80,150,300,500),
                    sim=1:5)
ratio <- grid$ratio[sim_idx]
np <- 150
n_dummy <- floor(np * ratio)
sim <- grid$sim[sim_idx]

# other parameters (likely don't need to adjust for now)
num_types <- 3
num_combos <- num_types - 1
num_points_per_type <- rep(np,num_types)

num_points_gen <- mean(num_points_per_type)
num_pt_groups <- 1
size_im <- 1500
potentials <- make_rbfs(n_basis_functions = 3, max_dist = 75, basis_function_sigma = 15)
mean_alpha <- log(np/size_im^2)
sigma_beta_global <- 0.5
sigma_alpha_global <- 0.5
sigma_beta_indiv <- 0.1
sigma_alpha_indiv <- 0.5
sigma_beta_local <- 0.1
sigma_alpha_local <- 0.5
scale_sigmas <- 5
grainsize <- 1
intensity_resolution <- 128


# derived parameters
seed <- as.integer(2024 + sim_idx)
coef_seed <- as.integer(2024 + sim)
print(seed)
print(coef_seed)
num_images <- images_per_pt * num_pts
sample_to_indiv <- rep(1:num_pts,each = images_per_pt)
num_pts_per_group <- ceiling(num_pts/num_pt_groups)
indiv_to_group <- rep(1:num_pt_groups,each = num_pts_per_group,length.out=num_pts)
num_pot <- length(potentials)
scale_sigma_betas <- seq(5,1,length.out=num_pot)
# file_data_stan <- paste0(path,"data_stan_sim",sim,"_ratio_",ratio,"_np_",np,".json")
# file_ground_truth <- paste0(path,"ground_truth_sim",sim,"_ratio_",ratio,"_np_",np,".rds")
# file_pats <- paste0(path,"pats_sim",sim,"_ratio_",ratio,"_np_",np,".rds")


# make logistic parameters
params <- make_simulation_parameters(mean_alpha,
                                  sigma_beta_global,
                                  # sigma_alpha_global,
                                  sigma_beta_indiv,
                                  # sigma_alpha_indiv,
                                  sigma_beta_local,
                                  # sigma_alpha_local,
                                  scale_sigmas,
                                  num_pt_groups,
                                  num_types,
                                  num_combos,
                                  num_pot,
                                  indiv_to_group,
                                  num_pts,
                                  num_images,
                                  sample_to_indiv,
                                  coef_seed)
betas_local <- params$betas_local[-1,]

set.seed(seed)
W <- owin(c(0,size_im),c(0,size_im))
area <- size_im^2
num_images <- num_pts * images_per_pt

i <- 1

if(num_types == 2) {
  pat <- rpoispp(lambda=num_points_gen/area,win = W)
  marks(pat) <- factor("t1")
} else {
  pat <- rmpoispp(lambda=rep(num_points_gen/area,num_types-1),win = W)
}

dens <- lapply(1:(num_types-1),\(j) {
  start <- (j-1)*num_pot+1
  stop <- j*num_pot
  coeffs <- betas_local[start:stop,i]
  custom_kernel <- Vectorize(function(x, y) {
    d <- sqrt(x^2 + y^2)  # Compute distance
    
    # Compute weighted sum of basis functions
    kernel_value <- sum(sapply(seq_along(coeffs), function(i) {
      coeffs[i] * potentials[[i]](d)
    }))
    
    return(kernel_value)
  })
  
  subs <- unmark(subset(pat,marks == j))
  dens <- smooth_density_fft(subs, custom_kernel, resolution = 256)
  # plot(dens)
  # points(subs$x,subs$y,pch=".",cex=3,col="blue")
})

ddf <- lapply(1:length(dens),\(j) as.data.frame(dens[[j]]) %>% mutate(type=j)) %>% bind_rows()
dens <- Reduce("+",dens)

ddf <- dens %>%
  as.data.frame() %>%
  mutate(type=num_types) %>%
  bind_rows(ddf)

subs <- lapply(1:(num_types-1),\(j) {
  sub <- unmark(subset(pat,marks == j))
  as.data.frame(sub) %>%
    mutate(type=j)
}) %>%
  bind_rows() %>%
  mutate(base=TRUE)

lambda_integral <- sum(exp(dens$v)) * (dens$xstep * dens$ystep)  # Approximate integral using grid summation

# Compute the necessary beta0[i] to achieve expected np points
beta0 <- log(150 / lambda_integral)
pat2 <- rpoispp(lambda = exp(dens + beta0))

pat2 %>%
  as.data.frame() %>%
  mutate(type=num_types,base=FALSE) %>%
  bind_rows(subs) %>%
  mutate(type=factor(type)) %>%
  mutate(type=fct_recode(type,
                         "Source 1"="1",
                         "Source 2"="2",
                         "Target"="3")) -> subs
 
ddf %>%
  as_tibble() %>%
  mutate(type=factor(type)) %>%
  mutate(type=fct_recode(type,
                         "Source 1"="1",
                         "Source 2"="2",
                         "Target"="3")) %>%
  ggplot(aes(x,y)) +
  geom_raster(aes(fill=value)) +
  geom_point(data=subs,aes(shape=base),color="green",size=1) +
  facet_wrap(~type,ncol=3) +
  scale_fill_viridis_c(option="magma") +
  coord_fixed(ratio = 1) +  # Ensures square panels
  # theme_classic() +
  labs(x="X",y="Y") +
  guides(shape="none",fill="none")
fsave("example_intensity",height=3)

x_seq <- seq(0,100,0.1)
weights <- c(0.8,-0.4,0.3)
rbfs <- make_rbfs(n_basis_functions = 3, max_dist = 75, basis_function_sigma = 15)

lapply(1:(num_types-1),\(j) {
  start <- (j-1)*num_pot+1
  stop <- j*num_pot
  weights <- betas_local[start:stop,i]
  pot_mat <- lapply(1:length(rbfs),\(i) rbfs[[i]](x_seq)*weights[i]) %>% do.call(cbind,.) %>% rowSums()
  tibble(x=x_seq,y=pot_mat,type=j)
  
}) %>%
  bind_rows() %>%
  mutate(type=factor(type)) %>%
  ggplot(aes(x,y)) +
  geom_line(aes(color=type,group=type),linewidth=1.5) +
  geom_hline(yintercept=0,linetype="dashed") +
  labs(x="Distance (microns)",y="Log-intensity",color="Type")
fsave("example_SICs")

# asummetric demo
set.seed(123)
# Adjustable parameters
n_tumor <- 20
n_clusters <- 5
cd8_per_tumor <- 5  # <-- adjust this to control how many CD8+ T cells cluster around each tumor

# Tumor cells: random spatial scatter
tumor <- data.frame(
  x = runif(n_tumor, 0.2, 0.8),
  y = runif(n_tumor, 0.2, 0.8),
  type = "Tumor"
)

# CD8+ T cells clustered around selected tumor cells
cd8 <- tumor %>%
  sample_n(n_clusters) %>%
  rowwise() %>%
  do({
    center_x <- .$x
    center_y <- .$y
    data.frame(
      x = rnorm(cd8_per_tumor, mean = center_x, sd = 0.03),
      y = rnorm(cd8_per_tumor, mean = center_y, sd = 0.03),
      type = "CD8+ T"
    )
  }) %>% bind_rows()

cells <- bind_rows(tumor, cd8)

# Plot
ggplot(cells, aes(x = x, y = y, color = type, shape = type)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Tumor" = "firebrick", "CD8+ T" = "forestgreen")) +
  scale_shape_manual(values = c("Tumor" = 17, "CD8+ T" = 16)) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "top") +
  coord_fixed() +
  labs(
    # title = "Asymmetric Spatial Interaction",
    # subtitle = "CD8+ T cells cluster near tumors, but not vice versa",
    x = NULL, y = NULL
  )
fsave("example_asymmetry",width=5)


# kcross demo

library(patchwork)

set.seed(123)

# 1. Generate asymmetric pattern: B clustered around A
win <- owin(c(0,1), c(0,1))

# Type A: base points (e.g., Tumor)
n_A <- 10
points_A <- runifpoint(n_A, win = win)

# Type B: clustered around A
n_B_per_A <- 5
sigma <- 0.03
points_B <- do.call(rbind, lapply(1:n_A, function(i) {
  center <- c(points_A$x[i], points_A$y[i])
  x <- rnorm(n_B_per_A, mean = center[1], sd = sigma)
  y <- rnorm(n_B_per_A, mean = center[2], sd = sigma)
  cbind(x, y)
}))
points_B <- points_B[points_B[,1] > 0 & points_B[,1] < 1 & points_B[,2] > 0 & points_B[,2] < 1,]

# Combine into multi-type ppp
marks_A <- rep("A", n_A)
marks_B <- rep("B", nrow(points_B))
multi_pp <- ppp(
  x = c(points_A$x, points_B[,1]),
  y = c(points_A$y, points_B[,2]),
  marks = factor(c(marks_A, marks_B)),
  window = win
)

# 2. Compute Kcross from A to B
K_AB <- Lcross(multi_pp, i = "A", j = "B")
plot(K_AB)

# 3. Point pattern plot
pp_df <- as.data.frame(multi_pp)

p1 <- ggplot(pp_df, aes(x = x, y = y, color = marks, shape = marks)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("A" = "firebrick", "B" = "forestgreen")) +
  scale_shape_manual(values = c("A" = 17, "B" = 16)) +
  theme_minimal(base_size = 14) +
  coord_fixed() +
  theme(legend.position = "top") +
  labs(
    # title = "Asymmetric Point Pattern",
    # subtitle = "Type B clustered around Type A",
    color = "Cell Type", shape = "Cell Type"
  )

# 4. Kcross plot
k_df <- data.frame(r = K_AB$r, K = K_AB$border)

p2 <- ggplot(k_df, aes(x = r, y = K)) +
  geom_line(color = "black", size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal(base_size = 14) +
  labs(
    # title = "Lcross: A → B",
    x = "Distance r",
    y = expression(L[AB](r))
  )
p2
# 5. Combine plots
combined_plot <- p1 + p2
print(combined_plot)
fsave("lcross_demo")



# newer cross demo

# Load required libraries
library(spatstat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Set seed for reproducibility
set.seed(123)

# Create study window
win <- square(1)

# NEW APPROACH: Create asymmetry through inhomogeneous intensity
# B points create "attraction fields" but are themselves random

# Step 1: Generate type B points (completely random)
n_B <- 50
B_points <- runifpoint(n_B, win)

# Step 2: Create an intensity surface for A points based on B locations
# A points will have higher probability near B points

# Create a fine grid for intensity calculation
grid_res <- 100
x_seq <- seq(0, 1, length.out = grid_res)
y_seq <- seq(0, 1, length.out = grid_res)
grid <- expand.grid(x = x_seq, y = y_seq)

# Calculate intensity at each grid point based on distance to nearest B point
min_intensity <- 10  # Background intensity
max_intensity <- 100  # Peak intensity near B points
attraction_range <- 0.15  # How far the attraction extends

grid$intensity <- min_intensity
for(i in 1:nrow(grid)) {
  # Find distance to nearest B point
  distances <- sqrt((grid$x[i] - B_points$x)^2 + (grid$y[i] - B_points$y)^2)
  min_dist <- min(distances)
  
  # Create intensity that decreases with distance from B points
  if(min_dist <= attraction_range) {
    attraction_strength <- exp(-min_dist / (attraction_range/3))
    grid$intensity[i] <- min_intensity + (max_intensity - min_intensity) * attraction_strength
  }
}

# Convert to intensity image
intensity_im <- im(matrix(grid$intensity, nrow = grid_res), 
                   xcol = x_seq, yrow = y_seq)

# Step 3: Generate A points using the inhomogeneous intensity
# This creates the asymmetry: A is attracted to B, but B doesn't "know" about A
A_points <- rpoispp(intensity_im)
n_A <- A_points$n

# Step 4: Combine into marked point pattern
all_x <- c(A_points$x, B_points$x)
all_y <- c(A_points$y, B_points$y)
marks <- factor(c(rep("A", n_A), rep("B", n_B)))

# Create marked point pattern
pp <- ppp(all_x, all_y, window = win, marks = marks)

# Step 5: Calculate cross-type summary statistics
# Gcross: nearest neighbor distance distributions
G_A_to_B <- Gcross(pp, i = "A", j = "B", correction = "best")
G_B_to_A <- Gcross(pp, i = "B", j = "A", correction = "best")

# Also calculate Kcross for comparison
K_A_to_B <- Kcross(pp, i = "A", j = "B", correction = "best")
K_B_to_A <- Kcross(pp, i = "B", j = "A", correction = "best")

# Step 6: Create plots

# Plot 1: Point pattern with intensity surface
plot_data <- data.frame(
  x = all_x,
  y = all_y,
  type = marks
)

# Add intensity surface as background
intensity_df <- data.frame(
  x = rep(x_seq, each = grid_res),
  y = rep(y_seq, grid_res),
  intensity = as.vector(intensity_im$v)
)

p1 <- ggplot() +
  geom_raster(data = intensity_df, aes(x = x, y = y, fill = intensity), alpha = 0.3) +
  scale_fill_gradient(low = "white", high = "yellow", name = "A intensity") +
  geom_point(data = plot_data, aes(x = x, y = y, color = type, shape = type), 
             size = 2.5, alpha = 0.9) +
  scale_color_manual(values = c("A" = "#E31A1C", "B" = "#1F78B4")) +
  scale_shape_manual(values = c("A" = 16, "B" = 17)) +
  labs(title = "Asymmetric Spatial Pattern",
       subtitle = "A (red) attracted to B (blue), B random",
       x = "X coordinate", y = "Y coordinate") +
  theme_minimal() +
  theme(legend.position = "bottom",
        aspect.ratio = 1) +
  coord_fixed()

# Plot 2: Gcross A -> B (should show attraction)
gcross_A_to_B_df <- data.frame(
  r = G_A_to_B$r,
  G = G_A_to_B$km,
  G_theo = G_A_to_B$theo
)

p2 <- ggplot(gcross_A_to_B_df, aes(x = r)) +
  geom_line(aes(y = G_theo), color = "gray50", linetype = "dashed", size = 1) +
  geom_line(aes(y = G), color = "#E31A1C", size = 1.2) +
  labs(title = "G-cross: A → B",
       subtitle = "Should show attraction\n(above random line)",
       x = "Distance r", 
       y = "G(r)") +
  theme_minimal() +
  ylim(0, 1)

# Plot 3: Gcross B -> A (should be closer to random)
gcross_B_to_A_df <- data.frame(
  r = G_B_to_A$r,
  G = G_B_to_A$km,
  G_theo = G_B_to_A$theo
)

p3 <- ggplot(gcross_B_to_A_df, aes(x = r)) +
  geom_line(aes(y = G_theo), color = "gray50", linetype = "dashed", size = 1) +
  geom_line(aes(y = G), color = "#1F78B4", size = 1.2) +
  labs(title = "G-cross: B → A",
       subtitle = "Should be closer to random\n(near gray line)",
       x = "Distance r", 
       y = "G(r)") +
  theme_minimal() +
  ylim(0, 1)

# Combine plots side by side
combined_plot <- p1 | p2 | p3

# Display the combined plot
print(combined_plot)

# Additional analysis: Calculate and compare summary statistics
cat("\n=== ASYMMETRY ANALYSIS ===\n")
cat("Generation method: Inhomogeneous Poisson process\n")
cat("- A points have higher intensity near B points\n")
cat("- B points are completely random\n\n")

# Calculate mean nearest neighbor distances
subset_A <- subset(pp, marks == "A")
subset_B <- subset(pp, marks == "B")

if(n_A > 0 && n_B > 0) {
  mean_dist_A_to_B <- mean(nncross(subset_A, subset_B)$dist)
  mean_dist_B_to_A <- mean(nncross(subset_B, subset_A)$dist)
  
  cat("Mean nearest neighbor distances:\n")
  cat(sprintf("A to nearest B: %.4f\n", mean_dist_A_to_B))
  cat(sprintf("B to nearest A: %.4f\n", mean_dist_B_to_A))
  cat(sprintf("Ratio (A->B)/(B->A): %.2f\n", mean_dist_A_to_B / mean_dist_B_to_A))
  
  # Calculate G-function values at a specific distance to quantify asymmetry
  test_distance <- 0.1
  idx <- which.min(abs(G_A_to_B$r - test_distance))
  
  G_AB_at_test <- G_A_to_B$iso[idx]
  G_BA_at_test <- G_B_to_A$iso[idx]
  G_theo_at_test <- G_A_to_B$theo[idx]
  
  cat(sprintf("\nG-function values at r = %.2f:\n", test_distance))
  cat(sprintf("G_A->B: %.3f\n", G_AB_at_test))
  cat(sprintf("G_B->A: %.3f\n", G_BA_at_test))
  cat(sprintf("Random: %.3f\n", G_theo_at_test))
  cat(sprintf("Asymmetry ratio: %.2f\n", G_AB_at_test / G_BA_at_test))
}

# Optional: Create Kcross comparison
kcross_A_to_B_df <- data.frame(
  r = K_A_to_B$r,
  K = K_A_to_B$iso,
  K_theo = K_A_to_B$theo
)

kcross_B_to_A_df <- data.frame(
  r = K_B_to_A$r,
  K = K_B_to_A$iso,
  K_theo = K_B_to_A$theo
)

p4 <- ggplot(kcross_A_to_B_df, aes(x = r)) +
  geom_line(aes(y = K_theo), color = "gray50", linetype = "dashed", size = 1) +
  geom_line(aes(y = K), color = "#E31A1C", size = 1.2) +
  labs(title = "K-cross: A → B",
       subtitle = "Should show aggregation",
       x = "Distance r", 
       y = "K(r)") +
  theme_minimal()

p5 <- ggplot(kcross_B_to_A_df, aes(x = r)) +
  geom_line(aes(y = K_theo), color = "gray50", linetype = "dashed", size = 1) +
  geom_line(aes(y = K), color = "#1F78B4", size = 1.2) +
  labs(title = "K-cross: B → A",
       subtitle = "Should be closer to random",
       x = "Distance r", 
       y = "K(r)") +
  theme_minimal()

# Show Kcross plots too
kcross_combined <- p1 | p4 | p5
cat("\n\nK-cross functions:\n")
print(kcross_combined)

