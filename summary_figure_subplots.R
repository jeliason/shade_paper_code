library(spatstat)
library(tidyverse)
library(patchwork)
source("utils.R")

fsave <- \(fname,height=5,width=5,...) {
  ggsave(paste0(figures_folder,fname), height=height, width=width, units="in",...)
}
figures_folder <- "figures/summary_figures/"

seed <- 2026
set.seed(seed)
W <- owin(c(0,1500),c(0,1500))
area <- 1500^2
num_points_gen <- 100
num_types <- 4
num_pot <- 3
betas <- matrix(c(-0.4,0.3,0.1,
           0.2,-0.3,0.2,
           0.3,0.3,-0.2),ncol=3)
potentials <- make_rbfs(n_basis_functions = 3, max_dist = 80, basis_function_sigma = 40)

if(num_types == 2) {
  pat <- rpoispp(lambda=num_points_gen/area,win = W)
  marks(pat) <- factor("t1")
} else {
  pat <- rmpoispp(lambda=rep(num_points_gen/area,num_types-1),win = W)
}

dens_list <- lapply(1:(num_types-1),\(j) {
  coeffs <- betas[j,]
  custom_kernel <- Vectorize(function(x, y) {
    d <- sqrt(x^2 + y^2)  # Compute distance
    
    # Compute weighted sum of basis functions
    kernel_value <- sum(sapply(seq_along(coeffs), function(i) {
      coeffs[i] * potentials[[i]](d)
    }))
    
    return(kernel_value)
  })
  
  subs <- unmark(subset(pat,marks == j))
  dens <- smooth_density_fft(subs, custom_kernel, resolution = 128)
})

dens <- Reduce("+",dens_list)
lambda_integral <- sum(exp(dens$v)) * (dens$xstep * dens$ystep)  # Approximate integral using grid summation

# Compute the necessary beta0[i] to achieve expected np points
beta0 <- log(num_points_gen / lambda_integral)
pat2 <- rpoispp(lambda = exp(dens + beta0))
marks(pat2) <- "4"
# plot(dens)
# points(pat2$x,pat2$y,pch=".",cex=3,col="blue")


pat_final <- superimpose(pat,pat2)

marks(pat_final) <- fct_recode(marks(pat_final),
                               "Source 1"="1",
                               "Source 2"="2",
                               "Source 3"="3",
                               "Target"="4")

plot(pat_final)


# Combine both into one data frame
cells <- as.data.frame(pat_final) %>%
  mutate(marks = factor(marks,levels=c("Source 1","Source 2","Source 3","Target")))

# Define custom colors and shapes
custom_colors <- c("Source 1" = "steelblue", 
                   "Source 2" = "forestgreen", 
                   "Source 3" = "goldenrod", 
                   "Target"   = "darkred")

custom_shapes <- c("Source 1" = 16,  # solid circle
                   "Source 2" = 17,  # solid triangle
                   "Source 3" = 15,  # solid square
                   "Target"   = 18)  # solid diamond

# Create the plot
ggplot(cells, aes(x = x, y = y, color = marks, shape = marks)) +
  geom_point(size = 1, alpha = 0.9) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = custom_shapes) +
  theme_void(base_size = 10) +
  coord_fixed() +
  labs(x = NULL, y = NULL, color = "", shape = "") +
  theme(
    # legend.position = "top",
    plot.background = element_rect(fill = "white", color = NA),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
fsave("cells.pdf",height = 1.75,width = 2.9, device = cairo_pdf)

###### SIFs
# Assume pat is a data.frame with x and y columns for cell locations
dens_df <- lapply(1:length(dens_list),\(i) as.data.frame(dens_list[[i]]) %>% mutate(type=names(custom_colors)[i])) %>% bind_rows()
# Convert to data frame and ensure marks is a factor
patfinal_df <- as.data.frame(pat_final) %>%
  mutate(type = factor(marks, levels = c("Source 1", "Source 2", "Source 3", "Target")))

# Create the plot
dens_df %>%
  filter(type != "Target") %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = value)) +
  scico::scale_fill_scico(palette = "imola") +
  geom_point(data = patfinal_df %>% filter(type != "Target"), aes(x = x, y = y, shape = marks), 
             size = 1, alpha = 0.9) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = custom_shapes) +
  coord_fixed() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void(base_size = 10) +
  facet_wrap(~type) +
  labs(fill = "Spatial Interaction Field", color = "", shape = "") +
  guides(color="none",shape="none") +
  theme(
    legend.position = "top",
    # strip.text = element_blank(),
    legend.text = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
fsave("sifs_ind.pdf",height = 1.38,width = 2.9,device=cairo_pdf)

x_seq <- seq(0,100,0.1)
rbfs <- potentials
lapply(1:(num_types-1),\(j) {
  weights <- betas[j,]
  pot_mat <- lapply(1:length(rbfs),\(i) rbfs[[i]](x_seq)*weights[i]) %>% do.call(cbind,.) %>% rowSums()
  tibble(x=x_seq,y=pot_mat,type=j)
  
}) %>%
  bind_rows() %>%
  mutate(type=factor(type)) %>%
  ggplot(aes(x,y)) +
  geom_line(aes(color=type,group=type),linewidth=1.5) +
  geom_hline(yintercept=0,linetype="dashed") +
  facet_wrap(~type) +
  theme_bw(base_size = 20) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(1, 0, 1, 0), "lines")  # top, right, bottom, left
  ) +
  guides(color="none") +
  labs(x="Distance",y="",color="Type")
fsave("sics_ind.pdf",height = 3,width = 6, device = cairo_pdf)


as.data.frame(dens) %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = value)) +
  scico::scale_fill_scico(palette = "imola") +
  geom_point(data = patfinal_df %>% filter(type == "Target"), aes(x = x, y = y, shape = marks), 
             size = 1, alpha = 0.9) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = custom_shapes) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void(base_size=10) +
  labs(fill = "Conditional Intensity", color = "", shape = "") +
  guides(color="none",shape="none") +
  theme(
    legend.position = "top",
    strip.text = element_blank(),
    legend.text = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
fsave("sif_dep.pdf",height = 2.41,width = 2.9,device=cairo_pdf)

#### Benefits of SHADE
x_seq <- seq(0,100,0.1)
rbfs <- potentials
cohort_weights <- betas[1,]
cohort_y <- lapply(1:length(rbfs),\(i) rbfs[[i]](x_seq)*cohort_weights[i]) %>% do.call(cbind,.) %>% rowSums()

# Generate patient-level weights by perturbing the cohort weights
set.seed(42)
num_patients <- 6
patient_weights <- replicate(num_patients, cohort_weights + rnorm(length(cohort_weights), sd = 0.05), simplify = FALSE)

# Compute SICs for each patient
patient_dfs <- lapply(1:num_patients, function(p_id) {
  weights <- patient_weights[[p_id]]
  y <- lapply(1:length(rbfs), function(i) rbfs[[i]](x_seq) * weights[i]) %>%
    do.call(cbind, .) %>%
    rowSums()
  tibble(x = x_seq, y = y, patient = paste0("Patient ", p_id))
}) %>%
  bind_rows()

# Cohort-level data
cohort_df <- tibble(x = x_seq, y = cohort_y)

# Plot
ggplot() +
  geom_line(data = patient_dfs, aes(x = x, y = y, group = patient),
            linetype = "dotted",color="#F8766D", linewidth = 0.5) +
  geom_line(data = cohort_df, aes(x = x, y = y),color="#F8766D", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # annotate("text", x = 65, y = 0.4, label = "Cohort-level SIC", hjust = 0, size = 4.5) +
  # annotate("text", x = 65, y = 0.15, label = "Patient-level SICs", hjust = 0, size = 4.5, color = "gray40") +
  labs(x = "Distance", y = "Log-intensity") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x=element_blank(),
        axis.text.y = element_blank())
fsave("multi_example.pdf",height = 1.6,width = 2.9,device=cairo_pdf)


# Simulate posterior samples of weights (e.g., from MCMC)
set.seed(123)
n_samples <- 100
cohort_weights <- betas[2,]
weight_samples <- replicate(n_samples, cohort_weights + rnorm(length(cohort_weights), sd = 0.05), simplify = FALSE)

# Compute SIC for each sample
sic_samples <- lapply(weight_samples, function(w) {
  y <- lapply(1:length(rbfs), function(i) rbfs[[i]](x_seq) * w[i]) %>%
    do.call(cbind, .) %>%
    rowSums()
  tibble(x = x_seq, y = y)
}) %>%
  bind_rows(.id = "sample")

# Summarize mean and 95% credible interval
summary_df <- sic_samples %>%
  group_by(x) %>%
  summarise(
    y_mean = mean(y),
    y_lower = quantile(y, 0.025),
    y_upper = quantile(y, 0.975)
  )

# Plot
ggplot(summary_df, aes(x = x)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), fill = "#7CAE00", alpha = 0.3) +
  geom_line(aes(y = y_mean), color = "#7CAE00", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Distance", y = "Log-intensity") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x=element_blank(),
        axis.text.y = element_blank())
fsave("uq.pdf",height = 1.6,width = 2.9,device=cairo_pdf)

# Simulated SICs for two cohorts
cohort_1_weights <- betas[1,]
cohort_1_y <- lapply(1:length(rbfs),\(i) rbfs[[i]](x_seq)*cohort_1_weights[i]) %>% do.call(cbind,.) %>% rowSums()

cohort_2_weights <- betas[3,]
cohort_2_y <- lapply(1:length(rbfs),\(i) rbfs[[i]](x_seq)*cohort_2_weights[i]) %>% do.call(cbind,.) %>% rowSums()

# Example credible intervals (just faked here for demonstration)
lower_1 <- cohort_1_y - 0.1
upper_1 <- cohort_1_y + 0.1
lower_2 <- cohort_2_y - 0.08
upper_2 <- cohort_2_y + 0.08

# Combined with uncertainty
sic_df_ci <- tibble(
  x = rep(x_seq, 2),
  y = c(cohort_1_y, cohort_2_y),
  y_lower = c(lower_1, lower_2),
  y_upper = c(upper_1, upper_2),
  cohort = rep(c("Cohort A", "Cohort B"), each = length(x_seq))
)

# Plot with ribbon
ggplot(sic_df_ci, aes(x = x, y = y, color = cohort, fill = cohort)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Cohort A" = "#F8766D", "Cohort B" = "#00BFC4")) +
  scale_fill_manual(values = c("Cohort A" = "#F8766D", "Cohort B" = "#00BFC4")) +
  labs(x = "Distance", y = "Log-intensity", color = "Cohort", fill = "Cohort") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top") +
  theme(axis.text.x=element_blank(),
        axis.text.y = element_blank())
fsave("group_comp.pdf",height = 1.6,width = 2.9,device=cairo_pdf)
