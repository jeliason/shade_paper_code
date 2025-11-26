library(tidyverse)
library(latex2exp)
library(spatstat)
library(SHADE)

# Load utility functions and constants
source("utils.R")

theme_set(theme_bw(base_size=14, base_family='Helvetica')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

fsave <- \(fname,height=5,width=8) {
  ggsave(paste0(figures_folder,fname,".pdf"),device=cairo_pdf, height=height, width=width, units="in")
}
figures_folder <- "manuscript/images/demo_plots/"


# Figure 2
x_seq <- seq(0,100,0.1)
weights <- c(0.8,-0.4,0.3)
rbfs <- make_rbfs(n_basis_functions = 3, max_dist = 75, basis_function_sigma = 15)
pot_mat <- lapply(1:length(rbfs),\(i) rbfs[[i]](x_seq)*weights[i]) %>% do.call(cbind,.) %>% rowSums()
tibble(x=x_seq,y=pot_mat) %>%
  filter(x >= MIN_INTERACTION_RADIUS) %>%
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


# Supplementary Figure S2
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

# Supplementary Figure S1
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
  filter(x >= MIN_INTERACTION_RADIUS) %>%
  ggplot(aes(x,y)) +
  geom_line(aes(color=type,group=type),linewidth=1.5) +
  geom_hline(yintercept=0,linetype="dashed") +
  labs(x="Distance (microns)",y="Log-intensity",color="Type")
fsave("example_SICs")
