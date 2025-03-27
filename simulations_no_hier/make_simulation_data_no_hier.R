library(spatstat)
library(tidyverse)
library(Matrix)
devtools::load_all()
source("scripts/util_sims_no_alpha.R")

# set environment
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV == "laptop") {
  path <- "./no_hier/"
  num_pts_per_group <- 10
  images_per_pt <- 2
  sim_idx <- 3
} else {
  path <- "/nfs/turbo/umms-ukarvind/joelne/BISTRO/no_hier/"
  args <- commandArgs(trailingOnly=TRUE)
  sim_idx <- as.numeric(args[1])
  num_pts_per_group <- 20
  images_per_pt <- 4
  print(sim_idx)
}

# parameters to adjust
grid <- expand.grid(sim=1:100)
sim <- grid$sim[sim_idx]

# other parameters (likely don't need to adjust for now)
num_types <- 3
ratio <- 2
np <- 150

num_points_per_type <- rep(np,num_types)
n_dummy <- floor(np * ratio)
num_combos <- num_types - 1
num_points_per_type <- rep(np,num_types)

num_points_gen <- mean(num_points_per_type)
num_pt_groups <- 2
num_pts <- num_pts_per_group * num_pt_groups
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
seed <- as.integer(3024 + sim_idx)
coef_seed <- as.integer(3024 + sim)
print(seed)
print(coef_seed)
num_images <- images_per_pt * num_pts
sample_to_indiv <- rep(1:num_pts,each = images_per_pt)
num_pts_per_group <- ceiling(num_pts/num_pt_groups)
indiv_to_group <- rep(1:num_pt_groups,each = num_pts_per_group,length.out=num_pts)
num_pot <- length(potentials)
scale_sigma_betas <- seq(5,1,length.out=num_pot)
file_data_stan <- paste0(path,"data_stan_sim",sim,".json")
file_ground_truth <- paste0(path,"ground_truth_sim",sim,".rds")
file_pats <- paste0(path,"pats_sim",sim,".rds")


# make logistic parameters
params <- make_acyclic_parameters(mean_alpha,
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


pats <- lapply(1:num_images,\(i) {
  # predict the last type from all of the others
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
    dens <- smooth_density_fft(subs, custom_kernel, resolution = 128)
    # plot(dens)
    # points(subs$x,subs$y,pch=".",cex=3,col="blue")
  })
  
  # lapply(dens,plot)
  dens <- Reduce("+",dens)
  # plot(dens)
  # plot(dens + beta0[i])
  # plot(exp(dens + beta0[i]))
  lambda_integral <- sum(exp(dens$v)) * (dens$xstep * dens$ystep)  # Approximate integral using grid summation
  
  # Compute the necessary beta0[i] to achieve expected np points
  beta0 <- log(np / lambda_integral)
  pat2 <- rpoispp(lambda = exp(dens + beta0))
  marks(pat2) <- factor(num_types)
  print(i)
  
  pat <- superimpose(pat,pat2)
  list(pat=pat,beta0=beta0)
})

plot(pats[[8]]$pat)
table(pats[[7]]$pat$marks)

saveRDS(pats,file_pats)

beta0_local <- sapply(pats,\(o) o$beta0)
params$betas_local[1,] <- beta0_local

beta0_indiv <- sapply(1:num_pts,\(i) mean(beta0_local[which(sample_to_indiv == i)]))
beta0_global <- sapply(1:num_pt_groups,\(i) mean(beta0_indiv[which(indiv_to_group == i)]))

params$betas_indiv[1,] <- beta0_indiv
params$betas_global[1,] <- beta0_global

data_lists <- lapply(1:num_images,\(i) {
  pat <- pats[[i]]$pat
  Q <- make_quadrature(pat,n_dummy = n_dummy)
  
  data_list <- make_acyclic_data(Q,focal_cell=as.character(num_types),potentials = potentials, verbose = FALSE)
  
  print(i)
  data_list
})


is_cell_list <- lapply(data_lists,\(o) o$response)
oset_list <- lapply(data_lists,\(o) o$offset)
is_cell <- unlist(is_cell_list)
oset <- unlist(oset_list)

x_cells <- do.call(rbind,lapply(data_lists,\(o) o$data))
n_cells = nrow(x_cells)
d_cells = ncol(x_cells)
sample_id <- lapply(1:num_images,\(i) {
  rep(i,length(is_cell_list[[i]]))
}) %>% unlist()

y_start_stop <- lapply(1:num_images,\(i) {
  idx <- sort(which(sample_id == i))
  c(idx[1],idx[length(idx)])
}) %>% do.call(rbind,.)

x_cells <- as(x_cells,"RsparseMatrix")

u <- as.integer(x_cells@p + 1)
w <- x_cells@x
v <- x_cells@j + 1

data_start_stop <- lapply(1:num_images,\(i) {
  row_start <- y_start_stop[i,1]
  row_stop <- y_start_stop[i,2]
  
  data_start <- u[row_start]
  data_end <- u[row_stop + 1] - 1
  c(data_start,data_end)
}) %>% do.call(rbind,.)

data_stan <- list(
  num_indiv = num_pts,
  num_types = num_types,
  num_pot = num_pot,
  num_pt_groups = num_pt_groups,
  n_cells = n_cells,
  d_cells = d_cells,
  is_cell = is_cell,
  oset = -oset,
  n_samples = num_images,
  sample_to_indiv = sample_to_indiv,
  indiv_to_group = indiv_to_group,
  sample_id = sample_id,
  
  mean_alpha = mean_alpha,
  scale_sigmas = scale_sigmas,
  scale_sigma_betas = scale_sigma_betas,
  scale_sigma_alpha = scale_sigmas,
  
  grainsize = grainsize,
  
  y_start_stop = y_start_stop,
  data_start_stop = data_start_stop,
  n_nz = length(w),
  w = w,
  v = v,
  u = u
)

print("writing json...")
write_json_chunked(data_stan,file_data_stan,chunk_size = 1e6)
print("json written.")

saveRDS(params,file_ground_truth)
