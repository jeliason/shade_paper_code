library(spatstat)
library(tidyverse)
library(Matrix)
library(SHADE)
# install.packages("../SHADE/",type="source",repos=NULL)
library(posterior)
library(cmdstanr)
library(ggdist)

theme_set(theme_bw(base_size=14, base_family='Helvetica')+
            theme(
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
            )
)

fsave <- \(fname,height=5,width=8) {
  ggsave(paste0(figures_folder,fname,".pdf"),device = cairo_pdf, height=height, width=width, units="in")
}
figures_folder <- "./figures/demo_plots/"

#' Make dispersions dataframe
#'
#' @param Q logistic regression quadrature scheme, as generated by make_quadrature
#' @param potentials List of potential functions of interaction between two points, takes as input the Euclidean distance between points.
#'
#' @return sparse matrix of dispersions
#' @export
make_presence_dispersions <- function(Q,
                                      potentials,
                                      verbose=FALSE) {
  
  if (!is.factor(spatstat.geom::marks(Q$data))) {
    stop("The marks on the points in Q must be factors.")
  }
  
  m <- spatstat.geom::marks(Q)
  types <- levels(m)
  
  m_data <- spatstat.geom::marks(Q$data)
  m_all <- spatstat.geom::marks(Q)
  types_grid <- expand.grid(t1=types,t2=types) %>%
    filter(t1 != t2)
  
  pd <- spatstat.geom::crossdist(spatstat.geom::union.quad(Q),Q$data)
  
  p <- progressr::progressor(length(potentials))
  data <- lapply(1:length(potentials),\(i) {
    potential <- potentials[[i]]
    nm_pot <- names(potentials)[i]
    pot <- matrix(potential(pd),nrow = nrow(pd))

    diag(pot) <- 0
    
    pot_ix_jx <- lapply(1:nrow(types_grid),\(j) {
      t1 <- types_grid$t1[j]
      t2 <- types_grid$t2[j]
      ix <- which(m_all == t1)
      jx <- which(m_data == t2)
      
      as.matrix(pot[ix,jx])
    })

    dispersions <- lapply(1:nrow(types_grid),\(j) {
      t1 <- types_grid$t1[j]
      ix <- which(m_all == t1)
      vec <- rep(0,nrow(pot))
      vec[ix] <- rowSums(pot_ix_jx[[j]])
      vec <- Matrix::Matrix(vec,sparse = TRUE)
    }) %>%
      do.call(cbind,.)
    
    colnames(dispersions) <- paste0(nm_pot,"_",types_grid$t1,"_",types_grid$t2)
    
    p()
    return(dispersions)
    
  }) %>%
    do.call(cbind,.)
  # Enforce symmetry: average (A->B) and (B->A), and keep only one
  keep_cols <- c()
  
  for(i in 1:num_pot) {
    nm_pot <- names(potentials)[i]
    for (j in 1:nrow(types_grid)) {
      t1 <- as.character(types_grid$t1[j])
      t2 <- as.character(types_grid$t2[j])
      
      if (t1 <= t2) {  # Canonical direction
        col1 <- paste0(nm_pot, "_", t1, "_", t2)
        col2 <- paste0(nm_pot, "_", t2, "_", t1)
        
        avg_col <- (data[, col1] + data[, col2]) / 2
        data[, col1] <- avg_col
        keep_cols <- c(keep_cols, col1)
      }
    }
  }
  # Subset to keep only canonical (unique) symmetric columns
  data <- data[, keep_cols, drop = FALSE]
  nms <- colnames(data)
  xy_part <- stringr::str_extract(nms, "_[A-Za-z0-9+ ]+_[A-Za-z0-9+ ]+")
  ord <- order(xy_part)
  data <- data[,ord]
  
  return(data)
}

# set environment
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV != "HPC") {
  path <- "./sim_dummy_points/data/"
  num_pts <- 5
  images_per_pt <- 1
} else {
  path <- "./sim_dummy_points/data/"
  args <- commandArgs(trailingOnly=TRUE)
  sim_idx <- as.numeric(args[1])
  num_pts <- 40
  images_per_pt <- 2
  print(sim_idx)
}

ratio <- 2
np <- 100
n_dummy <- floor(np * ratio)

# other parameters (likely don't need to adjust for now)
num_types <- 2
num_combos <- num_types - 1
num_points_per_type <- rep(np,num_types)

num_points_gen <- mean(num_points_per_type)
num_pt_groups <- 0
size_im <- 1500
potentials <- make_rbfs(n_basis_functions = 3, max_dist = 75, basis_function_sigma = 15)
mean_alpha <- log(np/size_im^2)
sigma_beta_local <- 1
sigma_alpha_local <- 0.5
scale_sigmas <- 5
grainsize <- 1
intensity_resolution <- 128


# derived parameters
seed <- as.integer(2028)
coef_seed <- as.integer(2028)
print(seed)
print(coef_seed)
num_images <- images_per_pt * num_pts
num_pot <- length(potentials)
scale_sigma_betas <- seq(5,1,length.out=num_pot)
types_grid <- expand.grid(t1=num_types,t2=1:(num_types-1))
num_combos <- nrow(types_grid)
combo_pot <- num_combos * num_pot

# betas_local <- matrix(rnorm(num_images*combo_pot,mean=0,sd=sigma_beta_local),ncol=num_images,nrow=combo_pot)
betas_local <- matrix(rep(c(0,2,1),times=num_images),ncol=num_images)
alphas_local <- rnorm(num_images*(num_types-1),mean=mean_alpha,sd=sigma_alpha_local)

# betas_local <- rbind(as.vector(alphas_local),betas_local)
betas_local
# make logistic parameters
# params <- make_simulation_parameters(mean_alpha,
#                                      sigma_beta_global,
#                                      sigma_beta_indiv,
#                                      sigma_beta_local,
#                                      scale_sigmas,
#                                      num_pt_groups,
#                                      num_types,
#                                      num_combos,
#                                      num_pot,
#                                      indiv_to_group,
#                                      num_pts,
#                                      num_images,
#                                      sample_to_indiv,
#                                      coef_seed)
# betas_local <- params$betas_local[-1,]

set.seed(seed)
W <- owin(c(0,size_im),c(0,size_im))
area <- size_im^2
num_images <- num_pts * images_per_pt


outs <- lapply(1:num_images,\(i) {
  # predict the last type from all of the others
  if(num_types == 2) {
    pat <- rpoispp(lambda=num_points_gen/area,win = W)
    marks(pat) <- factor("1")
  } else {
    pat <- rmpoispp(lambda=rep(num_points_gen/area,num_types-1),win = W)
  }
  
  dens <- lapply(1:(num_types-1),\(j) {
    start <- (j-1)*num_pot+1
    stop <- j*num_pot
    coeffs <- betas_local[start:stop,i]
    print(coeffs)
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
  
  dens <- Reduce("+",dens)
  lambda_integral <- sum(exp(dens$v)) * (dens$xstep * dens$ystep)  # Approximate integral using grid summation
  
  # Compute the necessary beta0[i] to achieve expected np points
  beta0 <- log(np / lambda_integral)
  pat2 <- rpoispp(lambda = exp(dens + beta0))
  marks(pat2) <- factor(num_types)
  print(i)
  
  pat <- superimpose(pat,pat2)
  list(pat=pat,beta0=beta0)
})

pats <- lapply(outs,\(o) o$pat)
pat <- pats[[3]]
plot(pat)

gx <- alltypes(pat,Gcross)
plot(gx)

jx <- alltypes(pat,Jcross)
plot(jx)

beta0_local <- sapply(outs,\(o) o$beta0)

if (is.null(scale_sigma_betas)) {
  scale_sigma_betas <- seq(5, 1, length.out = n_basis_functions)
}

Qs <- lapply(pats, make_quadrature, n_dummy = n_dummy)

type <- "1"

offset_asymm <- lapply(Qs, function(Q) {
  log(spatstat.geom::intensity(Q$dummy)) |>
    tibble::enframe() |>
    dplyr::right_join(tibble::tibble(name = spatstat.geom::marks(Q)), by = "name") |>
    dplyr::filter(name == type) |>
    dplyr::pull(value)
})

data_lists_asymm <- lapply(Qs, function(Q) make_data(Q, potentials, type, verbose = FALSE))

stan_datas_asymm <- lapply(1:length(data_lists_asymm),\(i) {
  d <- data_lists_asymm[[i]]
  stan_data <- list(
    num_types = num_types,
    num_pot = num_pot,
    n_cells = nrow(d$data),
    d_cells = ncol(d$data),
    is_cell = d$response,
    oset = -offset_asymm[[i]],
    x_cells = as.matrix(d$data),
    mean_alpha = mean_alpha,
    scale_sigma_betas = scale_sigma_betas,
    scale_sigma_alpha = scale_sigmas,
    scale_sigmas = scale_sigmas,
    beta_start = 1
  )
})

mod <- cmdstan_model("sim_asymm_vs_symm/SHADE_single.stan")
fits_asymm <- lapply(stan_datas_asymm,\(stan_data) {
  fit <- mod$variational(stan_data, 
                         draws = 1e3, 
                         refresh = 5)
})


# RMSE

estimate <- lapply(1:length(stan_datas_asymm), \(i) {
  fit <- fits_asymm[[i]]
  tryCatch({
    estimate <- as_draws_rvars(fit$draws(variables = "beta_local"))$beta_local
    estimate <- estimate[-1] # remove alphas
  },error=\(e) return(rep(NA,combo_pot)))
}) %>%
  do.call(cbind,.)

estimate

true <- betas_local

rmse <- sqrt(rvar_mean((estimate - true)^2,na.rm=TRUE))
rmse_close <- sqrt(rvar_mean(((estimate - true)^2)[seq(1,combo_pot,by=num_pot),],na.rm=TRUE))
rmse_mid <- sqrt(rvar_mean(((estimate - true)^2)[seq(2,combo_pot,by=num_pot),],na.rm=TRUE))
rmse_far <- sqrt(rvar_mean(((estimate - true)^2)[seq(3,combo_pot,by=num_pot),],na.rm=TRUE))

asymm_df <- tibble(rmse=c(rmse,rmse_close,rmse_mid,rmse_far),
             level="beta_local",
             scale=c("all","close","mid","far"))
  
asymm_df

## symmetric model

offset_symm <- lapply(Qs, function(Q) {
  log(spatstat.geom::intensity(Q$dummy)) |>
    tibble::enframe() |>
    dplyr::right_join(tibble::tibble(name = spatstat.geom::marks(Q)), by = "name") |>
    dplyr::pull(value)
})

data_lists_symm <- lapply(Qs,\(Q) {
  dispersions <- make_presence_dispersions(Q,potentials,verbose)
  response <- as.numeric(spatstat.geom::is.data(Q))
  
  beta0 = spatstat.geom::marks(Q)
  onehot <- Matrix::sparse.model.matrix(~0+as.factor(beta0))
  colnames(onehot) <- paste0("beta0_",levels(beta0))
  data <- cbind(onehot,dispersions)
  
  list(data=data,response=response)
})

stan_datas_symm <- lapply(1:length(data_lists_symm),\(i) {
  d <- data_lists_symm[[i]]
  stan_data <- list(
    num_types = num_types,
    num_pot = num_pot,
    n_cells = nrow(d$data),
    d_cells = ncol(d$data),
    is_cell = d$response,
    oset = -offset_symm[[i]],
    x_cells = as.matrix(d$data),
    mean_alpha = mean_alpha,
    scale_sigma_betas = scale_sigma_betas,
    scale_sigma_alpha = scale_sigmas,
    scale_sigmas = scale_sigmas,
    beta_start = 2
  )
})

mod <- cmdstan_model("sim_asymm_vs_symm/SHADE_single.stan")
fits_symm <- lapply(stan_datas_symm,\(stan_data) {
  fit <- mod$variational(stan_data, 
                         draws = 1e3, 
                         refresh = 5)
})

estimate <- lapply(1:length(stan_datas_symm), \(i) {
  fit <- fits_symm[[i]]
  estimate <- as_draws_rvars(fit$draws(variables = "beta_local"))$beta_local
  estimate <- estimate[-c(1,2)] # remove alphas
}) %>%
  do.call(cbind,.)

estimate <- estimate/2
estimate

rmse <- sqrt(rvar_mean((estimate - true)^2))
rmse_close <- sqrt(rvar_mean(((estimate - true)^2)[seq(1,combo_pot,by=num_pot),]))
rmse_mid <- sqrt(rvar_mean(((estimate - true)^2)[seq(2,combo_pot,by=num_pot),]))
rmse_far <- sqrt(rvar_mean(((estimate - true)^2)[seq(3,combo_pot,by=num_pot),]))

symm_df <- tibble(rmse=c(rmse,rmse_close,rmse_mid,rmse_far),
             level="beta_local",
             scale=c("all","close","mid","far"))

symm_df
asymm_df

symm_df %>%
  mutate(Model="Symmetric") %>%
  bind_rows(asymm_df %>% mutate(Model = "Asymmetric")) %>%
  mutate(scale = factor(scale,levels=c("all","close","mid","far"))) %>%
  ggplot(aes(scale,ydist=rmse,color=Model,fill=Model)) +
  stat_halfeye() +
  labs(y="RMSE")
fsave("symm_vs_asymm")
