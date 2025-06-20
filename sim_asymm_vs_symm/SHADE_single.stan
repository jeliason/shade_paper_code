data {
  // Basics
  int<lower=1> num_types;               // number of cell types (or interaction targets)
  int<lower=1> num_pot;                 // number of potentials
  int<lower=0> n_cells;                 // number of cells
  int<lower=0> d_cells;                 // number of spatial features (input dimension)

  // Outcome
  array[n_cells] int<lower=0,upper=1> is_cell;
  vector[n_cells] oset;

  matrix[n_cells,d_cells] x_cells;
  int<lower=1> beta_start; // index of the first beta parameter

  // Hyperparameters
  real mean_alpha;
  vector<lower=0>[num_pot] scale_sigma_betas;
  real<lower=0> scale_sigma_alpha;
  // real<lower=0> scale_sigmas;
}

transformed data {

  int num_combos = num_types - 1;
  array[num_combos,num_pot] int beta_idx;
  for(i in 1:num_combos) {
    for(j in 1:num_pot) {
      beta_idx[i,j] = (i-1) * num_pot + j + beta_start;
    }
  }
}

parameters {
  vector[d_cells] beta_local;
  vector<lower=0>[num_pot] sigma_beta_local;  
  real<lower=0> tau_alpha_local;

}

model {
  // --- Sample priors ---
  tau_alpha_local ~ normal(scale_sigma_alpha, 10);
  beta_local[1:(beta_start)] ~ normal(mean_alpha, tau_alpha_local);
  // sigma_beta_local ~ normal(0,10);
  // beta_local[(num_types+1):num_elements(beta_local)] ~ normal(0,sigma_beta_local);

  for (j in 1:num_pot) {
    sigma_beta_local[j] ~ normal(scale_sigma_betas[j], scale_sigma_betas[j]);
    to_vector(beta_local[beta_idx[:,j]]) ~ normal(0, sigma_beta_local[j]);
  }

  vector[n_cells] Xb = x_cells * beta_local;

  target += bernoulli_logit_lupmf(is_cell | Xb + oset);
}
