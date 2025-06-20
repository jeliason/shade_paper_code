data {
  int<lower=0> N;                    // number of observations
  int<lower=0> K;                    // number of predictors
  array[N] int<lower=0,upper=1> y;   // binary outcomes
  matrix[N, K] X;                    // predictor matrix
  vector[N] oset;                  // offset vector
  
  // New data for predictions
  int<lower=0> N_new;                // number of new observations
  matrix[N_new, K] X_new;            // new predictor matrix
  vector[N_new] offset_new;          // new offset vector
}

parameters {
  real alpha;                        // intercept
  vector[K] beta;                    // coefficients for predictors
}

model {
  // Priors
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  
  // Likelihood
  y ~ bernoulli_logit(alpha + X * beta + oset);
}

generated quantities {
  // Linear predictor for new data
  vector[N_new] linear_pred_new = alpha + X_new * beta + offset_new;
  
  // Predicted probabilities for new data
  vector[N_new] y_pred_new = inv_logit(linear_pred_new);
  
  // Posterior predictive samples for new data
  array[N_new] int y_rep_new;
  for (n in 1:N_new) {
    y_rep_new[n] = bernoulli_logit_rng(linear_pred_new[n]);
  }
}
