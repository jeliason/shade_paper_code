data {
  int<lower=0> N;                    // number of observations
  int<lower=0> K;                    // number of predictors
  array[N] int<lower=0,upper=1> y;   // binary outcomes
  matrix[N, K] X;                    // predictor matrix
  vector[N] oset;                  // offset vector
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
