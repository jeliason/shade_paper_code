functions {
	real partial_sum_lpmf(array[] int slice_samples,
                        int sample_start_no_use, int sample_end_no_use,
												array[] int is_cell,
												array[,] int y_start_stop,
												array[,] int data_start_stop,
												vector w, 
												array[] int v,
												array[] int u,
												int n_cols,
                        vector oset,
                        matrix beta_local) {
    real partial_sum = 0.0;


    int n_samples = num_elements(slice_samples);

    for (s_idx in 1:n_samples) {
			int s = slice_samples[s_idx];
			int y_start = y_start_stop[s,1];
      int y_stop = y_start_stop[s,2];

			int n = y_stop - y_start + 2;
  		array[n] int slice_row_ptr = u[y_start:(y_stop + 1)];
			for (i in 1:n) {
				slice_row_ptr[i] = slice_row_ptr[i] - u[y_start] + 1;
			}

			int data_start = data_start_stop[s,1];
      int data_stop = data_start_stop[s,2];

			int n_rows = n - 1;
      // matrix vector multiplication

      vector[n_rows] Xb = csr_matrix_times_vector(n_rows,n_cols,w[data_start:data_stop],v[data_start:data_stop],slice_row_ptr,beta_local[:,s]);

      partial_sum += bernoulli_logit_lupmf(is_cell[y_start:y_stop] | Xb + oset[y_start:y_stop]);
    }

    return partial_sum;
  }
}


data {
  int<lower=1> num_types;
  int<lower=1> num_pot;
  int<lower=0> n_cells;
  int<lower=0> d_cells;
  array[n_cells] int<lower=0,upper=1> is_cell;
  vector[n_cells] oset;

  int<lower=1> n_samples;
  array[n_samples,2] int<lower=1,upper=n_cells> y_start_stop;

	int n_nz;

	array[n_samples,2] int<lower=1,upper=n_nz> data_start_stop;

	vector[n_nz] w;
	array[n_nz] int<lower=1,upper=n_cells> v;
	array[(n_cells+1)] int<lower=1,upper=(n_nz+1)> u;


  array[n_cells] int<lower=1,upper=n_samples> sample_id;

	real mean_alpha;
  vector<lower=0>[num_pot] scale_sigma_betas;
	real<lower=0> scale_sigma_alpha;
	real<lower=0> scale_sigmas;
  int<lower=1> grainsize;

}

transformed data {

  int num_combos = num_types - 1;
  array[num_combos,num_pot] int beta_idx;
  int beta_start = 1;
  for(i in 1:num_combos) {
    for(j in 1:num_pot) {
      beta_idx[i,j] = (i-1) * num_pot + j + beta_start;
    }
  }

  array[n_samples] int samples_list;
  for(s in 1:n_samples) {
    samples_list[s] = s;
  }
}

parameters {
  matrix[d_cells,n_samples] beta_local;
  vector<lower=0>[num_pot] sigma_beta_local;
  real<lower=0> tau_alpha_local;
}

model {
	for(j in 1:num_pot) {
		sigma_beta_local[j] ~ normal(scale_sigma_betas[j],scale_sigma_betas[j]);

		to_vector(beta_local[beta_idx[:,j],:]) ~ normal(0,sigma_beta_local[j]);
	}
	tau_alpha_local ~ normal(scale_sigma_alpha,10);
	to_vector(beta_local[1,:]) ~ normal(mean_alpha,tau_alpha_local);

	target += reduce_sum(partial_sum_lpmf, samples_list,
                     grainsize,
										 is_cell,
										 y_start_stop,
										 data_start_stop,
										 w,
										 v,
										 u,
										 d_cells,
										 oset,
										 beta_local);

}
