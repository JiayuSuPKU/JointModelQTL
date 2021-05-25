data {
  int<lower=0> I; // number of samples
  vector[I] G; // genotype
  vector<lower=0>[I] log1p_T; // log1p of total read counts
}

parameters {
  real<lower=0> b; // baseline expression
  real r; // cis-regulated effect
  real<lower=0> sigma_t;
}

model {
  vector[I] mu_t = b + 0.5 * G * r;
  log1p_T ~ normal(mu_t, sigma_t);
}

generated quantities {
  vector[I] mu_t = b + 0.5 * G * r;
  vector[I] log_lik = rep_vector(0, I);
  real sum_log_lik;

  for (i in 1:I){
    log_lik[i] += normal_lpdf(log1p_T[i] | mu_t[i], sigma_t);
  }

  sum_log_lik = sum(log_lik);
}
