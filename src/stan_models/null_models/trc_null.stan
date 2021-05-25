data {
  int<lower=0> I; // number of samples
  vector<lower=0>[I] log1p_T; // log1p of total read counts
}

parameters {
  real<lower=0> b; // baseline expression
  real<lower=0> sigma_t;
}

model {
  log1p_T ~ normal(b, sigma_t);
}

generated quantities {
  vector[I] log_lik = rep_vector(0, I);
  real sum_log_lik;

  for (i in 1:I){
    log_lik[i] += normal_lpdf(log1p_T[i] | b, sigma_t);
  }

  sum_log_lik = sum(log_lik);
}
