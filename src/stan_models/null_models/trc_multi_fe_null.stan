data {
  int<lower=0> I; // number of samples
  int<lower=0> N; // number of confounding factors
  vector<lower=0>[I] log1p_T; // log1p of total read counts
  matrix[I,N] X; // confounding factors
}

parameters {
  real<lower=0> b; // baseline expression
  real<lower=0> sigma_t; // variance of total read counts (log-scale)
  vector[N] beta; // effect size library size 
}

model {
  vector[I] mu_t;
  mu_t = b + X * beta;

  log1p_T ~ normal(mu_t, sigma_t);
}

generated quantities {
  vector[I] mu_t = b + X * beta;
  vector[I] log_lik = rep_vector(0, I);
  real sum_log_lik;

  for (i in 1:I){
    log_lik[i] += normal_lpdf(log1p_T[i] | mu_t[i], sigma_t);
  }

  sum_log_lik = sum(log_lik);
}
