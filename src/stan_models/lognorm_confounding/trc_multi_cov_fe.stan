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
