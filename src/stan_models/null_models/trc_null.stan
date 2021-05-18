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
