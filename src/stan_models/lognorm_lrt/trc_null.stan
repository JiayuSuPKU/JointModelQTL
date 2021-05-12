data {
  int<lower=0> I; // number of samples
  // int<lower=0> J; // number of genes
  // int<lower=0> K; // number of test snps
  // int<lower=0> L; // number of exonic snps 
  vector<lower=0>[I] log1p_T; // log1p of total read counts
}

parameters {
  real<lower=0> b; // baseline expression
  real<lower=0> sigma_t;
}

model {
  log1p_T ~ normal(b, sigma_t);
}
