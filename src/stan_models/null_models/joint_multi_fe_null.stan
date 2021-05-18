data {
  int<lower=0> I; // number of samples
  int<lower=0> N; // number of confounding factors
  vector<lower=0>[I] log1p_T; // log1p of total read counts
  vector[I] Is_ase_het; // whether the gene region is heterogeneous
  vector[I] logit_pi_alt; // logit of the proportion of alt reads
  matrix[I,N] X; // confounding factors
}

parameters {
  real<lower=0> b; // baseline expression
  real<lower=0> sigma_t; // variance of total read counts (log-scale)
  real<lower=0> sigma_a;
  vector[N] beta; // effect size library size 
}

model {
  vector[I] mu_t;
  mu_t = b + X * beta;

  for (i in 1:I){
    if (Is_ase_het[i] == 1)
      logit_pi_alt[i] ~ normal(0, sigma_a);
  }
  
  log1p_T ~ normal(mu_t, sigma_t);
}
