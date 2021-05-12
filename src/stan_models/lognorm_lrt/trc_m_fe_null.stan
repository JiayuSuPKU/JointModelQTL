functions {
   real cis_reg_effect(real G, real r){
     // G genotype
     // r cis-regulatory effect, also log read ratio
      if (G == 2)
        return r;
      else if (G == 1)
        return log(1 + exp(r)) - log(2);
      else
        return 0; 
   }
}

data {
  int<lower=0> I; // number of samples
  // int<lower=0> J; // number of genes
  // int<lower=0> K; // number of test snps
  // int<lower=0> L; // number of exonic snps 
  int<lower=0> N; // number of confounding factors
  vector[I] G; // genotype
  vector<lower=0>[I] log1p_T; // log1p of total read counts
  matrix[I,N] X; // confounding factors
}

parameters {
  real<lower=0> b; // baseline expression
  real r; // cis-regulated effect
  real<lower=0> sigma_t; // variance of total read counts (log-scale)
  vector[N] beta; // effect size library size 
}

model {
  vector[I] mu_t;

  mu_t = b + X * beta;

  for (i in 1:I){
    mu_t[i] += cis_reg_effect(G[i], r);
  }

  log1p_T ~ normal(mu_t, sigma_t);
}
