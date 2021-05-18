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
  vector[I] G; // genotype
  vector<lower=0>[I] log1p_T; // log1p of total read counts
}

parameters {
  real<lower=0> b; // baseline expression
  real r; // cis-regulated effect
  real<lower=0> sigma_t;
}

model {
  vector[I] mu_t;
  for (i in 1:I){
    mu_t[i] = b + cis_reg_effect(G[i], r);
  }
  log1p_T ~ normal(mu_t, sigma_t);
}

