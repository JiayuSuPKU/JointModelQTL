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
  vector[I] Lib_size; // relative library size
}

parameters {
  real<lower=0> b; // baseline expression
  real r; // cis-regulated effect
  real<lower=0> sigma_t; // variance of total read counts (log-scale)
  real beta; // effect size library size 
}

model {
  vector[I] mu_t;

  mu_t = b + Lib_size * beta;

  for (i in 1:I){
    mu_t[i] += cis_reg_effect(G[i], r);
  }

  log1p_T ~ normal(mu_t, sigma_t);
}
