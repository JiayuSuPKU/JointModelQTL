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
  int<lower=0> N_cond; // number of conditions
  int Cond[I]; // sample condition
  vector[I] G; // genotype
  vector<lower=0>[I] log1p_T; // log1p of total read counts
}

parameters {
  real<lower=0> b; // baseline expression
  vector[N_cond] R; // individual effect
  real<lower=0> sigma_t;
}

model {
  vector[I] mu_t;
  for (i in 1:I){
    mu_t[i] = b + cis_reg_effect(G[i], R[Cond[i]]);
  }
  log1p_T ~ normal(mu_t, sigma_t);
}

generated quantities {
  vector[I] mu_t;
  vector[I] log_lik = rep_vector(0, I);
  real sum_log_lik;

  for (i in 1:I){
    mu_t[i] = b + cis_reg_effect(G[i], R[Cond[i]]);
    log_lik[i] += normal_lpdf(log1p_T[i] | mu_t[i], sigma_t);
  }

  sum_log_lik = sum(log_lik);
}

