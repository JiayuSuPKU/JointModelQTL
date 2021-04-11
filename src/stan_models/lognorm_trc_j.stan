functions {
   vector cis_reg_effect(real G, vector r){
     // G genotype
     // r cis-regulatory effect, also log read ratio
      if (G == 2)
        return r;
      else if (G == 1)
        return log(1 + exp(r)) - log(2);
      else
        return rep_vector(0, dims(r)[1]); 
   }
}

data {
  int<lower=0> I; // number of samples
  int<lower=0> J; // number of genes
  // int<lower=0> K; // number of test snps
  // int<lower=0> L; // number of exonic snps 
  real G[I]; // genotype
  vector<lower=0>[J] log1p_T[I]; // log1p of total read counts
}

parameters {
  vector<lower=0>[J] b; // baseline expression
  vector[J] r; // cis-regulated effect
  vector<lower=0>[J] sigma_t;
}

model {
  vector[J] mu_t[I];
  for (i in 1:I){
    mu_t[i] = b + cis_reg_effect(G[i], r);
  }
  for (j in 1:J){
    log1p_T[,j] ~ normal(mu_t[,j], sigma_t[j]);
  }
  // log1p_T ~ multi_normal(mu_t, diag_matrix(sigma_t));
}

