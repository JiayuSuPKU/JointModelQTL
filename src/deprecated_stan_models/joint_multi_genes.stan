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
  real G[I]; // genotype
  vector[J] P[I]; // phasing
  vector<lower=0>[J] log1p_T[I]; // log1p of total read counts
  vector[J] Is_ase_het[I]; // whether the gene region is heterogeneous
  vector[J] logit_pi_alt[I]; // logit of the proportion of alt reads
}

parameters {
  vector<lower=0>[J] b; // baseline expression
  vector[J] r; // cis-regulated effect
  vector<lower=0>[J] sigma_t;
  vector<lower=0>[J] sigma_a;
}

model {
  vector[J] mu_t[I];
  vector[J] mu_a[I];
  
  for (i in 1:I){
    mu_t[i] = b + cis_reg_effect(G[i], r);
    mu_a[i] = P[i] .* r;
  }

  for (j in 1:J){
    log1p_T[,j] ~ normal(mu_t[,j], sigma_t[j]);
    for (i in 1:I){
      if (Is_ase_het[i,j] == 1)
        logit_pi_alt[i,j] ~ normal(mu_a[i,j], sigma_a[j]);
    }
  }
}

