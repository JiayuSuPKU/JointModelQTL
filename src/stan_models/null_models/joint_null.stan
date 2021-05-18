data {
  int<lower=0> I; // number of samples
  vector<lower=0>[I] log1p_T; // log1p of total read counts
  vector[I] Is_ase_het; // whether the gene region is heterogeneous
  vector[I] logit_pi_alt; // logit of the proportion of alt reads
}

parameters {
  real<lower=0> b; // baseline expression
  real<lower=0> sigma_t;
  real<lower=0> sigma_a;
}

model {
  for (i in 1:I){
    if (Is_ase_het[i] == 1)
      logit_pi_alt[i] ~ normal(0, sigma_a);
  }
  log1p_T ~ normal(b, sigma_t);
}
