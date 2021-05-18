data {
  int<lower=0> I; // number of samples
  int<lower=0> N_indiv; // number of individuals
  vector<lower=0>[I] log1p_T; // log1p of total read counts
  vector[I] Is_ase_het; // whether the gene region is heterogeneous
  vector[I] logit_pi_alt; // logit of the proportion of alt reads
  int Ori[I]; // sample origin
}

parameters {
  real<lower=0> b; // baseline expression
  real<lower=0> sigma_t; // variance of total read counts (log-scale)
  real<lower=0> sigma_a;
  vector[N_indiv] beta; // individual effect
}

model {
  vector[I] mu_t;

  for (i in 1:I){
    mu_t[i] = b + beta[Ori[i]];
    if (Is_ase_het[i] == 1)
      logit_pi_alt[i] ~ normal(0, sigma_a);
  }

  beta ~ normal(0, 1); // prior on beta
  
  log1p_T ~ normal(mu_t, sigma_t);
}
