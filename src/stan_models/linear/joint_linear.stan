data {
  int<lower=0> I; // number of samples
  vector[I] G; // genotype
  vector[I] P; // phasing
  vector<lower=0>[I] log1p_T; // log1p of total read counts
  vector[I] Is_ase_het; // whether the gene region is heterogeneous
  vector[I] logit_pi_alt; // logit of the proportion of alt reads
}

parameters {
  real<lower=0> b; // baseline expression
  real r; // cis-regulated effect
  real<lower=0> sigma_t;
  real<lower=0> sigma_a;
}

model {
  vector[I] mu_t = b + 0.5 * G * r;
  vector[I] mu_a = P * r;

  for (i in 1:I){
    if (Is_ase_het[i] == 1)
      logit_pi_alt[i] ~ normal(mu_a[i], sigma_a);
  }
  log1p_T ~ normal(mu_t, sigma_t);
}

generated quantities {
  vector[I] mu_t = b + 0.5 * G * r;
  vector[I] mu_a = P * r;
  vector[I] log_lik = rep_vector(0, I);
  real sum_log_lik;

  for (i in 1:I){
    log_lik[i] += normal_lpdf(log1p_T[i] | mu_t[i], sigma_t);
    if (Is_ase_het[i] == 1){
      log_lik[i] += normal_lpdf(logit_pi_alt[i] | mu_a[i], sigma_a);
    }
  }

  sum_log_lik = sum(log_lik);
}