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

generated quantities {
  vector[I] log_lik = rep_vector(0, I);
  real sum_log_lik;

  for (i in 1:I){
    log_lik[i] += normal_lpdf(beta[Ori[i]] | 0, 1);
    log_lik[i] += normal_lpdf(log1p_T[i] | b + beta[Ori[i]], sigma_t);
    if (Is_ase_het[i] == 1){
      log_lik[i] += normal_lpdf(logit_pi_alt[i] | 0, sigma_a);
    }
  }

  sum_log_lik = sum(log_lik);
}
