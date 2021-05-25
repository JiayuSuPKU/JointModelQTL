data {
  int<lower=0> I; // number of samples
  vector[I] Is_ase_het; // whether the gene region is heterogeneous
  vector[I] logit_pi_alt; // logit of the proportion of alt reads
}

parameters {
  real<lower=0> sigma_a;
}

model {
  for (i in 1:I){
    if (Is_ase_het[i] == 1){
      target += normal_lpdf(logit_pi_alt[i] | 0, sigma_a);
    }
  }
}

generated quantities {
  vector[I] log_lik = rep_vector(0, I);
  real sum_log_lik;

  for (i in 1:I){
    if (Is_ase_het[i] == 1){
      log_lik[i] += normal_lpdf(logit_pi_alt[i] | 0, sigma_a);
    }
  }

  sum_log_lik = sum(log_lik);
}
