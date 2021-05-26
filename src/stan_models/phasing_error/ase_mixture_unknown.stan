data {
  int<lower=0> I; // number of samples
  vector[I] P; // phasing
  vector[I] Is_ase_het; // whether the gene region is heterogeneous
  vector[I] logit_pi_alt; // logit of the proportion of alt reads
}

parameters {
  real<lower=0> r; // cis-regulated effect
  real<lower=0> sigma_a;
  real<lower=0,upper=0.5> p_error; // unknown phasing error, will be the same for each sample
}

model {
  vector[I] mu_a;
  mu_a = P * r;

  for (i in 1:I){
    if (Is_ase_het[i] == 1){
      target += log_mix(1 - p_error,
        normal_lpdf(logit_pi_alt[i] | mu_a[i], sigma_a),
        normal_lpdf(logit_pi_alt[i] | -mu_a[i], sigma_a));
    }
  }
}

generated quantities {
  vector[I] mu_a = P * r;
  vector[I] log_lik = rep_vector(0, I);
  real sum_log_lik;

  
  for (i in 1:I){
    if (Is_ase_het[i] == 1){
      log_lik[i] += log_mix(1 - p_error,
        normal_lpdf(logit_pi_alt[i] | mu_a[i], sigma_a),
        normal_lpdf(logit_pi_alt[i] | -mu_a[i], sigma_a));
    }
  }

  sum_log_lik = sum(log_lik);
}
