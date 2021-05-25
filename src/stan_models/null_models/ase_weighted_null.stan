functions {
   real scale_sigma(real sigma_a, real A_ref, real A_alt, real ceiling){
      return sigma_a ./ fmin(A_ref + A_alt, ceiling);
   }
}

data {
  int<lower=0> I; // number of samples
  vector[I] Is_ase_het; // whether the gene region is heterogeneous
  vector[I] logit_pi_alt; // logit of the proportion of alt reads
  vector[I] A_ref; // ref allele read counts
  vector[I] A_alt; // alt allele read counts
}

parameters {
  real<lower=0> sigma_a;
}

model {
  for (i in 1:I){
    if (Is_ase_het[i] == 1){
      target += normal_lpdf(logit_pi_alt[i] | 0, scale_sigma(sigma_a, A_ref[i], A_alt[i], 100));
    }
  }
}

generated quantities {
  vector[I] log_lik = rep_vector(0, I);
  real sum_log_lik;

  for (i in 1:I){
    if (Is_ase_het[i] == 1){
      log_lik[i] += normal_lpdf(logit_pi_alt[i] | 0, scale_sigma(sigma_a, A_ref[i], A_alt[i], 100));
    }
  }

  sum_log_lik = sum(log_lik);
}
