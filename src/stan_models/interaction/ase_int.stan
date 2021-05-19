data {
  int<lower=0> I; // number of samples
  int<lower=0> N_cond; // number of conditions
  int Cond[I]; // sample condition
  vector[I] P; // phasing
  vector[I] Is_ase_het; // whether the gene region is heterogeneous
  vector[I] logit_pi_alt; // logit of the proportion of alt reads
}

parameters {
  vector[N_cond] R; // individual effect
  real<lower=0> sigma_a;
}

model {
  vector[I] mu_a;
  for (i in 1:I){
    mu_a[i] = P[i] * R[Cond[i]];
  }

  for (i in 1:I){
    if (Is_ase_het[i] == 1){
      target += normal_lpdf(logit_pi_alt[i] | mu_a[i], sigma_a);
    }
  }
}
