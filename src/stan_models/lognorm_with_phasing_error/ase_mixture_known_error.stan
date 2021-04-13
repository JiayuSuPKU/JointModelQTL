functions {
   real cis_reg_effect(real G, real r){
     // G genotype
     // r cis-regulatory effect, also log read ratio
      if (G == 2)
        return r;
      else if (G == 1)
        return log(1 + exp(r)) - log(2);
      else
        return 0; 
   }
}

data {
  int<lower=0> I; // number of samples
  // int<lower=0> J; // number of genes
  // int<lower=0> K; // number of test snps
  // int<lower=0> L; // number of exonic snps 
  vector[I] G; // genotype
  vector[I] P; // phasing
  vector[I] Is_ase_het; // whether the gene region is heterogeneous
  vector[I] logit_pi_alt; // logit of the proportion of alt reads
  vector<lower=0, upper=1>[I] P_error; // phasing error rate for each sample
}

parameters {
  real<lower=0> r; // cis-regulated effect
  real<lower=0> sigma_a;
}

model {
  vector[I] mu_a;
  mu_a = P * r;

  for (i in 1:I){
    if (Is_ase_het[i] == 1){
      target += log_mix(1 - P_error[i],
        normal_lpdf(logit_pi_alt[i] | mu_a[i], sigma_a),
        normal_lpdf(logit_pi_alt[i] | -mu_a[i], sigma_a));
    }
  }
}

