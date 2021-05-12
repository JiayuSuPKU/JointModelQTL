data {
  int<lower=0> I; // number of samples
  // int<lower=0> J; // number of genes
  // int<lower=0> K; // number of test snps
  // int<lower=0> L; // number of exonic snps 
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
