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
  vector<lower=0>[I] log1p_T; // log1p of total read counts
  matrix[I, I] Cov_ori; // sample origin index matrix, Cov_ori[i, j] := 1 if i, j are from the same individual
}

parameters {
  real<lower=0> b; // baseline expression
  real r; // cis-regulated effect
  real<lower=0> sigma_sample; // per-sample variance 
  real<lower=0> sigma_indiv; // per-individual variance
}

model {
  vector[I] mu_t; // mean
  matrix[I, I] sigma; // covariance matrix
  matrix[I, I] L;  // cholesky factor


  for (i in 1:I){
    mu_t[i] = b + cis_reg_effect(G[i], r);
  }

  sigma = diag_matrix(rep_vector(sigma_sample, I));
  sigma = sigma + sigma_indiv * Cov_ori; 
  L=cholesky_decompose(sigma);

  log1p_T ~ multi_normal_cholesky(mu_t, L);
}
