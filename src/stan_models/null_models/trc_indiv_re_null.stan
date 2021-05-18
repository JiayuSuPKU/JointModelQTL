data {
  int<lower=0> I; // number of samples
  int<lower=0> N_indiv; // number of individuals
  vector<lower=0>[I] log1p_T; // log1p of total read counts
  int Ori[I]; // sample origin
}

parameters {
  real<lower=0> b; // baseline expression
  real<lower=0> sigma_t; // variance of total read counts (log-scale)

  vector[N_indiv] beta; // individual effect
}

model {
  vector[I] mu_t;

  for (i in 1:I){
    mu_t[i] = b + beta[Ori[i]];
  }

  beta ~ normal(0, 1); // prior on beta
  
  log1p_T ~ normal(mu_t, sigma_t);
}
