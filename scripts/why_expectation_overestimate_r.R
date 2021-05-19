## The expectation model has a tendency of overestimating R
stan_models <- load_models()
ase_expect_known <- stan_model(file = "./src/deprecated_stan_models/ase_expect_known.stan")

# helper function to estimate and extract r
estimateR.p_error <- function(n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
                              phi = 3, theta = 30, baseline = 3, r = 1.5, p_error = 0.1) {
  Y <- simulateCisEffect.s2s(
    n_i = n_i, maf = maf, prob_ref = prob_ref, prob_as = prob_as,
    phi = phi, theta = theta, baseline = baseline, r = r, p_error = p_error
  )

  ase <- optimizing(stan_models$basic$ase, data = Y$data, algorithm = "LBFGS", as_vector = F)
  ase_mix_known <- optimizing(
    stan_models$phasing_error$ase_mixture_known,
    data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F
  )
  ase_expect_known <- optimizing(
    ase_expect_known,
    data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F
  )

  out <- data.frame(
    r_est = c(ase$par$r, ase_mix_known$par$r, ase_expect_known$par$r),
    model = c("baseline", "mixture", "expectation"),
    r = r,
    p_error = p_error
  )

  return(out)
}

r_est_p_error <- data.frame()

for (p_error in c(0, 0.05, 0.1, 0.2, 0.3, 0.4)) {
  for (r in seq(1, 2, 0.03)) {
    r_est_p_error <- rbind(
      r_est_p_error, estimateR.p_error(p_error = p_error, r = r)
    )
  }
}

r_est_p_error %>%
  mutate(p_error = factor(p_error)) %>%
  ggplot(aes(x = r, y = r_est, group = p_error, color = p_error)) +
  geom_point() +
  geom_ribbon(stat = "smooth", method = "lm", se = TRUE, alpha = 0.15, aes(color = NULL)) +
  geom_line(stat = "smooth", method = "lm", size = 1.5) +
  ylim(0, 3) +
  facet_wrap(~model) +
  labs(title = "n_i=1000, maf=0.1, r=1.5, phi=3, theta=30") +
  theme_bw() +
  scale_color_npg()


## Explanation: The mixture model works through re-scaling observations, while
## the expectation model works through re-scaling R

k_scale <- function(x, p){
  return(
    (1 - p - p*exp(-2*x))/(1 - p + p*exp(-2*x))
  )
}

expect <- function(r_prime, p){
  er = ((1 - p)*exp(r_prime) - p)/(1 - p - p*exp(r_prime))
  return(log(er))
}


df_p_error <- data.frame()
for (p in c(0, 0.05, 0.1, 0.2, 0.3, 0.4)){
  x = seq(-5, 5, length.out = 1000)
  k = k_scale(x = x, p = p)
  r_prime = seq(0, 2, length.out = 1000)
  r = expect(r_prime = r_prime, p = p)

  df_p_error <- rbind(df_p_error, data.frame(x, k, r_prime, r, p))
}

df_p_error %>% mutate(p = factor(p)) %>%
  ggplot(aes(x = x, y = k, group = p, color = p)) +
  geom_line() +
  labs(x = "x*mu/sigma^2", color = "p_error") +
  scale_color_npg() +
  theme_bw()

df_p_error %>% mutate(p = factor(p)) %>%
  ggplot(aes(x = r_prime, y = r, group = p, color = p)) +
  geom_line() +
  labs(x = "R\'", y = "R", color = "p_error") +
  scale_color_npg() +
  theme_bw()
