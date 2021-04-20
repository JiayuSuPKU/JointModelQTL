library(tidyverse)
library(ggsci)
library(ggrepel)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## L-BFGS
Y <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3, r = 1.5, p_error = 0
)

# config stan models
lognorm_ase <- stan_model(file = "./src/stan_models/lognorm_ase.stan")
lognorm_ase_mix_known <- stan_model(
  file = "./src/stan_models/lognorm_with_phasing_error/ase_mixture_known_error.stan"
)
lognorm_ase_mix_unknown <- stan_model(
  file = "./src/stan_models/lognorm_with_phasing_error/ase_mixture_unknown_error.stan"
)
lognorm_ase_expect_known <- stan_model(
  file = "./src/stan_models/lognorm_with_phasing_error/ase_expect_known_error.stan"
)


lbfgs_ase <- optimizing(
  lognorm_ase, data = Y$data, algorithm = "LBFGS", as_vector = F,
  draw = 10, importance_resampling = T)
lbfgs_ase_mix_known <- optimizing(
  lognorm_ase_mix_known,
  data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F
)
lbfgs_ase_mix_unknown <- optimizing(
  lognorm_ase_mix_unknown,
  data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F
)
lbfgs_ase_expect_known <- optimizing(
  lognorm_ase_expect_known,
  data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F
)

# runtime comparisons
system.time({
  optimizing(
    lognorm_ase_mix_known,
    data = Y$data, algorithm = "LBFGS"
  )
})
system.time({
  stan(
    file = "./src/stan_models/lognorm_with_phasing_error/ase_mixture_unknown_error.stan",
    data = Y$data
  )
})

# estimate r
estimateR.p_error <- function(n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
                              phi = 3, theta = 30, baseline = 3, r = 1.5, p_error = 0.1) {
  Y <- simulateCisEffect.s2s(
    n_i = n_i, maf = maf, prob_ref = prob_ref, prob_as = prob_as,
    phi = phi, theta = theta, baseline = baseline, r = r, p_error = p_error
  )

  ase <- optimizing(lognorm_ase, data = Y$data, algorithm = "LBFGS", as_vector = F)
  ase_mix_known <- optimizing(
    lognorm_ase_mix_known,
    data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F
  )
  ase_expect_known <- optimizing(
    lognorm_ase_expect_known,
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
  geom_ribbon(stat='smooth', method = "lm", se=TRUE, alpha=0.15, aes(color = NULL)) +
  geom_line(stat='smooth', method = "lm", size = 1.5) +
  ylim(0, 3) +
  facet_wrap(~ model) +
  labs(title = "n_i=1000, maf=0.1, r=1.5, phi=3, theta=30") +
  theme_bw() +
  scale_color_npg()
