library(tidyverse)
library(ggsci)
library(ggrepel)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

Y.0 <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 100, baseline = 3, r = 1.5,
  origin = NULL, origin.effect = NULL
)

# add individual effects
Y.1 <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 100, baseline = 3, r = 1.5,
  origin = sample(rep(1:100, each = 10), 1000, replace = F), origin.effect = rnorm(100)
)

fit_trc.0 <- stan(file = "./src/stan_models/lognorm_trc.stan", data = Y.0$data)
fit_trc.1 <- stan(file = "./src/stan_models/lognorm_trc.stan", data = Y.1$data)
fit_trc_indiv <- stan(file = "./src/stan_models/lognorm_confounding/trc_individual.stan", data = Y.1$data)

data.frame(
  r_est = c(
    extract(fit_trc.0, "r")[[1]],
    extract(fit_trc.1, "r")[[1]],
    extract(fit_trc_indiv, "r")[[1]]),
  model = rep(c("baseline", "+ individual effect", "trc_indiv"), each = 4000)
) %>% ggplot(aes(x = r_est, group = model, color = model)) +
  geom_vline(xintercept = 1.5, color = "black") +
  geom_density(size = 1) +
  labs(title = "n=1000, maf=0.1, b=3, phi=3") +
  theme_bw()

# add library size
Y.2 <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 100, baseline = 3, r = 1.5,
  lib.size = rnorm(1000)
)

fit_trc.2 <- stan(file = "./src/stan_models/lognorm_trc.stan", data = Y.2$data)

fit_trc_indiv <- stan(file = "./src/stan_models/lognorm_confounding/trc_individual.stan", data = Y.1$data)

