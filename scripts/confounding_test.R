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
fit_trc_libsize <- stan(file = "./src/stan_models/lognorm_confounding/trc_lib_size.stan", data = Y.2$data)

g1 <- data.frame(
  r_est = c(
    extract(fit_trc.2, "r")[[1]],
    extract(fit_trc_libsize, "r")[[1]]),
  model = rep(c("baseline", "trc_x"), each = 4000)
) %>% ggplot(aes(x = r_est, group = model, color = model)) +
  geom_vline(xintercept = 1.5, color = "black") +
  geom_density(size = 1) +
  labs(title = "n_confounder=1\nn_sample=1000, maf=0.1, b=3, phi=3, theta=30") +
  theme_bw()

# add multiple confounding factors
confounder = cbind(
  rbinom(1000, 1, prob = 0.1),
  rbinom(1000, 1, prob = 0.3),
  rbinom(1000, 1, prob = 0.5),
  rnorm(1000),
  rnorm(1000, 0, 3),
  rnorm(1000, 1, 1),
  sample(1:10, size = 1000, replace = T),
  runif(1000)
)
confounder.effect = rep(1, 8)

Y.multi <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 100, baseline = 3, r = 1.5,
  confounder = confounder, confounder.effect = confounder.effect
)

fit_trc <- stan(file = "./src/stan_models/lognorm_trc.stan", data = Y.multi$data)
fit_trc_x <- stan(
  file = "./src/stan_models/lognorm_confounding/trc_fixed_effects.stan",
  data = Y.multi$data
)
traceplot(fit_trc_x, pars = "beta")

g2 <- data.frame(
  r_est = c(
    extract(fit_trc, "r")[[1]],
    extract(fit_trc_x, "r")[[1]]),
  model = rep(c("baseline", "trc_x"), each = 4000)
) %>% ggplot(aes(x = r_est, group = model, color = model)) +
  geom_vline(xintercept = 1.5, color = "black") +
  geom_density(size = 1) +
  labs(title = "n_confounder=10\nn_sample=1000, maf=0.1, b=3, phi=3, theta=30") +
  theme_bw()

# runtime analysis

runtime_n <- data.frame()
for (n in seq(1, 50, 5)){
  Y <- simulateCisEffect.s2s(
    n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
    phi = 3, theta = 100, baseline = 3, r = 1.5,
    confounder = matrix(rnorm(1000*n), ncol = n), confounder.effect = rnorm(n)
  )

  fit_trc_x <- stan(
    file = "./src/stan_models/lognorm_confounding/trc_fixed_effects.stan",
    data = Y$data
  )

  runtime_n <- rbind(runtime_n, get_elapsed_time(fit_trc_x) %>% colMeans())
}

runtime_n <- runtime_n %>%
  `names<-`(c("warmup", "sample")) %>% mutate(n = seq(1, 50, 5))

g3 <- runtime_n %>%
  pivot_longer(!n, names_to = "stage", values_to = "runtime") %>%
  ggplot(aes(x = n, y = runtime, color = stage)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x) +
  scale_color_aaas() +
  labs(x = "Number of confounders (fixed effect)", y = "Run time (seconds)", color = "Stage") +
  theme_bw()

cowplot::plot_grid(g1, g2, g3, nrow = 1)
