library(tidyverse)
library(ggsci)
library(ggrepel)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## random effects

# add individual effects
Y.0 <- simulateCisEffect.s2s(
  n_i = 100, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3, r = 1.5,
  origin = NULL, origin.effect = NULL
)
Y.1 <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3, r = 1.5,
  origin = sample(rep(1:100, each = 10), 1000, replace = F), origin.effect = rnorm(100)
)

fit_trc.0 <- stan(file = "./src/stan_models/lognorm_trc.stan", data = Y.0$data)
fit_trc.1 <- stan(file = "./src/stan_models/lognorm_trc.stan", data = Y.1$data)
fit_trc_indiv_re <- stan(file = "./src/stan_models/lognorm_confounding/trc_indiv_re.stan", data = Y.1$data)

data.frame(
  r_est = c(
    extract(fit_trc.0, "r")[[1]],
    extract(fit_trc.1, "r")[[1]],
    extract(fit_trc_indiv_re, "r")[[1]]
  ),
  model = rep(c("n_i=100", "n_i=1000, n_indiv=100", "trc_indiv_re"), each = 4000)
) %>% ggplot(aes(x = r_est, group = model, color = model)) +
  geom_vline(xintercept = 1.5, color = "black") +
  geom_density(size = 1) +
  labs(title = "n=1000, maf=0.1, b=3, phi=3") +
  scale_color_npg() +
  theme_bw()

estimateR <- function(n_indiv, n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
                      phi = 3, theta = 30, baseline = 3, r = 1.5) {
  Y.0 <- simulateCisEffect.s2s(
    n_i = n_indiv, maf = maf, prob_ref = prob_ref, prob_as = prob_as,
    phi = phi, theta = theta, baseline = baseline, r = r
  )
  Y.1 <- simulateCisEffect.s2s(
    n_i = n_i, maf = maf, prob_ref = prob_ref, prob_as = prob_as,
    phi = phi, theta = theta, baseline = baseline, r = r,
    origin = sample(rep(1:n_indiv, each = 1000 / n_indiv), 1000, replace = F),
    origin.effect = rnorm(n_indiv)
  )

  fit_trc.0 <- stan(file = "./src/stan_models/lognorm_trc.stan", data = Y.0$data)
  fit_trc.1 <- stan(file = "./src/stan_models/lognorm_trc.stan", data = Y.1$data)
  fit_trc_indiv_re <- stan(file = "./src/stan_models/lognorm_confounding/trc_indiv_re.stan", data = Y.1$data)

  out <- data.frame(
    r_est = c(
      extract(fit_trc.0, "r")[[1]],
      extract(fit_trc.1, "r")[[1]],
      extract(fit_trc_indiv_re, "r")[[1]]
    ),
    model = rep(c("No replicates", "Replicates", "Replicates + random effect"), each = 4000),
    r = r,
    n_indiv = n_indiv
  )

  return(out)
}

r_est_n_indiv <- lapply(
  c(10, 50, 100, 200, 500, 1000), estimateR
) %>%
  do.call(what = "rbind")

r_est_n_indiv %>% ggplot(aes(x = r_est, group = model, color = model)) +
  facet_wrap(~n_indiv, labeller = label_both) +
  xlim(0, 3) +
  geom_vline(xintercept = 1.5, color = "black") +
  geom_density(size = 1) +
  labs(title = "n=1000, maf=0.1, b=3, r=1.5, phi=3") +
  scale_color_aaas() +
  theme_bw()

# runtime analysis
runtime_n_indiv <- data.frame()

for (n_indiv in c(10, 50, 100, 200, 500, 1000)) {
  Y <- simulateCisEffect.s2s(
    n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
    phi = 3, theta = 100, baseline = 3, r = 1.5,
    origin = sample(rep(1:n_indiv, each = 1000 / n_indiv), 1000, replace = F),
    origin.effect = rnorm(n_indiv)
  )

  fit_trc_indiv_re <- stan(
    file = "./src/stan_models/lognorm_confounding/trc_indiv_re.stan", data = Y$data
  )

  runtime_n_indiv <- rbind(runtime_n_indiv, get_elapsed_time(fit_trc_indiv_re) %>% colMeans())
}

runtime_n_indiv <- runtime_n_indiv %>%
  `names<-`(c("warmup", "sample")) %>%
  mutate(n_indiv = c(10, 50, 100, 200, 500, 1000))

runtime_n_indiv %>%
  pivot_longer(!n_indiv, names_to = "stage", values_to = "runtime") %>%
  ggplot(aes(x = n_indiv, y = runtime, color = stage)) +
  geom_point() +
  geom_line() +
  scale_color_aaas() +
  labs(
    title = "Fixed total sample size (n=1000)",
    x = "Number of individuals", y = "Run time (seconds)", color = "Stage"
  ) +
  theme_bw()

## re through covariance structure

Y.1 <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3, r = 1,
  origin = sample(rep(1:100, each = 10), 1000, replace = F), origin.effect = rnorm(100)
)

# 1000 transitions would take 3100s, too slow
# fit_trc_indiv_suez <- stan(file = "./src/stan_models/lognorm_confounding/trc_indiv_suez.stan", data = Y.1$data)
trc <- stan_model(file = "./src/stan_models/lognorm_trc.stan")
trc_indiv_re <- stan_model(file = "./src/stan_models/lognorm_confounding/trc_indiv_re.stan")
trc_indiv_suez <- stan_model(file = "./src/stan_models/lognorm_confounding/trc_indiv_suez.stan")

lbfgs_trc <- optimizing(
  trc,
  data = Y.1$data,
  algorithm = "LBFGS", verbos = F, as_vector = F
)
lbfgs_trc_indiv_re <- optimizing(
  trc_indiv_re,
  data = Y.1$data,
  algorithm = "LBFGS", verbos = F, as_vector = F
)
lbfgs_trc_indiv_suez <- optimizing(
  trc_indiv_suez,
  data = Y.1$data,
  algorithm = "LBFGS", verbos = F, as_vector = F
)

# runtime comparisons
system.time({
  optimizing(
    trc,
    data = Y.1$data, algorithm = "LBFGS"
  )
})
system.time({
  optimizing(
    trc_indiv_re,
    data = Y.1$data,
    algorithm = "LBFGS", verbos = F, as_vector = F
  )
})
system.time({
  optimizing(
    trc_indiv_suez,
    data = Y.1$data,
    algorithm = "LBFGS", verbos = F, as_vector = F
  )
})

p1 = data.frame(
  t = c(0.020, 0.025, 4.834, 0.005, 0.071, 4.725, 0.013, 0.025, 4.297, 0.004, 0.022, 5.485),
  model = rep(c("trc", "trc_re", "trc_cov"), times = 4)
) %>% mutate(model = factor(model, c("trc", "trc_re", "trc_cov"))) %>%
  group_by(model) %>%
  mutate(avg_t = mean(t)) %>%
  ggplot(aes(x = model)) +
  geom_bar(aes(y = avg_t, fill = model), stat = "identity", position=position_dodge()) +
  geom_point(aes(y = t)) +
  geom_text(aes(y = avg_t + .5, label = avg_t)) +
  scale_fill_aaas() +
  labs(
    title = "L-BFGS Runtime (n = 1000, n_indiv = 100)",
    x = "", y = "Elapsed Time (s)") +
  theme_bw() + theme(legend.position = "none")


# estimation comparison
estimateR.lbfgs <- function(n_indiv = 100, n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
                            phi = 3, theta = 30, baseline = 3, r = 1.5){

  Y.0 <- simulateCisEffect.s2s(
    n_i = n_indiv, maf = maf, prob_ref = prob_ref, prob_as = prob_as,
    phi = phi, theta = theta, baseline = baseline, r = r
  )
  Y.1 <- simulateCisEffect.s2s(
    n_i = n_i, maf = maf, prob_ref = prob_ref, prob_as = prob_as,
    phi = phi, theta = theta, baseline = baseline, r = r,
    origin = sample(rep(1:n_indiv, each = 1000 / n_indiv), 1000, replace = F),
    origin.effect = rnorm(n_indiv)
  )

  f_trc.0 <- optimizing(trc, data = Y.0$data, algorithm = "LBFGS", verbos = F, as_vector = F)
  f_trc.1 <- optimizing(trc, data = Y.1$data, algorithm = "LBFGS", verbos = F, as_vector = F)
  f_trc_indiv_re <- optimizing(trc_indiv_re, data = Y.1$data, algorithm = "LBFGS", verbos = F, as_vector = F)
  f_trc_indiv_suez <- optimizing(trc_indiv_suez, data = Y.1$data, algorithm = "LBFGS", verbos = F, as_vector = F)

  out <- data.frame(
    r_est = sapply(
      list(f_trc.0, f_trc.1, f_trc_indiv_re, f_trc_indiv_suez),
      function(x) x$par$r
    ),
    aic = sapply(
      list(f_trc.0, f_trc.1, f_trc_indiv_re, f_trc_indiv_suez),
      function(x) -2*x$value
    ) + 2 * c(0, 0, n_indiv, 1),
    model = c("No replicates", "Replicates", "Replicates + random effect", "Replicates + covariance"),
    r = r,
    n_indiv = n_indiv
  )

  return(out)
}

# compare aic
aic_lbfgs <- data.frame()
for (n_indiv in c(10, 50, 100, 200, 500, 1000)){
  aic_lbfgs <- rbind(aic_lbfgs, estimateR.lbfgs(n_indiv = n_indiv, r = 1))
}

p2 = aic_lbfgs %>%
  mutate(model = factor(
    model, levels = c("No replicates", "Replicates",
                      "Replicates + random effect", "Replicates + covariance")
    )
  ) %>%
  ggplot(aes(x = n_indiv, y = aic, color = model, group = model)) +
  geom_point() + geom_line() +
  scale_color_manual(values = c("gray", pal_aaas()(3))) +
  ylim(-100, 2000) +
  labs(title = "L-BFGS Model Selection (n = 1000)", y = "AIC") +
  theme_bw()

cowplot::plot_grid(p1, p2, rel_widths = c(2, 3))

# compare r_est
r_est_lbfgs <- data.frame()

for (n_indiv in c(10, 50, 100, 200, 500, 1000)){
  for (r in seq(0, 1, 0.05)){
    r_est_lbfgs <- rbind(r_est_lbfgs, estimateR.lbfgs(n_indiv = n_indiv, r = r))
  }
}

r_est_lbfgs %>%
  ggplot(aes(x = r, y = r_est, color = model, group = model)) +
  facet_wrap(~ n_indiv, labeller = label_both) +
  geom_point() +
  ylim(-2, 2) +
  geom_ribbon(stat='smooth', method = "lm", se=TRUE, alpha=0.15, aes(color = NULL)) +
  geom_line(stat='smooth', method = "lm", size = 1.5) +
  scale_color_npg() +
  labs(
    title = "Fixed total sample size (n=1000)",
    x = "True R", y = "Estimated R", color = "Models"
  ) +
  theme_bw()


## fixed effects
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
    extract(fit_trc_libsize, "r")[[1]]
  ),
  model = rep(c("baseline", "trc_x"), each = 4000)
) %>% ggplot(aes(x = r_est, group = model, color = model)) +
  geom_vline(xintercept = 1.5, color = "black") +
  geom_density(size = 1) +
  labs(title = "n_confounder=1\nn_sample=1000, maf=0.1, b=3, phi=3, theta=30") +
  theme_bw()

# add multiple confounding factors
confounder <- cbind(
  rbinom(1000, 1, prob = 0.1),
  rbinom(1000, 1, prob = 0.3),
  rbinom(1000, 1, prob = 0.5),
  rnorm(1000),
  rnorm(1000, 0, 3),
  rnorm(1000, 1, 1),
  sample(1:10, size = 1000, replace = T),
  runif(1000)
)
confounder.effect <- rep(1, 8)

Y.multi <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 100, baseline = 3, r = 1.5,
  confounder = confounder, confounder.effect = confounder.effect
)

fit_trc <- stan(file = "./src/stan_models/lognorm_trc.stan", data = Y.multi$data)
fit_trc_x <- stan(
  file = "./src/stan_models/lognorm_confounding/trc_multi_cov_fe.stan",
  data = Y.multi$data
)
traceplot(fit_trc_x, pars = "beta")

g2 <- data.frame(
  r_est = c(
    extract(fit_trc, "r")[[1]],
    extract(fit_trc_x, "r")[[1]]
  ),
  model = rep(c("baseline", "trc_x"), each = 4000)
) %>% ggplot(aes(x = r_est, group = model, color = model)) +
  geom_vline(xintercept = 1.5, color = "black") +
  geom_density(size = 1) +
  labs(title = "n_confounder=10\nn_sample=1000, maf=0.1, b=3, phi=3, theta=30") +
  theme_bw()

# runtime analysis

runtime_n <- data.frame()
for (n in seq(1, 50, 5)) {
  Y <- simulateCisEffect.s2s(
    n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
    phi = 3, theta = 100, baseline = 3, r = 1.5,
    confounder = matrix(rnorm(1000 * n), ncol = n), confounder.effect = rnorm(n)
  )

  fit_trc_x <- stan(
    file = "./src/stan_models/lognorm_confounding/trc_fixed_effects.stan",
    data = Y$data
  )

  runtime_n <- rbind(runtime_n, get_elapsed_time(fit_trc_x) %>% colMeans())
}

runtime_n <- runtime_n %>%
  `names<-`(c("warmup", "sample")) %>%
  mutate(n = seq(1, 50, 5))

g3 <- runtime_n %>%
  pivot_longer(!n, names_to = "stage", values_to = "runtime") %>%
  ggplot(aes(x = n, y = runtime, color = stage)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x) +
  scale_color_aaas() +
  labs(x = "Number of confounders (fixed effect)", y = "Run time (seconds)", color = "Stage") +
  theme_bw()

cowplot::plot_grid(g1, g2, g3, nrow = 1)
