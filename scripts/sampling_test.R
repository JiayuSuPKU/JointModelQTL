library(tidyverse)
library(ggsci)
library(ggrepel)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan_models <- load_models()

### J = 1, K = 1, L = 1

Y <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 100, baseline = 3, r = 0.5
)

results <- estimateCisRegEffects(
  data = Y$data, stan_models = stan_models,
  model = "joint", method = "sampling",
  significance_test = F
)

## sensitivity of R estimation

# helper function to test the sensitivity of the estimation of R posterior
sensitivityRPosterior <- function(weighted_ase = FALSE, phasing_error = "null",
                                  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
                                  phi = 3, theta = 30, baseline = 3, r = 0.5, ...) {

  Y <- simulateCisEffect.s2s(
    n_i = n_i, maf = maf, prob_ref = prob_ref, prob_as = prob_as,
    phi = phi, theta = theta, baseline = baseline, r = r, ...)

  trc <- estimateCisRegEffects(
    data = Y$data, stan_models = stan_models,
    model = "trc", method = "sampling",
    return_posterior = "full", significance_test = F,
    weighted_ase = weighted_ase, phasing_error = phasing_error
  )

  joint <- estimateCisRegEffects(
    data = Y$data, stan_models = stan_models,
    model = "joint", method = "sampling",
    return_posterior = "full", significance_test = F,
    weighted_ase = weighted_ase, phasing_error = phasing_error
  )

  ase <- estimateCisRegEffects(
    data = Y$data, stan_models = stan_models,
    model = "ase", method = "sampling",
    return_posterior = "full", significance_test = F,
    weighted_ase = weighted_ase, phasing_error = phasing_error
  )

  return(rbind(joint, trc, ase))
}

# sensitivity of r_est on theta (ase variance)
r_est_theta <- data.frame()
for (theta in c(100, 30, 10, 3)) {
  r <- sensitivityRPosterior(
    n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
    phi = 3, theta = theta, baseline = 3, r = 0.5
  )

  r <- r %>% mutate(theta = theta)
  r_est_theta <- rbind(r_est_theta, r)
}

r_est_theta %>% ggplot(aes(x = r_est, group = model, color = model)) +
  facet_wrap(~theta, labeller = label_both) +
  geom_vline(xintercept = 0.5, color = "black") +
  geom_density(size = 1) +
  labs(title = "n=1000, maf=0.1, b=3, phi=3") +
  theme_bw() +
  scale_color_npg()


# sensitivity of r_est on sample size and maf
r_est_ni_maf <- data.frame()
for (n_i in c(100, 300, 1000)) {
  for (maf in c(0.05, 0.1, 0.3)) {
    r <- sensitivityRPosterior(
      n_i = n_i, maf = maf, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5
    )
    r <- r %>% mutate(n_i = n_i, maf = maf)

    r_est_ni_maf <- rbind(r_est_ni_maf, r)
  }
}

r_est_ni_maf %>% ggplot(aes(x = r_est, group = model, color = model)) +
  facet_grid(
    maf ~ n_i,
    labeller = labeller(maf = label_both, n_i = label_both)
  ) +
  geom_vline(xintercept = 0.5, color = "black") +
  geom_density(size = 1) +
  labs(title = "b=3, phi=3, theta=30") +
  theme_bw() +
  scale_color_npg()


# sensitivity of r_est on cis-QTL effect size
r_est_r <- data.frame()
for (r_true in c(0, 0.2, 0.4, 0.6, 0.8, 1)) {
  r <- sensitivityRPosterior(
    n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
    phi = 3, theta = 30, baseline = 3, r = r_true
  )
  r <- r %>% mutate(r_true = r_true)

  r_est_r <- rbind(r_est_r, r)
}

r_est_r %>% ggplot(aes(x = r_est, group = model, color = model)) +
  facet_wrap(~r_true, labeller = labeller(r_true = label_both), nrow = 2) +
  geom_vline(aes(xintercept = r_true), color = "black") +
  geom_density(size = 1) +
  labs(title = "n_i=1000, maf=0.1, b=3, phi=3, theta=30") +
  theme_bw() +
  scale_color_npg()


# sensitivity of r_est on baseline expression
r_est_b <- data.frame()
for (b in c(0, 1, 3, 5, 7, 9)) {
  r <- sensitivityRPosterior(
    n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
    phi = 3, theta = 30, baseline = b, r = 0.5
  )
  r <- r %>% mutate(b = b)

  r_est_b <- rbind(r_est_b, r)
}

r_est_b %>% ggplot(aes(x = r_est, group = model, color = model)) +
  facet_wrap(~b, labeller = labeller(b = label_both)) +
  geom_vline(xintercept = 0.5, color = "black") +
  geom_density(size = 1) +
  labs(title = "n_i=1000, maf=0.1, phi=3, theta=30") +
  theme_bw() +
  scale_color_npg()


# weighted ase for low baseline expression and small effect size
r_est_weighted <- data.frame()

for (b in c(1, 3)) {
  for (r in c(0.1, 0.25, 0.4, 0.55)) {
    r_normal <- sensitivityRPosterior(
      n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = b, r = r
    )
    r_normal <- r_normal %>% mutate(weighted = 0, baseline = b, r = r)

    r_w <- sensitivityRPosterior(
      weighted_ase = T,
      n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = b, r = r
    )
    r_w <- r_w %>% mutate(weighted = 1, baseline = b, r = r)

    r_est_weighted <- rbind(r_est_weighted, r_normal, r_w)
  }
}

g1 <- r_est_weighted %>% filter(weighted == 0) %>%
  ggplot(aes(x = r_est, group = model, color = model)) +
  facet_grid(r ~ baseline,
    labeller = labeller(baseline = label_both, r = label_both)
  ) +
  geom_vline(aes(xintercept = r), color = "black") +
  geom_density(size = 1) +
  labs(title = "Basic models (n_i=1000, maf=0.1, phi=3, theta=30)") +
  theme_bw() +
  scale_color_npg()

g2 <- r_est_weighted %>% filter(weighted == 1) %>%
  ggplot(aes(x = r_est, group = model, color = model)) +
  facet_grid(r ~ baseline,
             labeller = labeller(baseline = label_both, r = label_both)
  ) +
  geom_vline(aes(xintercept = r), color = "black") +
  geom_density(size = 1) +
  labs(title = "Weighted models") +
  theme_bw() +
  scale_color_npg()

cowplot::plot_grid(g1, g2)


## runtime analysis

# helper function for testing runtime under multiple scenarios
getRunTime <- function(n_j, n_i = 1000, maf = 0.1, prob_ref = 0.5,
                       gene_pars = list(
                         list(prob_as = 0.5, phi = 3, theta = 30, baseline = 3, r = 0.5)
                       ), p_error = NULL) {
  if (n_j == 1) {
    Y <- simulateCisEffect.s2s(
      n_i = n_i, maf = maf, prob_ref = prob_ref,
      prob_as = gene_pars[[1]][["prob_as"]], phi = gene_pars[[1]][["phi"]],
      theta = gene_pars[[1]][["theta"]], baseline = gene_pars[[1]][["baseline"]],
      r = gene_pars[[1]][["r"]], p_error = p_error
    )

    fit_trc <- sampling(stan_models$basic$trc, data = Y$data, refresh = 0)
    fit_trcase <- sampling(stan_models$basic$joint, data = Y$data, refresh = 0)
  } else {
    Y <- simulateCisEffect.m2s(
      n_i = n_i, n_j = n_j, maf = maf, prob_ref = prob_ref,
      gene_pars = gene_pars
    )

    fit_trc <- stan(file = "./src/deprecated_stan_models/trc_multi_genes.stan", data = Y$data)
    fit_trcase <- stan(file = "./src/deprecated_stan_models/joint_multi_genes.stan", data = Y$data)
  }

  df.time <- rbind(get_elapsed_time(fit_trc), get_elapsed_time(fit_trcase)) %>%
    data.frame(row.names = NULL) %>%
    mutate(model = rep(c("trc", "trcase"), each = 4), n_j = n_j)
}

# n_i (J = 1, K = 1, L = 1)
runtime_n_i_100 <- getRunTime(n_j = 1, n_i = 100)
runtime_n_i_400 <- getRunTime(n_j = 1, n_i = 400)
runtime_n_i_1600 <- getRunTime(n_j = 1, n_i = 1600)
runtime_n_i_3200 <- getRunTime(n_j = 1, n_i = 3200)
runtime_n_i_6400 <- getRunTime(n_j = 1, n_i = 6400)

runtime_n_i <- rbind(
  runtime_n_i_100, runtime_n_i_400, runtime_n_i_1600,
  runtime_n_i_3200, runtime_n_i_6400
) %>% mutate(n_i = rep(c(100, 400, 1600, 3200, 6400), each = 8))

runtime_n_i.label <- runtime_n_i %>%
  group_by(model, n_i) %>%
  summarise(t = mean(warmup + sample))

runtime_n_i %>% ggplot(aes(x = n_i, y = warmup + sample, color = model)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x) +
  geom_label_repel(
    data = runtime_n_i.label, max.overlaps = 30,
    aes(x = n_i, y = t, label = round(t, 2))
  ) +
  scale_color_npg() +
  labs(x = "Number of samples", y = "Run time (warmup + sample, seconds)", color = "Model") +
  theme_bw()

# n_j (J > 1, K = 1, L = 1)

runtime_n_j_1 <- getRunTime(n_j = 1)
runtime_n_j_2 <- getRunTime(n_j = 2)
runtime_n_j_4 <- getRunTime(n_j = 4)
runtime_n_j_16 <- getRunTime(n_j = 16)
runtime_n_j_32 <- getRunTime(n_j = 32)
runtime_n_j_64 <- getRunTime(n_j = 64)

runtime_n_j <- rbind(
  runtime_n_j_1, runtime_n_j_2, runtime_n_j_4,
  runtime_n_j_16, runtime_n_j_32, runtime_n_j_64
)
runtime_n_j.label <- runtime_n_j %>%
  group_by(model, n_j) %>%
  summarise(t = mean(warmup + sample))

runtime_n_j %>% ggplot(aes(x = n_j, y = warmup + sample, color = model)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x) +
  geom_label_repel(
    data = runtime_n_j.label, max.overlaps = 30,
    aes(x = n_j, y = t, label = round(t, 2))
  ) +
  scale_color_npg() +
  labs(x = "Number of genes", y = "Run time (warmup + sample, seconds)", color = "Model") +
  theme_bw()


##  models that include phasing error

# runtime analysis

Y <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3, r = 0.5, p_error = 0.1
)

# basic model
fit_ase <- sampling(stan_models$basic$ase, data = Y$data)
# mixture model with known p_error for each sample
fit_ase_mix_known <- sampling(stan_models$phasing_error$ase_mixture_known, data = Y$data)
# mixture model with unknown p_error for each sample
fit_ase_mix_unknown <- sampling(stan_models$phasing_error$ase_mixture_unknown, data = Y$data)
# expectation model with known p_error for each sample
fit_ase_expect_known <- stan(
  file = "./src/deprecated_stan_models/ase_expect_known.stan",
  data = Y$data
)

runtime_p_error <- lapply(
  list(fit_ase, fit_ase_mix_known, fit_ase_mix_unknown, fit_ase_expect_known),
  get_elapsed_time
) %>%
  do.call(what = "rbind") %>%
  as.data.frame() %>%
  mutate(model = rep(c(
    "basic", "mix_known_p_error", "mix_unknown_p_error",
    "expect_known_p_error"
  ), each = 4)) %>%
  group_by(model) %>%
  summarize(warmup = mean(warmup), sample = mean(sample)) %>%
  pivot_longer(!model, names_to = "stage", values_to = "runtime")

runtime_p_error %>% ggplot(aes(x = model, y = runtime, fill = stage)) +
  geom_bar(position = "stack", stat = "identity", width = 0.5) +
  labs(x = "ASE models", y = "Runtime (s)") +
  scale_fill_npg() +
  theme_bw()


# sensitivity of r_est on phasing error
r_est_p_error <- data.frame()

for (p_error in c(0.05, 0.15, 0.25)) {
  for (r in c(0.3, 0.6, 0.9)) {
    r_normal <- sensitivityRPosterior(p_error = p_error, r = r)
    r_normal <- r_normal %>% mutate(mixture = "null", mean_p_error = p_error, r = r)

    r_p_k <- sensitivityRPosterior(phasing_error = "known", p_error = p_error, r = r)
    r_p_k <- r_p_k %>% mutate(mixture = "known", mean_p_error = p_error, r = r)

    r_p_u <- sensitivityRPosterior(phasing_error = "unknown", p_error = p_error, r = r)
    r_p_u <- r_p_u %>% mutate(mixture = "unknown", mean_p_error = p_error, r = r)

    r_est_p_error <- rbind(r_est_p_error, r_normal, r_p_k, r_p_u)
  }
}

g3 <- r_est_p_error %>% filter(mixture == "null") %>%
  ggplot(aes(x = r_est, group = model, color = model)) +
  facet_grid(mean_p_error ~ r,
    labeller = labeller(mean_p_error = label_both, r = label_both)
  ) +
  geom_vline(aes(xintercept = r), color = "black") +
  geom_density(size = 1) +
  labs(title = "Basic Models (n_i=1000, maf=0.1, b=3, phi=3, theta=30)") +
  theme_bw() +
  scale_color_npg()

g4 <- r_est_p_error %>% filter(mixture == "known") %>%
  ggplot(aes(x = r_est, group = model, color = model)) +
  facet_grid(mean_p_error ~ r,
             labeller = labeller(mean_p_error = label_both, r = label_both)
  ) +
  geom_vline(aes(xintercept = r), color = "black") +
  geom_density(size = 1) +
  labs(title = "Mixture Models (known)") +
  theme_bw() +
  scale_color_npg()

g5 <- r_est_p_error %>% filter(mixture == "unknown") %>%
  ggplot(aes(x = r_est, group = model, color = model)) +
  facet_grid(mean_p_error ~ r,
             labeller = labeller(mean_p_error = label_both, r = label_both)
  ) +
  geom_vline(aes(xintercept = r), color = "black") +
  geom_density(size = 1) +
  labs(title = "Mixture Models (unknown)") +
  theme_bw() +
  scale_color_npg()

cowplot::plot_grid(g3, g4, g5, nrow = 1)

