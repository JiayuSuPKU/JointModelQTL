source("./R/simulation.R")

library(rstan)
library(dplyr)
library(ggsci)
library(ggrepel)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

Y <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 100, baseline = 3, r = 1.5
)

fit_trc <- stan(file = "./src/stan_models/lognorm_trc.stan", data = Y$data)
fit_trcase <- stan(file = "./src/stan_models/lognorm_trcase.stan", data = Y$data)

estimateR <- function(...) {
  # helper function for testing under multiple scenarios

  Y <- simulateCisEffect.s2s(...)

  fit_trc <- stan(file = "./src/stan_models/lognorm_trc.stan", data = Y$data)
  fit_trcase <- stan(file = "./src/stan_models/lognorm_trcase.stan", data = Y$data)

  r_est_trc <- extract(fit_trc, "r")[[1]]
  r_est_trcase <- extract(fit_trcase, "r")[[1]]

  return(
    data.frame(
      r_est = c(r_est_trc, r_est_trcase),
      model = rep(c("trc", "trcase"), each = length(r_est_trc))
    )
  )
}

# theta

r_est_theta <- data.frame()

for (theta in c(100, 30, 10, 3)) {
  r <- estimateR(
    n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
    phi = 3, theta = theta, baseline = 3, r = 1.5
  )

  r <- r %>% mutate(theta = theta)
  r_est_theta <- rbind(r_est_theta, r)
}

r_est_theta %>% ggplot(aes(x = r_est, group = model, color = model)) +
  facet_wrap(~theta, labeller = labeller(r_true = label_both)) +
  geom_vline(xintercept = 1.5, color = "black") +
  geom_density(size = 1) +
  labs(title = "n=1000, maf=0.1, b=3, phi=3") +
  theme_bw() +
  scale_color_npg()

# sample size and maf

r_est_ni_maf <- data.frame()

for (n_i in c(100, 300, 1000)) {
  for (maf in c(0.05, 0.1, 0.3)) {
    r <- estimateR(
      n_i = n_i, maf = maf, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 1.5
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
  geom_vline(xintercept = 1.5, color = "black") +
  geom_density(size = 1) +
  labs(title = "b=3, phi=3, theta=30") +
  theme_bw() +
  scale_color_npg()

# effect size

r_est_r <- data.frame()

for (r_true in c(0.75, 1, 1.25, 1.5, 1.75, 2)) {
  r <- estimateR(
    n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
    phi = 3, theta = 30, baseline = 3, r = r_true
  )
  r <- r %>% mutate(r_true = r_true)

  r_est_r <- rbind(r_est_r, r)
}

r_est_r %>% ggplot(aes(x = r_est, group = model, color = model)) +
  facet_wrap(~r_true, labeller = labeller(r_true = label_both)) +
  geom_vline(aes(xintercept = r_true), color = "black") +
  geom_density(size = 1) +
  labs(title = "n_i=1000, maf=0.1, b=3, phi=3, theta=30") +
  theme_bw() +
  scale_color_npg()

# baseline expression

r_est_b <- data.frame()

for (b in c(0, 1, 3, 5, 7, 9)) {
  r <- estimateR(
    n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
    phi = 3, theta = 30, baseline = b, r = 1.5
  )
  r <- r %>% mutate(b = b)

  r_est_b <- rbind(r_est_b, r)
}

r_est_b %>% ggplot(aes(x = r_est, group = model, color = model)) +
  facet_wrap(~b, labeller = labeller(b = label_both)) +
  geom_vline(xintercept = 1.5, color = "black") +
  geom_density(size = 1) +
  labs(title = "n_i=1000, maf=0.1, r=1.5, phi=3, theta=30") +
  theme_bw() +
  scale_color_npg()


# weighted ase

estimateR.weighted <- function(n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
                               phi = 3, theta = 30, baseline = 1, r = 1.5) {
  Y <- simulateCisEffect.s2s(
    n_i = n_i, maf = maf, prob_ref = prob_ref, prob_as = prob_as,
    phi = phi, theta = theta, baseline = baseline, r = r
  )

  fit_ase <- stan(file = "./src/stan_models/lognorm_ase.stan", data = Y$data)
  fit_weighted_ase <- stan(
    file = "./src/stan_models/lognorm_weighted_ase.stan", data = Y$data
  )

  return(
    data.frame(
      r_est = c(extract(fit_ase, "r")[[1]], extract(fit_weighted_ase, "r")[[1]]),
      model = rep(c("ase", "weighted ase"), each = 4000),
      baseline = rep(baseline, 8000),
      r = rep(r, 8000)
    )
  )
}

r_est_weighted <- data.frame()

for (baseline in c(1, 3, 5)){
  for (r in c(1.1, 1.3, 1.5)){
    r_est_weighted <- rbind(
      r_est_weighted, estimateR.weighted(baseline = baseline, r = r))
  }
}

r_est_weighted %>% ggplot(aes(x = r_est, group = model, color = model)) +
  facet_grid(baseline ~ r,
             labeller = labeller(baseline = label_both, r = label_both)) +
  geom_vline(aes(xintercept = r), color = "black") +
  geom_density(size = 1) +
  labs(title = "n_i=1000, maf=0.1, r=1.5, phi=3, theta=30") +
  theme_bw() +
  scale_color_npg()

# runtime analysis

estimateTime <- function(n_j, n_i = 1000, maf = 0.1, prob_ref = 0.5,
                         gene_pars = list(
                           list(prob_as = 0.5, phi = 3, theta = 30, baseline = 3, r = 1.5)
                         ), p_error = NULL) {
  # helper function for testing under multiple scenarios

  if (n_j == 1) {
    Y <- simulateCisEffect.s2s(
      n_i = n_i, maf = maf, prob_ref = prob_ref,
      prob_as = gene_pars[[1]][["prob_as"]], phi = gene_pars[[1]][["phi"]],
      theta = gene_pars[[1]][["theta"]], baseline = gene_pars[[1]][["baseline"]],
      r = gene_pars[[1]][["r"]], p_error = p_error
    )

    fit_trc <- stan(file = "./src/stan_models/lognorm_trc.stan", data = Y$data)
    fit_trcase <- stan(file = "./src/stan_models/lognorm_trcase.stan", data = Y$data)
  } else {
    Y <- simulateCisEffect.m2s(
      n_i = n_i, n_j = n_j, maf = maf, prob_ref = prob_ref,
      gene_pars = gene_pars
    )

    fit_trc <- stan(file = "./src/stan_models/lognorm_trc_j.stan", data = Y$data)
    fit_trcase <- stan(file = "./src/stan_models/lognorm_trcase_j.stan", data = Y$data)
  }

  df.time <- rbind(get_elapsed_time(fit_trc), get_elapsed_time(fit_trcase)) %>%
    data.frame(row.names = NULL) %>%
    mutate(model = rep(c("trc", "trcase"), each = 4), n_j = n_j)
}

# J = 1, K = 1, L = 1
runtime_n_i_100 <- estimateTime(n_j = 1, n_i = 100)
runtime_n_i_400 <- estimateTime(n_j = 1, n_i = 400)
runtime_n_i_1600 <- estimateTime(n_j = 1, n_i = 1600)
runtime_n_i_3200 <- estimateTime(n_j = 1, n_i = 3200)
runtime_n_i_6400 <- estimateTime(n_j = 1, n_i = 6400)

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

# J > 1, K = 1, L = 1

Y <- simulateCisEffect.m2s(
  n_i = 1000, n_j = 10, maf = 0.1, prob_ref = 0.5,
  gene_pars = list(list(prob_as = 0.5, phi = 3, theta = 30, baseline = 3, r = 1.5))
)

fit_trc <- stan(file = "./src/stan_models/lognorm_trc_j.stan", data = Y$data)
fit_trcase <- stan(file = "./src/stan_models/lognorm_trcase_j.stan", data = Y$data)

summary(fit_trcase)
get_elapsed_time(fit_trcase)

runtime_n_j_1 <- estimateTime(n_j = 1)
runtime_n_j_2 <- estimateTime(n_j = 2)
runtime_n_j_4 <- estimateTime(n_j = 4)
runtime_n_j_16 <- estimateTime(n_j = 16)
runtime_n_j_32 <- estimateTime(n_j = 32)
runtime_n_j_64 <- estimateTime(n_j = 64)

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

# phasing error

Y <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 100, baseline = 3, r = 1.5, p_error = 0.1
)

fit_ase <- stan(
  file = "./src/stan_models/lognorm_ase.stan",
  data = Y$data
)

fit_ase_mix_known <- stan(
  file = "./src/stan_models/lognorm_with_phasing_error/ase_mixture_known_error.stan",
  data = Y$data
) # p_error for each sample

fit_ase_mix_unknown <- stan(
  file = "./src/stan_models/lognorm_with_phasing_error/ase_mixture_unknown_error.stan",
  data = Y$data
) # not possible to have p_error for each sample

fit_ase_expect_known <- stan(
  file = "./src/stan_models/lognorm_with_phasing_error/ase_expect_known_error.stan",
  data = Y$data
) # p_error for each sample

fit_ase_expect_unknown <- stan(
  file = "./src/stan_models/lognorm_with_phasing_error/ase_expect_unknown_error.stan",
  data = Y$data
) # not possible to have p_error for each sample
