source("./R/simulation.R")

library(rstan)
library(dplyr)
library(ggsci)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

Y <- simulateCisEffectSingle(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 100, baseline = 3, r = 1.5)

fit_trc <- stan(file = "./src/stan_models/lognorm_trc.stan", data = Y$data)
fit_trcase <- stan(file = "./src/stan_models/lognorm_trcase.stan", data = Y$data)

estimateR <- function(...) {
  # helper function for testing under multiple scenarios

  Y <- simulateCisEffectSingle(...)

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


# J > 1, K = 1, L = 1

Y <- simulateCisEffectMulti(
  n_i = 1000, n_j = 10, maf = 0.1, prob_ref = 0.5,
  gene_pars = list(list(prob_as = 0.5, phi = 3, theta = 30, baseline = 3, r = 1.5)))

fit_trc <- stan(file = "./src/stan_models/lognorm_trc_j.stan", data = Y$data)
fit_trcase <- stan(file = "./src/stan_models/lognorm_trcase_j.stan", data = Y$data)

summary(fit_trcase)
