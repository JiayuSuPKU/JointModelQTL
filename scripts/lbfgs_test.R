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
  lognorm_ase,
  data = Y$data, algorithm = "LBFGS", as_vector = F,
  draw = 10, importance_resampling = T
)
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
  geom_ribbon(stat = "smooth", method = "lm", se = TRUE, alpha = 0.15, aes(color = NULL)) +
  geom_line(stat = "smooth", method = "lm", size = 1.5) +
  ylim(0, 3) +
  facet_wrap(~model) +
  labs(title = "n_i=1000, maf=0.1, r=1.5, phi=3, theta=30") +
  theme_bw() +
  scale_color_npg()

## likelihood ratio test
lognorm_trc <- stan_model(file = "./src/stan_models/lognorm_trc.stan")
lognorm_ase <- stan_model(file = "./src/stan_models/lognorm_ase.stan")
lognorm_trcase <- stan_model(file = "./src/stan_models/lognorm_trcase.stan")

lognorm_trc_null <- stan_model(file = "./src/stan_models/lognorm_lrt/trc_null.stan")
lognorm_ase_null <- stan_model(file = "./src/stan_models/lognorm_lrt/ase_null.stan")
lognorm_trcase_null <- stan_model(file = "./src/stan_models/lognorm_lrt/trcase_null.stan")

# main estimate function
estimateCisRegEffect <- function(data) {
  mle_trc <- optimizing(lognorm_trc, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)
  mle_ase <- optimizing(lognorm_ase, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)
  mle_trcase <- optimizing(lognorm_trcase, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)

  mle_trc_null <- optimizing(lognorm_trc_null, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)
  mle_ase_null <- optimizing(lognorm_ase_null, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)
  mle_trcase_null <- optimizing(lognorm_trcase_null, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)

  chi_sq <- -2 * c(
    mle_trc_null$value - mle_trc$value,
    mle_ase_null$value - mle_ase$value,
    mle_trcase_null$value - mle_trcase$value
  )
  df <- c(1, 1, 1)
  p_val <- 1 - pchisq(q = chi_sq, df = df)
  out <- data.frame(
    r_est = sapply(list(mle_trc, mle_ase, mle_trcase), function(x) x$par$r),
    model = c("trc", "ase", "trcase"),
    chi_sq = chi_sq,
    df = df,
    p_val = p_val
  )

  return(out)
}

Y <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3, r = 0.5
)
estimateCisRegEffect(Y$data)

# power analysis
design_power <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(b = c(1, 3, 5))
) %>%
  merge(y = data.frame(r = seq(0.01, 1.01, 0.05)))

power_test <- apply(design_power, 1, function(x) {
  lapply(1:200, function(i) {
    Y <- simulateCisEffect.s2s(
      n_i = x[1], baseline = x[2], r = x[3],
      maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
    )
    out <- estimateCisRegEffect(Y$data) %>%
      mutate(n_i = x[1], baseline = x[2], r_true = x[3])
    return(out)
  }) %>% do.call(what = "rbind")
}) %>% do.call(what = "rbind")


df_power <- power_test %>%
  group_by(n_i, baseline, r_true, model) %>%
  mutate(filter = p_val <= 0.05, power = sum(filter) / 200)

df_power %>% ggplot(aes(x = r_true, y = power, color = model, group = model)) +
  geom_point() +
  geom_line() +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  labs(title = "Power at 5% significance level", y = "Power", x = "R (log FC)") +
  scale_color_npg() +
  theme_bw()

# false positive analysis
design_fp <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(b = c(1, 3, 5))
)


fp_test <- apply(design_fp[rep(1:9, each = 100), ], 1, function(x) {
  lapply(1:100, function(i) {
    Y <- simulateCisEffect.s2s(
      n_i = x[1], baseline = x[2], r = 0,
      maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
    )
    estimateCisRegEffect(Y$data)
  }) %>%
    do.call(what = "rbind") %>%
    group_by(model) %>%
    summarise(fp = mean(p_val <= 0.05)) %>%
    mutate(n_i = x[1], baseline = x[2], r_true = 0)
}) %>% do.call(what = "rbind")

fp_test %>% ggplot(aes(x = model, y = fp, color = model)) +
  geom_hline(yintercept = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  labs(title = "False Positive Rate at 5% significance level", y = "False Positive", x = "Model") +
  scale_color_npg() +
  theme_bw()

## weighted ase
lognorm_ase_w <- stan_model(file = "./src/stan_models/lognorm_weighted_ase.stan")
lognorm_ase_w_null <- stan_model(file = "./src/stan_models/lognorm_lrt/ase_w_null.stan")

estimateCisRegEffect.weighted_ase <- function(data) {
  mle_ase <- optimizing(lognorm_ase, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)
  mle_ase_w <- optimizing(lognorm_ase_w, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)

  mle_ase_null <- optimizing(lognorm_ase_null, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)
  mle_ase_w_null <- optimizing(lognorm_ase_w_null, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)

  chi_sq <- -2 * c(
    mle_ase_null$value - mle_ase$value,
    mle_ase_w_null$value - mle_ase_w$value
  )
  df <- c(1, 1)
  p_val <- 1 - pchisq(q = chi_sq, df = df)
  out <- data.frame(
    r_est = sapply(list(mle_ase, mle_ase_w), function(x) x$par$r),
    model = c("ase", "weighted ase"),
    chi_sq = chi_sq,
    df = df,
    p_val = p_val
  )

  return(out)
}

# power
power_test.ase_w <- apply(design_power, 1, function(x) {
  lapply(1:200, function(i) {
    Y <- simulateCisEffect.s2s(
      n_i = x[1], baseline = x[2], r = x[3],
      maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
    )
    estimateCisRegEffect.weighted_ase(Y$data)
  }) %>%
    do.call(what = "rbind") %>%
    group_by(model) %>%
    summarise(power = mean(p_val <= 0.05)) %>%
    mutate(n_i = x[1], baseline = x[2], r_true = x[3])
}) %>% do.call(what = "rbind")


power_test.ase_w %>%
  ggplot(aes(x = r_true, y = power, color = model, group = model)) +
  geom_point() +
  geom_line() +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  labs(title = "Power at 5% significance level", y = "Power", x = "R (log FC)") +
  scale_color_aaas() +
  theme_bw()

# irregular case
set.seed(2021)
Y <- simulateCisEffect.s2s(
  n_i = 500, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3, r = 1
)
estimateCisRegEffect.weighted_ase(Y$data)

mle_ase <- optimizing(lognorm_ase, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
mle_ase_w <- optimizing(lognorm_ase_w, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)

fit_ase <- sampling(lognorm_ase, data = Y$data)
fit_ase_w <- sampling(lognorm_ase_w, data = Y$data)
data.frame(
  r_est = c(extract(fit_ase, pars = "r")[[1]], extract(fit_ase_w, pars = "r")[[1]]),
  model = rep(c("ase", "weighted ase"), each = 4000)
) %>%
  ggplot(aes(x = r_est, group = model, color = model)) +
  geom_density() +
  geom_vline(xintercept = 1) +
  scale_color_aaas() +
  labs(title = "Posterior distribution of R (n_i = 500, b = 3)") +
  theme_bw()

## phasing error
lognorm_ase_mix_known <- stan_model(file = "./src/stan_models/lognorm_with_phasing_error/ase_mixture_known_error.stan")
lognorm_ase_mix_unknown <- stan_model(file = "./src/stan_models/lognorm_with_phasing_error/ase_mixture_unknown_error.stan")

estimateCisRegEffect.p_error <- function(data) {
  mle_ase <- optimizing(lognorm_ase, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)
  mle_ase_mix_k <- optimizing(lognorm_ase_mix_known, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)
  mle_ase_mix_u <- optimizing(lognorm_ase_mix_unknown, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)

  mle_ase_null <- optimizing(lognorm_ase_null, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)

  chi_sq <- -2 * c(
    mle_ase_null$value - mle_ase$value,
    mle_ase_null$value - mle_ase_mix_k$value,
    mle_ase_null$value - mle_ase_mix_u$value
  )
  df <- c(1, 1, 1)
  p_val <- 1 - pchisq(q = chi_sq, df = df)
  out <- data.frame(
    r_est = sapply(list(mle_ase, mle_ase_mix_k, mle_ase_mix_u), function(x) x$par$r),
    model = c("ase", "ase_mix_known", "ase_mix_unknown"),
    chi_sq = chi_sq,
    df = df,
    p_val = p_val
  )

  return(out)
}

Y <- simulateCisEffect.s2s(
  n_i = 500, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3, r = 0.1, p_error = 0.1
)
estimateCisRegEffect.p_error(Y$data)

# power
design_power.p_error <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(p_error = c(0.05, 0.15, 0.25))
) %>%
  merge(y = data.frame(r = seq(0.01, 1.01, 0.05)))

power_test.p_error <- apply(design_power.p_error, 1, function(x) {
  lapply(1:200, function(i) {
    Y <- simulateCisEffect.s2s(
      n_i = x[1], p_error = x[2], r = x[3],
      baseline = 3, maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
    )
    estimateCisRegEffect.p_error(Y$data)
  }) %>%
    do.call(what = "rbind") %>%
    group_by(model) %>%
    summarise(power = mean(p_val <= 0.05)) %>%
    mutate(n_i = x[1], p_error = x[2], r_true = x[3])
}) %>% do.call(what = "rbind")

g2 <- power_test.p_error %>%
  ggplot(aes(x = r_true, y = power, color = model, group = model)) +
  geom_point() +
  geom_line() +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  labs(title = "Power at 5% significance level", y = "Power", x = "R (log FC)") +
  scale_color_aaas() +
  theme_bw()

# false positive
design_fp.p_error <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(p_error = c(0.05, 0.15, 0.25))
)

fp_test.p_error <- apply(design_fp.p_error[rep(1:9, each = 100), ], 1, function(x) {
  lapply(1:100, function(i) {
    Y <- simulateCisEffect.s2s(
      n_i = x[1], p_error = x[2], r = 0, baseline = 3,
      maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
    )
    estimateCisRegEffect.p_error(Y$data)
  }) %>%
    do.call(what = "rbind") %>%
    group_by(model) %>%
    summarise(fp = mean(p_val <= 0.05)) %>%
    mutate(n_i = x[1], p_error = x[2], r_true = 0)
}) %>% do.call(what = "rbind")

g1 <- fp_test.p_error %>% ggplot(aes(x = model, y = fp, color = model)) +
  geom_hline(yintercept = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  labs(title = "False Positive Rate at 5% significance level", y = "False Positive", x = "Model") +
  scale_color_aaas() +
  theme_bw() +
  theme(legend.position = "none")

cowplot::plot_grid(g1, g2, rel_widths = c(2, 3))

## fixed-effect confounders
lognorm_trc_fe <- stan_model(file = "./src/stan_models/lognorm_confounding/trc_multi_cov_fe.stan")
lognorm_trc_fe_null <- stan_model(file = "./src/stan_models/lognorm_lrt/trc_m_fe_null.stan")

estimateCisRegEffect.fe <- function(data) {
  mle_trc <- optimizing(lognorm_trc, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)
  mle_trc_fe <- optimizing(lognorm_trc_fe, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)

  mle_trc_null <- optimizing(lognorm_trc_null, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)
  mle_trc_fe_null <- optimizing(lognorm_trc_fe_null, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)

  chi_sq <- -2 * c(
    mle_trc_null$value - mle_trc$value,
    mle_trc_fe_null$value - mle_trc_fe$value
  )
  df <- c(1, 1)
  p_val <- 1 - pchisq(q = chi_sq, df = df)
  out <- data.frame(
    r_est = sapply(list(mle_trc, mle_trc_fe), function(x) x$par$r),
    model = c("trc", "trc_fe"),
    chi_sq = chi_sq,
    df = df,
    p_val = p_val
  )

  return(out)
}

confounder <- lapply(1:5, function(i){
  rnorm(1000, runif(1, -1, 1), runif(1, 1, 3))
}) %>% do.call(what = "cbind")
confounder.effect <- rep(1, 5)


Y <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3, r = 0.2,
  confounder = confounder, confounder.effect = confounder.effect
)
estimateCisRegEffect.fe(Y$data)
fit_trc_fe <- sampling(lognorm_trc_fe, data = Y$data)

# power
design_power.fe <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(n_n = c(1, 5, 10))
) %>%
  merge(y = data.frame(r = seq(0.01, 1.01, 0.05)))

power_test.fe <- apply(design_power.fe, 1, function(x) {
  lapply(1:200, function(i) {

    confounder <- lapply(1:x[2], function(i){
      rnorm(x[1], runif(1, -1, 1), runif(1, 1, 10))
    }) %>% do.call(what = "cbind")
    confounder.effect <- rep(1, x[2])

    Y <- simulateCisEffect.s2s(
      n_i = x[1], r = x[3],
      confounder = confounder, confounder.effect = confounder.effect,
      baseline = 3, maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
    )

    estimateCisRegEffect.fe(Y$data)
  }) %>%
    do.call(what = "rbind") %>%
    group_by(model) %>%
    summarise(power = mean(p_val <= 0.05)) %>%
    mutate(n_i = x[1], n_n = x[2], r_true = x[3])
}) %>% do.call(what = "rbind")

g2 <- power_test.fe %>%
  ggplot(aes(x = r_true, y = power, color = model, group = model)) +
  geom_point() +
  geom_line() +
  facet_grid(n_i ~ n_n, labeller = label_both) +
  labs(title = "Power at 5% significance level", y = "Power", x = "R (log FC)") +
  scale_color_aaas() +
  theme_bw()

# false positive
design_fp.fe <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(n_n = c(1, 5, 10))
)

fp_test.fe <- apply(design_fp.fe[rep(1:9, each = 100), ], 1, function(x) {
  lapply(1:100, function(i) {

    confounder <- lapply(1:x[2], function(i){
      rnorm(x[1], runif(1, -1, 1), runif(1, 1, 10))
    }) %>% do.call(what = "cbind")
    confounder.effect <- rep(1, x[2])

    Y <- simulateCisEffect.s2s(
      n_i = x[1], r = 0,
      confounder = confounder, confounder.effect = confounder.effect,
      baseline = 3, maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
    )

    estimateCisRegEffect.fe(Y$data)
  }) %>%
    do.call(what = "rbind") %>%
    group_by(model) %>%
    summarise(fp = mean(p_val <= 0.05)) %>%
    mutate(n_i = x[1], n_n = x[2], r_true = 0)
}) %>% do.call(what = "rbind")

g1 <- fp_test.fe %>% ggplot(aes(x = model, y = fp, color = model)) +
  geom_hline(yintercept = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(n_i ~ n_n, labeller = label_both) +
  labs(title = "False Positive Rate at 5% significance level", y = "False Positive", x = "Model") +
  scale_color_aaas() +
  theme_bw() +
  theme(legend.position = "none")


cowplot::plot_grid(g1, g2, rel_widths = c(2, 3))
