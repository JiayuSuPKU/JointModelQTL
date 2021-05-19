library(tidyverse)
library(ggsci)
library(ggrepel)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan_models <- load_models()

## basic

Y <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3, r = 0.5
)

# runtime comparisons
system.time({
  optimizing(stan_models$basic$joint, data = Y$data)
})
system.time({
  sampling(stan_models$basic$joint, data = Y$data, refresh = 0)
})

# point estimator
results <- estimateCisRegEffects(
  data = Y$data, stan_models = stan_models,
  model = "joint", method = "optimizing",
  significance_test = F
)
# with likelihood ratio test
results <- estimateCisRegEffects(
  data = Y$data, stan_models = stan_models,
  model = "joint", method = "optimizing",
  significance_test = T
)

# helper functions
testPower <- function(design_power, stan_models, par = "baseline",
                      weighted_ase = FALSE, phasing_error = "null",
                      confounders = "null") {
  out <- apply(design_power, 1, function(x) {
    lapply(1:200, function(i) {
      sim_pars <- list(
        n_i = x[1], r = x[3],
        maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
      )
      sim_pars$baseline <- ifelse(par == "baseline", x[2], 3)
      sim_pars$p_error <- ifelse(par == "p_error", x[2], NULL)

      Y <- do.call(simulateCisEffect.s2s, sim_pars)

      joint <- estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = "joint", method = "optimizing",
        significance_test = T, confounders = confounders,
        weighted_ase = weighted_ase, phasing_error = phasing_error
      )

      trc <- estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = "trc", method = "optimizing",
        significance_test = T, confounders = confounders,
        weighted_ase = weighted_ase, phasing_error = phasing_error
      )

      ase <- estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = "ase", method = "optimizing",
        significance_test = T, confounders = confounders,
        weighted_ase = weighted_ase, phasing_error = phasing_error
      )

      return(rbind(joint, trc, ase) %>%
        mutate(n_i = x[1], !!par := x[2], r_true = x[3]))
    }) %>% do.call(what = "rbind")
  }) %>% do.call(what = "rbind")

  return(out)
}

testFalsePositive <- function(design_fp, stan_models, par = "baseline",
                              weighted_ase = FALSE, phasing_error = "null",
                              confounders = "null") {
  out <- apply(design_fp[rep(1:9, each = 100), ], 1, function(x) {
    sim_pars <- list(
      n_i = x[1], r = 0,
      maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
    )
    sim_pars$baseline <- ifelse(par == "baseline", x[2], 3)
    sim_pars$p_error <- ifelse(par == "p_error", x[2], NULL)

    lapply(1:100, function(i) {
      Y <- do.call(simulateCisEffect.s2s, sim_pars)

      joint <- estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = "joint", method = "optimizing",
        significance_test = T, confounders = confounders,
        weighted_ase = weighted_ase, phasing_error = phasing_error
      )

      trc <- estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = "trc", method = "optimizing",
        significance_test = T, confounders = confounders,
        weighted_ase = weighted_ase, phasing_error = phasing_error
      )

      ase <- estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = "ase", method = "optimizing",
        significance_test = T, confounders = confounders,
        weighted_ase = weighted_ase, phasing_error = phasing_error
      )

      return(rbind(joint, trc, ase))
    }) %>%
      do.call(what = "rbind") %>%
      group_by(model) %>%
      summarise(fp = mean(p_val <= 0.05)) %>%
      mutate(n_i = x[1], !!par := x[2], r_true = 0)
  }) %>% do.call(what = "rbind")

  return(out)
}

# power analysis
design_power <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(b = c(1, 3, 5))
) %>%
  merge(y = data.frame(r = seq(0.01, 1.01, 0.05)))

test_power <- testPower(design_power = design_power, stan_models = stan_models)
df_power <- test_power %>%
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

test_fp <- testFalsePositive(design_fp = design_fp, stan_models = stan_models)
test_fp %>% ggplot(aes(x = model, y = fp, color = model)) +
  geom_hline(yintercept = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  labs(title = "False Positive Rate at 5% significance level", y = "False Positive", x = "Model") +
  scale_color_npg() +
  theme_bw()


## weighted ase
# power
test_power.w <- testPower(design_power = design_power, stan_models = stan_models, weighted_ase = T)
df_power.w <- test_power.w %>%
  group_by(n_i, baseline, r_true, model) %>%
  mutate(filter = p_val <= 0.05, power = sum(filter) / 200)
df_power.w %>%
  ggplot(aes(x = r_true, y = power, color = model, group = model)) +
  geom_point() +
  geom_line() +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  labs(title = "Power at 5% significance level", y = "Power", x = "R (log FC)") +
  scale_color_aaas() +
  theme_bw()

# irregular case
set.seed(202105)
Y <- simulateCisEffect.s2s(
  n_i = 500, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3, r = 1
)
estimateCisRegEffects(Y$data, stan_models = stan_models, weighted_ase = T)

mle_joint <- optimizing(stan_models$basic$joint, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
mle_joint_w <- optimizing(stan_models$basic$joint_weighted, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)

fit_joint <- sampling(stan_models$basic$joint, data = Y$data)
fit_joint_w <- sampling(stan_models$basic$joint_weighted, data = Y$data)
data.frame(
  r_est = c(extract(fit_joint, pars = "r")[[1]], extract(fit_joint_w, pars = "r")[[1]]),
  model = rep(c("ase", "weighted ase"), each = 4000)
) %>%
  ggplot(aes(x = r_est, group = model, color = model)) +
  geom_density() +
  geom_vline(xintercept = 1) +
  scale_color_aaas() +
  labs(title = "Posterior distribution of R (n_i = 500, b = 3)") +
  theme_bw()


## phasing error
Y <- simulateCisEffect.s2s(
  n_i = 500, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3, r = 1, p_error = 0.3
)

# the joint model does not work very well
# it turns out that the problem lies in the calculation of likelihood
rbind(
  estimateCisRegEffects(Y$data, stan_models, model = "joint", phasing_error = "null"),
  estimateCisRegEffects(Y$data, stan_models, model = "joint", phasing_error = "known"),
  estimateCisRegEffects(Y$data, stan_models, model = "joint", phasing_error = "unknown")
) %>% mutate(phasing_error = c("null", "known", "unknown"), r_true = 1) %>%
  select(r_est, r_true, model, phasing_error, lr, df, p_val)

rbind(
  estimateCisRegEffects(
    Y$data, stan_models, model = "joint", method = "sampling", phasing_error = "null",
    significance_test = T, return_posterior = "mean"
  ),
  estimateCisRegEffects(
    Y$data, stan_models, model = "joint", method = "sampling", phasing_error = "known",
    significance_test = T, return_posterior = "mean"
  ),
  estimateCisRegEffects(
    Y$data, stan_models, model = "joint", method = "sampling", phasing_error = "unknown",
    significance_test = T, return_posterior = "mean"
  )
) %>% mutate(phasing_error = c("null", "known", "unknown"), r_true = 1) %>%
  select(r_est, r_true, model, phasing_error, bf)


rbind(
  estimateCisRegEffects(Y$data, stan_models, model = "ase", phasing_error = "null"),
  estimateCisRegEffects(Y$data, stan_models, model = "ase", phasing_error = "known"),
  estimateCisRegEffects(Y$data, stan_models, model = "ase", phasing_error = "unknown"),
  estimateCisRegEffects(Y$data, stan_models, model = "trc", phasing_error = "null"),
  estimateCisRegEffects(Y$data, stan_models, model = "trc", phasing_error = "known"),
  estimateCisRegEffects(Y$data, stan_models, model = "trc", phasing_error = "unknown")
) %>% mutate(phasing_error = rep(c("null", "known", "unknown"), 2), r_true = 1) %>%
  select(r_est, r_true, model, phasing_error, lr, df, p_val)

rbind(
  estimateCisRegEffects(
    Y$data, stan_models, model = "ase", method = "sampling", phasing_error = "null",
    significance_test = T, return_posterior = "mean"
  ),
  estimateCisRegEffects(
    Y$data, stan_models, model = "ase", method = "sampling", phasing_error = "known",
    significance_test = T, return_posterior = "mean"
  ),
  estimateCisRegEffects(
    Y$data, stan_models, model = "ase", method = "sampling", phasing_error = "unknown",
    significance_test = T, return_posterior = "mean"
  )
) %>% mutate(phasing_error = c("null", "known", "unknown"), r_true = 1) %>%
  select(r_est, r_true, model, phasing_error, bf)


# power
design_power.p_error <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(p_error = c(0.05, 0.15, 0.25))
) %>%
  merge(y = data.frame(r = seq(0.01, 1.01, 0.05)))

test_power.p_error_null <- testPower(
  design_power.p_error, stan_models,
  par = "p_error", phasing_error = "null"
)
test_power.p_error_k <- testPower(
  design_power.p_error, stan_models,
  par = "p_error", phasing_error = "known"
)
test_power.p_error_u <- testPower(
  design_power.p_error, stan_models,
  par = "p_error", phasing_error = "unknown"
)

g2 <- test_power.p_error_k %>%
  group_by(n_i, p_error, r_true, model) %>%
  mutate(filter = p_val <= 0.05, power = sum(filter) / 200) %>%
  ggplot(aes(x = r_true, y = power, color = model, group = model)) +
  geom_point() +
  geom_line() +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  labs(title = "Mixture Models (known): Power at 5% significance level", y = "Power", x = "R (log FC)") +
  scale_color_aaas() +
  theme_bw()

g3 <- test_power.p_error_u %>%
  group_by(n_i, p_error, r_true, model) %>%
  mutate(filter = p_val <= 0.05, power = sum(filter) / 200) %>%
  ggplot(aes(x = r_true, y = power, color = model, group = model)) +
  geom_point() +
  geom_line() +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  labs(title = "Mixture Models (unknown): Power at 5% significance level", y = "Power", x = "R (log FC)") +
  scale_color_aaas() +
  theme_bw()

cowplot::plot_grid(g2, g3)


# false positive
design_fp.p_error <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(p_error = c(0.05, 0.15, 0.25))
)


test_fp.p_error_null <- testFalsePositive(
  design_fp = design_fp.p_error, stan_models = stan_models,
  par = "p_error", phasing_error = "null"
)
test_fp.p_error_k <- testFalsePositive(
  design_fp = design_fp.p_error, stan_models = stan_models,
  par = "p_error", phasing_error = "known"
)
test_fp.p_error_u <- testFalsePositive(
  design_fp = design_fp.p_error, stan_models = stan_models,
  par = "p_error", phasing_error = "unknown"
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

g1 <- test_fp.p_error_null %>% ggplot(aes(x = model, y = fp, color = model)) +
  geom_hline(yintercept = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  labs(title = "Basic Models: FP at 5%", y = "False Positive", x = "Model") +
  scale_color_aaas() +
  theme_bw() +
  theme(legend.position = "none")

g2 <- test_fp.p_error_k %>% ggplot(aes(x = model, y = fp, color = model)) +
  geom_hline(yintercept = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  labs(title = "Mixture Models (known): FP at 5%", y = "False Positive", x = "Model") +
  scale_color_aaas() +
  theme_bw() +
  theme(legend.position = "none")

g3 <- test_fp.p_error_u %>% ggplot(aes(x = model, y = fp, color = model)) +
  geom_hline(yintercept = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  labs(title = "Mixture Models (unknown): FP at 5%", y = "False Positive", x = "Model") +
  scale_color_aaas() +
  theme_bw() +
  theme(legend.position = "none")

cowplot::plot_grid(g1, g2, g3, nrow = 1)


## fixed-effect confounders
confounder <- lapply(1:5, function(i) {
  rnorm(1000, runif(1, -1, 1), runif(1, 1, 3))
}) %>% do.call(what = "cbind")
confounder.effect <- rep(1, 5)


Y <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3, r = 0.5,
  confounder = confounder, confounder.effect = confounder.effect
)
estimateCisRegEffects(
  Y$data,
  stan_models = stan_models, model = "trc", confounders = "fixed"
)
estimateCisRegEffects(
  Y$data,
  stan_models = stan_models, model = "trc", confounders = "null"
)

fit_trc_fe <- sampling(lognorm_trc_fe, data = Y$data)

# power
design_power.fe <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(n_n = c(1, 5, 10))
) %>%
  merge(y = data.frame(r = seq(0.01, 1.01, 0.05)))

power_test.fe <- apply(design_power.fe, 1, function(x) {
  lapply(1:200, function(i) {
    confounder <- lapply(1:x[2], function(i) {
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
    confounder <- lapply(1:x[2], function(i) {
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

## random-effect confounders
lognorm_trc_re <- stan_model(file = "./src/stan_models/lognorm_confounding/trc_indiv_re.stan")
lognorm_trc_re_null <- stan_model(file = "./src/stan_models/lognorm_lrt/trc_indiv_re_null.stan")

estimateCisRegEffect.re <- function(data) {
  mle_trc <- optimizing(lognorm_trc, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)
  mle_trc_re <- optimizing(lognorm_trc_re, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)

  mle_trc_null <- optimizing(lognorm_trc_null, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)
  mle_trc_re_null <- optimizing(lognorm_trc_re_null, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)

  chi_sq <- -2 * c(
    mle_trc_null$value - mle_trc$value,
    mle_trc_re_null$value - mle_trc_re$value
  )
  df <- c(1, 1)
  p_val <- 1 - pchisq(q = chi_sq, df = df)
  out <- data.frame(
    r_est = sapply(list(mle_trc, mle_trc_re), function(x) x$par$r),
    model = c("trc", "trc_re"),
    chi_sq = chi_sq,
    df = df,
    p_val = p_val
  )

  return(out)
}

Y <- simulateCisEffect.s2s(
  n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3, r = 1,
  origin = sample(rep(1:100, each = 10), 1000, replace = F), origin.effect = rnorm(100)
)

estimateCisRegEffect.re(Y$data)

# power
design_power.re <- data.frame(
  n_i = c(100, 400, 400, 1000, 1000, 1000),
  n_indiv = c(100, 100, 200, 100, 200, 500)
) %>%
  merge(y = data.frame(r = seq(0.01, 1.01, 0.05)))

power_test.re <- apply(design_power.re, 1, function(x) {
  lapply(1:200, function(i) {
    origin <- sample(rep(1:x[2], each = x[1] / x[2]), x[1], replace = F)
    origin.effect <- rnorm(x[2])

    Y <- simulateCisEffect.s2s(
      n_i = x[1], r = x[3],
      origin = origin, origin.effect = origin.effect,
      baseline = 3, maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
    )

    estimateCisRegEffect.re(Y$data)
  }) %>%
    do.call(what = "rbind") %>%
    group_by(model) %>%
    summarise(power = mean(p_val <= 0.05)) %>%
    mutate(n_i = x[1], n_indiv = x[2], r_true = x[3])
}) %>% do.call(what = "rbind")

g2 <- power_test.re %>%
  ggplot(aes(x = r_true, y = power, color = model, group = model)) +
  geom_point() +
  geom_line() +
  facet_grid(n_i ~ n_indiv, labeller = label_both) +
  labs(title = "Power at 5% significance level", y = "Power", x = "R (log FC)") +
  scale_color_aaas() +
  theme_bw()

# false positive
design_fp.re <- data.frame(
  n_i = c(100, 400, 400, 1000, 1000, 1000),
  n_indiv = c(100, 100, 200, 100, 200, 500)
)

fp_test.re <- apply(design_fp.re[rep(1:6, each = 100), ], 1, function(x) {
  lapply(1:100, function(i) {
    origin <- sample(rep(1:x[2], each = x[1] / x[2]), x[1], replace = F)
    origin.effect <- rnorm(x[2])

    Y <- simulateCisEffect.s2s(
      n_i = x[1], r = 0,
      origin = origin, origin.effect = origin.effect,
      baseline = 3, maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
    )

    estimateCisRegEffect.re(Y$data)
  }) %>%
    do.call(what = "rbind") %>%
    group_by(model) %>%
    summarise(fp = mean(p_val <= 0.05)) %>%
    mutate(n_i = x[1], n_indiv = x[2], r_true = 0)
}) %>% do.call(what = "rbind")

g1 <- fp_test.re %>% ggplot(aes(x = model, y = fp, color = model)) +
  geom_hline(yintercept = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(n_i ~ n_indiv, labeller = label_both) +
  labs(title = "False Positive Rate at 5% significance level", y = "False Positive", x = "Model") +
  scale_color_aaas() +
  theme_bw() +
  theme(legend.position = "none")

cowplot::plot_grid(g1, g2, rel_widths = c(2, 3))



