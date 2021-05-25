library(tidyverse)
library(ggsci)
library(ggrepel)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load("./scripts/benchmark/basic.RData")
load("./scripts/benchmark/linear.RData")

stan_models <- load_models()

pal_col <- c(
  "ase" = pal_aaas()(10)[1],
  "joint" = pal_aaas()(10)[2],
  "trc" = pal_aaas()(10)[3],
  "joint_linear" = pal_npg()(10)[2],
  "trc_linear" = pal_npg()(10)[5],
  "trc_lm" = pal_aaas()(10)[9]
)

## 1) Estimation accuracy

# sampling
r_est_ac.s <- data.frame()
for (r in c(0, 0.3, 0.6, 0.9)) {
  Y <- simulateCisEffect.s2s(
    n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
    phi = 3, theta = 30, baseline = 3, r = r
  )

  trc <- estimateCisRegEffects(
    data = Y$data, stan_models = stan_models,
    model = "trc", method = "sampling", significance_test = F
  )
  joint <- estimateCisRegEffects(
    data = Y$data, stan_models = stan_models,
    model = "joint", method = "sampling", significance_test = F
  )
  trc_linear <- estimateCisRegEffects(
    data = Y$data, stan_models = stan_models,
    model = "trc", method = "sampling", significance_test = F, linear_trc = T
  )
  joint_linear <- estimateCisRegEffects(
    data = Y$data, stan_models = stan_models,
    model = "joint", method = "sampling", significance_test = F, linear_trc = T
  )

  r_est <- data.frame(
    r_est = c(trc$r_est, joint$r_est, trc_linear$r_est, joint_linear$r_est),
    model = rep(c("trc", "joint", "trc_linear", "joint_linear"), each = 4000),
    r_true = r
  )
  r_est_ac.s <- rbind(r_est_ac.s, r_est)
}

g1 <- r_est_ac.s %>%
  mutate(model = factor(model, c("trc", "joint", "trc_linear", "joint_linear"))) %>%
  ggplot(aes(x = r_est, group = model, color = model)) +
  facet_wrap(~r_true, labeller = label_both, ncol = 1) +
  geom_vline(aes(xintercept = r_true), color = "black") +
  geom_density(size = 1) +
  labs(x = "R_est", y = "Posterior density", title = "Accuracy: Posterior") +
  scale_color_manual(values = pal_col) +
  theme_bw()

# optimizing
r_est_ac.o <- data.frame()
for (r in c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) {
  out <- lapply(1:200, function(i) {
    Y <- simulateCisEffect.s2s(
      n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = r
    )

    trc <- estimateCisRegEffects(
      data = Y$data, stan_models = stan_models,
      model = "trc", method = "optimizing", significance_test = F
    )
    joint <- estimateCisRegEffects(
      data = Y$data, stan_models = stan_models,
      model = "joint", method = "optimizing", significance_test = F
    )
    trc_linear <- estimateCisRegEffects(
      data = Y$data, stan_models = stan_models,
      model = "trc", method = "optimizing", significance_test = F, linear_trc = T
    )
    joint_linear <- estimateCisRegEffects(
      data = Y$data, stan_models = stan_models,
      model = "joint", method = "optimizing", significance_test = F, linear_trc = T
    )
    lm.model <- lm(log1p_T ~ G, data = data.frame(log1p_T = Y$data$log1p_T, G = Y$data$G))

    data.frame(
      r_est = c(trc$r_est, joint$r_est, trc_linear$r_est, joint_linear$r_est, lm.model$coefficients["G"] * 2),
      model = c("trc", "joint", "trc_linear", "joint_linear", "trc_lm")
    )
  }) %>%
    do.call(what = "rbind") %>%
    mutate(r_true = r)

  r_est_ac.o <- rbind(r_est_ac.o, out)
}

g2 <- r_est_ac.o %>%
  mutate(model = factor(model, c("trc", "joint", "trc_linear", "joint_linear", "trc_lm"))) %>%
  ggplot(aes(x = r_true, y = r_est, group = interaction(model, r_true), fill = model)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(y = r_true), color = "black", shape = 17, size = 3) +
  labs(x = "R_true", y = "R_est", title = "Accuracy: MLE") +
  ylim(c(-0.5, 1.5)) +
  scale_fill_manual(values = pal_col) +
  theme_bw()

cowplot::plot_grid(g1, g2)


## 3) Run time

# sampling
runtime.s <- data.frame()
for (n_i in c(100, 400, 800, 1600, 3200, 6400)) {
  runtime <- lapply(1:3, function(i) {
    Y <- simulateCisEffect.s2s(
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5
    )

    runtime.trc <- system.time({
      estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = "trc", method = "sampling", significance_test = F
      )
    })
    runtime.joint <- system.time({
      estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = "joint", method = "sampling", significance_test = F
      )
    })
    runtime.trc_linear <- system.time({
      estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = "trc", method = "sampling", significance_test = F, linear_trc = T
      )
    })
    runtime.joint_linear <- system.time({
      estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = "joint", method = "sampling", significance_test = F, linear_trc = T
      )
    })

    rbind(runtime.trc, runtime.joint, runtime.trc_linear, runtime.joint_linear) %>%
      data.frame() %>%
      mutate(model = c("trc", "joint", "trc_linear", "joint_linear"), n_i = n_i)
  }) %>% do.call(what = "rbind")

  runtime.s <- rbind(runtime.s, runtime)
}

g1 <- runtime.s %>%
  ggplot(aes(x = n_i, y = elapsed, group = model, color = model)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, se = F) +
  scale_color_manual(values = pal_col) +
  ylim(c(0, 50)) +
  labs(
    x = "Number of samples", y = "Total Run Time (seconds)",
    color = "Model", title = "Runtime: Posterior"
  ) +
  theme_bw()


# optimizing
runtime.o <- data.frame()
for (n_i in c(100, 400, 800, 1600, 3200, 6400)) {
  runtime <- lapply(1:100, function(i) {
    Y <- simulateCisEffect.s2s(
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5
    )

    runtime.trc <- system.time({
      estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = "trc", method = "optimizing", significance_test = F
      )
    })
    runtime.joint <- system.time({
      estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = "joint", method = "optimizing", significance_test = F
      )
    })
    runtime.trc_linear <- system.time({
      estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = "trc", method = "optimizing", significance_test = F, linear_trc = T
      )
    })
    runtime.joint_linear <- system.time({
      estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = "joint", method = "optimizing", significance_test = F, linear_trc = T
      )
    })
    runtime.lm <- system.time({
      lm(log1p_T ~ G, data = data.frame(log1p_T = Y$data$log1p_T, G = Y$data$G))
    })

    rbind(runtime.trc, runtime.joint, runtime.trc_linear, runtime.joint_linear, runtime.lm) %>%
      data.frame() %>%
      mutate(model = c("trc", "joint", "trc_linear", "joint_linear", "trc_lm"), n_i = n_i)
  }) %>% do.call(what = "rbind")

  runtime.o <- rbind(runtime.o, runtime)
}

g2 <- runtime.o %>%
  ggplot(aes(x = n_i, y = elapsed, group = model, color = model)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, se = F) +
  scale_color_manual(values = pal_col) +
  # ylim(c(0, 0.3)) +
  labs(
    x = "Number of samples", y = "Total Run Time (seconds)",
    color = "Model", title = "Runtime: MLE"
  ) +
  theme_bw()

cowplot::plot_grid(g1, g2)


## 4) False positivity
design_fp <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(b = c(1, 3, 5))
)

# lm
fp_lm <- apply(design_fp[rep(1:9, each = 50), ], 1, function(x) {
  lapply(1:100, function(i) {
    sim_pars <- list(
      n_i = x["n_i"], r = 0, baseline = x["b"],
      maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
    )

    Y <- do.call(simulateCisEffect.s2s, sim_pars)
    lm.data <- data.frame(log1p_T = Y$data$log1p_T, G = Y$data$G)
    lm.model.0 <- lm(log1p_T ~ 1, data = lm.data)
    lm.model.1 <- lm(log1p_T ~ G, data = lm.data)

    data.frame(
      r_prime_est = summary(lm.model)$coef["G", "Estimate"],
      model = "trc_lm"
    ) %>% mutate(
      lr = -2 * (logLik(lm.model.0) - logLik(lm.model.1)),
      df = 1,
      p_val = 1 - pchisq(q = lr, df = df)
    )
  }) %>%
    do.call(what = "rbind") %>%
    mutate(n_i = x["n_i"], baseline = x["b"], r_true = 0)
}) %>%
  do.call(what = "rbind")
df_fp_lm <- fp_lm %>%
  mutate(rep_idx = rep(rep(1:50, each = 100), 9)) %>%
  group_by(n_i, baseline, model, rep_idx) %>%
  summarise(fp = mean(p_val <= 0.05))

# linear stan models
fp_linear.o <- testFalsePositive(
  design_fp = design_fp, stan_models = stan_models, n_rep = 50,
  models = c("trc", "joint"), linear_trc = T
)
df_fp_linear.o <- fp_linear.o %>%
  mutate(model = paste0(model, "_linear")) %>%
  mutate(rep_idx = rep(rep(1:50, each = 200), 9)) %>%
  group_by(n_i, baseline, model, rep_idx) %>%
  summarise(fp = mean(p_val <= 0.05))

g1 <- rbind(df_fp.o, df_fp_linear.o, df_fp_lm) %>%
  mutate(model = factor(model, levels = c("ase", "trc", "joint", "trc_linear", "joint_linear", "trc_lm"))) %>%
  ggplot(aes(x = model, y = fp, fill = model)) +
  geom_hline(yintercept = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  labs(title = "False Positivity: MLE + linear", y = "False Positive", x = "Model") +
  scale_fill_manual(values = pal_col) +
  theme_bw()


## 5) Power
design_power <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(b = c(1, 3, 5))
) %>% merge(y = data.frame(r = seq(0.01, 1.01, 0.05)))

# lm
power_lm <- apply(design_power, 1, function(x) {
  lapply(1:100, function(i) {
    sim_pars <- list(
      n_i = x["n_i"], r = x["r"], baseline = x["b"],
      maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
    )

    Y <- do.call(simulateCisEffect.s2s, sim_pars)
    lm.data <- data.frame(log1p_T = Y$data$log1p_T, G = Y$data$G)
    lm.model.0 <- lm(log1p_T ~ 1, data = lm.data)
    lm.model.1 <- lm(log1p_T ~ G, data = lm.data)

    data.frame(
      r_prime_est = summary(lm.model)$coef["G", "Estimate"],
      model = "trc_lm"
    ) %>% mutate(
      lr = -2 * (logLik(lm.model.0) - logLik(lm.model.1)),
      df = 1,
      p_val = 1 - pchisq(q = lr, df = df)
    )
  }) %>%
    do.call(what = "rbind") %>%
    mutate(n_i = x["n_i"], baseline = x["b"], r_true = x["r"])
}) %>% do.call(what = "rbind")
df_power_lm <- power_lm %>%
  group_by(n_i, baseline, model, r_true) %>%
  summarise(power = mean(p_val <= 0.05))

# linear stan models
power_linear.o <- testPower(
  design_power = design_power, stan_models = stan_models,
  models = c("trc", "joint"), linear_trc = T
)
df_power_linear.o <- power_linear.o %>%
  mutate(model = paste0(model, "_linear")) %>%
  group_by(n_i, baseline, r_true, model) %>%
  summarise(power = mean(p_val <= 0.05))

g2 <- rbind(df_power.o, df_power_linear.o, df_power_lm) %>%
  mutate(model = factor(model, levels = c("ase", "trc", "joint", "trc_linear", "joint_linear", "trc_lm"))) %>%
  ggplot(aes(x = r_true, y = power, color = model, group = model)) +
  geom_point() +
  geom_line() +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  labs(title = "Power: MLE + linear", y = "Power", x = "R (log FC)") +
  scale_color_manual(values = pal_col) +
  theme_bw()

cowplot::plot_grid(g1, g2)

## 6) ROC

# optimizing
roc <- rbind(
  power.o,
  power_linear.o %>% mutate(model = paste0(model, "_linear")),
  power_lm %>% rename(r_est = r_prime_est) %>% mutate(r_est = r_est * 2),
  fp.o %>% group_by(n_i, baseline, model) %>% sample_n(2100),
  fp_linear.o %>% mutate(model = paste0(model, "_linear")) %>% group_by(n_i, baseline, model) %>% sample_n(2100),
  fp_lm %>% rename(r_est = r_prime_est) %>% mutate(r_est = r_est * 2) %>% group_by(n_i, baseline, model) %>% sample_n(2100)
) %>%
  mutate(is_negative = r_true == 0)

df_roc <- lapply(seq(0, 1.01, 0.01), function(threshold){
  roc %>%
    mutate(is_pred_positive = p_val < threshold) %>%
    group_by(n_i, baseline, model) %>%
    summarise(
      fp = mean(is_pred_positive[is_negative]),
      tp = mean(is_pred_positive[!is_negative])
    ) %>%
    mutate(p_thresh= threshold)
}) %>% do.call(what = "rbind")

df_roc %>%
  mutate(model = factor(model, levels = c("ase", "trc", "joint", "trc_linear", "joint_linear", "trc_lm"))) %>%
  ggplot(aes(x = fp, y = tp, group = model, color = model)) +
  geom_line() +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  scale_color_manual(values = pal_col) +
  labs(title = "ROC curve: MLE + linear", y = "True Positive Rate", x = "False Positive Rate") +
  theme_bw()




