library(tidyverse)
library(ggsci)
library(ggrepel)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load("./scripts/benchmark/basic.RData")
load("./scripts/benchmark/phasing_error.RData")

stan_models <- load_models()

pal_col <- c(
  "ase" = pal_aaas()(10)[1],
  "joint" = pal_aaas()(10)[2],
  "trc" = pal_aaas()(10)[3],
  "ase_mix_k" = pal_npg()(10)[2],
  "joint_mix_k" = pal_npg()(10)[5],
  "ase_mix_u" = pal_npg()(10)[3],
  "joint_mix_u" = pal_npg()(10)[9]
)


## 2) Sensitivity (on sample size and baseline expression)

# sampling
r_est_sen.s <- data.frame()
for (p_error in c(0.1, 0.2, 0.3)) {
  for (n_i in c(100, 300, 500)) {
    r_est_null <- sensitivityR(
      models = c("ase", "joint"),
      method = "sampling", return_posterior = "full",
      phasing_error = "null",
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5, p_error = p_error
    ) %>% mutate(mixture = "null", p_error = p_error, n_i = n_i)

    r_est_known <- sensitivityR(
      models = c("ase", "joint"),
      method = "sampling", return_posterior = "full",
      phasing_error = "known",
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5, p_error = p_error
    ) %>% mutate(mixture = "known", p_error = p_error, n_i = n_i)

    r_est_unknown <- sensitivityR(
      models = c("ase", "joint"),
      method = "sampling", return_posterior = "full",
      phasing_error = "unknown",
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5, p_error = p_error
    ) %>% mutate(mixture = "unknown", p_error = p_error, n_i = n_i)

    r_est_sen.s <- rbind(r_est_sen.s, r_est_null, r_est_known, r_est_unknown)
  }
}

g1 <- r_est_sen.s %>%
  filter(mixture == "null") %>%
  mutate(model = factor(model, c("ase", "joint"))) %>%
  ggplot(aes(x = r_est, group = model, color = model)) +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  geom_vline(xintercept = 0.5, color = "black") +
  geom_density(size = 1) +
  ylim(c(0, 6)) +
  labs(x = "R_est", y = "Posterior density", title = "Sensitivity: Posterior (basic)") +
  scale_color_manual(values = pal_col) +
  theme_bw()

g2 <- r_est_sen.s %>%
  filter(mixture == "known") %>%
  mutate(model = paste0(model, "_mix_k")) %>%
  mutate(model = factor(model, c("ase_mix_k", "joint_mix_k"))) %>%
  ggplot(aes(x = r_est, group = model, color = model)) +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  geom_vline(xintercept = 0.5, color = "black") +
  geom_density(size = 1) +
  ylim(c(0, 6)) +
  labs(x = "R_est", y = "Posterior density", title = "(Mixture with known phasing error)") +
  scale_color_manual(values = pal_col) +
  theme_bw()

g3 <- r_est_sen.s %>%
  filter(mixture == "unknown") %>%
  mutate(model = paste0(model, "_mix_u")) %>%
  mutate(model = factor(model, c("ase_mix_u", "joint_mix_u"))) %>%
  ggplot(aes(x = r_est, group = model, color = model)) +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  geom_vline(xintercept = 0.5, color = "black") +
  geom_density(size = 1) +
  ylim(c(0, 6)) +
  labs(x = "R_est", y = "Posterior density", title = "(Mixture with unknown phasing error)") +
  scale_color_manual(values = pal_col) +
  theme_bw()

cowplot::plot_grid(g1, g2, g3, nrow = 1)

# optimizing
r_est_sen.o <- data.frame()
for (p_error in c(0.1, 0.2, 0.3)) {
  for (n_i in c(100, 300, 500)) {
    r_est_null <- sensitivityR(
      models = c("ase", "joint"),
      method = "optimizing", phasing_error = "null",
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5, p_error = p_error
    ) %>% mutate(mixture = "null", p_error = p_error, n_i = n_i)

    r_est_known <- sensitivityR(
      models = c("ase", "joint"),
      method = "optimizing", phasing_error = "known",
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5, p_error = p_error
    ) %>% mutate(mixture = "known", p_error = p_error, n_i = n_i)

    r_est_unknown <- sensitivityR(
      models = c("ase", "joint"),
      method = "optimizing", phasing_error = "unknown",
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5, p_error = p_error
    ) %>% mutate(mixture = "unknown", p_error = p_error, n_i = n_i)

    r_est_sen.o <- rbind(r_est_sen.o, r_est_null, r_est_known, r_est_unknown)
  }
}

g1 <- r_est_sen.o %>%
  filter(mixture == "null") %>%
  mutate(model = factor(model, c("ase", "trc", "joint"))) %>%
  ggplot(aes(x = model, y = r_est, group = model, fill = model)) +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0.5) +
  labs(x = "R_true", y = "R_est", title = "Sensitivity: MLE (basic)") +
  ylim(c(-0.5, 1)) +
  scale_fill_manual(values = pal_col) +
  theme_bw() +
  theme(axis.text.x = element_blank())

g2 <- r_est_sen.o %>%
  filter(mixture == "known") %>%
  mutate(model = paste0(model, "_mix_k")) %>%
  mutate(model = factor(model, c("ase_mix_k", "joint_mix_k"))) %>%
  ggplot(aes(x = model, y = r_est, group = model, fill = model)) +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0.5) +
  labs(x = "R_true", y = "R_est", title = "(Mixture with known phasing error)") +
  ylim(c(-0.5, 1)) +
  scale_fill_manual(values = pal_col) +
  theme_bw() +
  theme(axis.text.x = element_blank())

g3 <- r_est_sen.o %>%
  filter(mixture == "unknown") %>%
  mutate(model = paste0(model, "_mix_u")) %>%
  mutate(model = factor(model, c("ase_mix_u", "joint_mix_u"))) %>%
  ggplot(aes(x = model, y = r_est, group = model, fill = model)) +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0.5) +
  labs(x = "R_true", y = "R_est", title = "(Mixture with unknown phasing error)") +
  ylim(c(-0.5, 1)) +
  scale_fill_manual(values = pal_col) +
  theme_bw() +
  theme(axis.text.x = element_blank())

cowplot::plot_grid(g1, g2, g3, nrow = 1)


## 3) Run time

# sampling
runtime.s <- data.frame()
for (n_i in c(100, 400, 800, 1600, 3200, 6400)) {
  runtime <- lapply(1:3, function(i) {
    runtime_null <- getRunTime(
      models = c("ase", "joint"),
      method = "sampling", significance_test = F, phasing_error = "null",
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5, p_error = 0.2
    )
    runtime_k <- getRunTime(
      models = c("ase", "joint"),
      method = "sampling", significance_test = F, phasing_error = "known",
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5, p_error = 0.2
    ) %>% mutate(model = paste0(model, "_mix_k"))
    runtime_u <- getRunTime(
      models = c("ase", "joint"),
      method = "sampling", significance_test = F, phasing_error = "unknown",
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5, p_error = 0.2
    ) %>% mutate(model = paste0(model, "_mix_u"))

    rbind(runtime_null, runtime_k, runtime_u) %>% mutate(n_i = n_i)
  }) %>% do.call(what = "rbind")

  runtime.s <- rbind(runtime.s, runtime)
}

g1 <- runtime.s %>%
  mutate(model = factor(model, levels = c("ase", "joint", "ase_mix_k", "joint_mix_k", "ase_mix_u", "joint_mix_u"))) %>%
  ggplot(aes(x = n_i, y = elapsed, group = model, color = model)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, se = F) +
  scale_color_manual(values = pal_col) +
  ylim(c(0, 100)) +
  labs(
    x = "Number of samples", y = "Total Run Time (seconds)",
    color = "Model", title = "Runtime: Posterior + phasing error"
  ) +
  theme_bw()


# optimizing
runtime.o <- data.frame()
for (n_i in c(100, 400, 800, 1600, 3200, 6400)) {
  runtime <- lapply(1:100, function(i) {
    runtime_null <- getRunTime(
      models = c("ase", "joint"),
      method = "optimizing", significance_test = F, phasing_error = "null",
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5, p_error = 0.2
    )
    runtime_k <- getRunTime(
      models = c("ase", "joint"),
      method = "optimizing", significance_test = F, phasing_error = "known",
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5, p_error = 0.2
    ) %>% mutate(model = paste0(model, "_mix_k"))
    runtime_u <- getRunTime(
      models = c("ase", "joint"),
      method = "optimizing", significance_test = F, phasing_error = "unknown",
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5, p_error = 0.2
    ) %>% mutate(model = paste0(model, "_mix_u"))

    rbind(runtime_null, runtime_k, runtime_u) %>% mutate(n_i = n_i)
  }) %>% do.call(what = "rbind")

  runtime.o <- rbind(runtime.o, runtime)
}

g2 <- runtime.o %>%
  mutate(model = factor(model, levels = c("ase", "joint", "ase_mix_k", "joint_mix_k", "ase_mix_u", "joint_mix_u"))) %>%
  ggplot(aes(x = n_i, y = elapsed, group = model, color = model)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, se = F) +
  scale_color_manual(values = pal_col) +
  ylim(c(0, 0.4)) +
  labs(
    x = "Number of samples", y = "Total Run Time (seconds)",
    color = "Model", title = "Runtime: MLE + phasing error"
  ) +
  theme_bw()

cowplot::plot_grid(g1, g2)


## 4) False positivity

# optimizing
design_fp <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(p_error = c(0.1, 0.2, 0.3))
)

fp_p_error_null.o <- testFalsePositive(
  models = c("ase", "joint"), par = "p_error", phasing_error = "null",
  design_fp = design_fp, stan_models = stan_models, n_rep = 50
)
fp_p_error_k.o <- testFalsePositive(
  models = c("ase", "joint"), par = "p_error", phasing_error = "known",
  design_fp = design_fp, stan_models = stan_models, n_rep = 50
)
fp_p_error_u.o <- testFalsePositive(
  models = c("ase", "joint"), par = "p_error", phasing_error = "unknown",
  design_fp = design_fp, stan_models = stan_models, n_rep = 50
)

df_fp_p_error.o <- lapply(
  list(fp_p_error_null.o, fp_p_error_k.o, fp_p_error_u.o),
  function(fp) {
    fp %>%
      mutate(rep_idx = rep(rep(1:50, each = 200), 9)) %>%
      group_by(n_i, p_error, model, rep_idx) %>%
      summarise(fp = mean(p_val <= 0.05))
  }
) %>%
  do.call(what = "rbind") %>%
  ungroup() %>%
  mutate(mixture = rep(c("null", "known", "unknown"), each = 900)) %>%
  mutate(model = ifelse(mixture == "known", paste0(model, "_mix_k"), model)) %>%
  mutate(model = ifelse(mixture == "unknown", paste0(model, "_mix_u"), model))

g1 <- df_fp_p_error.o %>%
  mutate(model = factor(model, levels = c("ase", "joint", "ase_mix_k", "joint_mix_k", "ase_mix_u", "joint_mix_u"))) %>%
  ggplot(aes(x = model, y = fp, fill = model)) +
  geom_hline(yintercept = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  labs(title = "False Positivity: MLE + phasing error", y = "False Positive", x = "Model") +
  scale_fill_manual(values = pal_col) +
  theme_bw()


## 5) Power

# optimizing
design_power <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(p_error = c(0.1, 0.2, 0.3))
) %>% merge(y = data.frame(r = seq(0.01, 1.01, 0.05)))

power_p_error_null.o <- testPower(
  models = c("ase", "joint"), par = "p_error", phasing_error = "null",
  design_power = design_power, stan_models = stan_models
)
power_p_error_k.o <- testPower(
  models = c("ase", "joint"), par = "p_error", phasing_error = "known",
  design_power = design_power, stan_models = stan_models
)
power_p_error_u.o <- testPower(
  models = c("ase", "joint"), par = "p_error", phasing_error = "unknown",
  design_power = design_power, stan_models = stan_models
)

df_power_p_error.o <- lapply(
  list(power_p_error_null.o, power_p_error_k.o, power_p_error_u.o),
  function(power) {
    power %>%
      group_by(n_i, p_error, r_true, model) %>%
      summarise(power = mean(p_val <= 0.05))
  }
) %>%
  do.call(what = "rbind") %>%
  ungroup() %>%
  mutate(mixture = rep(c("null", "known", "unknown"), each = 378)) %>%
  mutate(model = ifelse(mixture == "known", paste0(model, "_mix_k"), model)) %>%
  mutate(model = ifelse(mixture == "unknown", paste0(model, "_mix_u"), model))

g2 <- df_power_p_error.o %>%
  mutate(model = factor(model, levels = c("ase", "joint", "ase_mix_k", "joint_mix_k", "ase_mix_u", "joint_mix_u"))) %>%
  ggplot(aes(x = r_true, y = power, color = model, group = model)) +
  geom_point() +
  geom_line() +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  labs(title = "Power: MLE + phasing error", y = "Power", x = "R (log FC)") +
  scale_color_manual(values = pal_col) +
  theme_bw()

cowplot::plot_grid(g1, g2)


## 6) ROC

# optimizing
roc <- rbind(
  power_p_error_null.o,
  power_p_error_k.o %>% mutate(model = paste0(model, "_mix_k")),
  power_p_error_u.o %>% mutate(model = paste0(model, "_mix_u")),
  fp_p_error_null.o %>% group_by(n_i, p_error, model) %>% sample_n(2100),
  fp_p_error_k.o %>% mutate(model = paste0(model, "_mix_k")) %>%
    group_by(n_i, p_error, model) %>% sample_n(2100),
  fp_p_error_u.o %>% mutate(model = paste0(model, "_mix_u")) %>%
    group_by(n_i, p_error, model) %>% sample_n(2100)
) %>%
  mutate(is_negative = r_true == 0)

df_roc <- lapply(seq(0, 1.01, 0.005), function(threshold) {
  roc %>%
    mutate(is_pred_positive = p_val < threshold) %>%
    group_by(n_i, p_error, model) %>%
    summarise(
      fp = mean(is_pred_positive[is_negative]),
      tp = mean(is_pred_positive[!is_negative])
    ) %>%
    mutate(p_thresh = threshold)
}) %>% do.call(what = "rbind")

df_roc %>%
  mutate(model = factor(model, levels = c("ase", "joint", "ase_mix_k", "joint_mix_k", "ase_mix_u", "joint_mix_u"))) %>%
  ggplot(aes(x = fp, y = tp, group = model, color = model)) +
  geom_line() +
  facet_grid(n_i ~ p_error, labeller = label_both) +
  scale_color_manual(values = pal_col) +
  labs(title = "ROC curve: MLE + phasing error", y = "True Positive Rate", x = "False Positive Rate") +
  theme_bw()
