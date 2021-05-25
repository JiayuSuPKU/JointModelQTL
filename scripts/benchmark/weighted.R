library(tidyverse)
library(ggsci)
library(ggrepel)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load("./scripts/benchmark/basic.RData")
load("./scripts/benchmark/weighted.RData")

stan_models <- load_models()

pal_col <- c(
  "ase" = pal_aaas()(10)[1],
  "joint" = pal_aaas()(10)[2],
  "trc" = pal_aaas()(10)[3],
  "joint_weighted" = pal_npg()(10)[2],
  "ase_weighted" = pal_npg()(10)[5]
)


## 2) Sensitivity (on sample size and baseline expression)

# sampling
r_est_sen.s <- data.frame()
for (b in c(1, 3, 5)) {
  for (n_i in c(100, 300, 500)) {
    r_est <- sensitivityR(
      models = c("ase", "joint"),
      method = "sampling", return_posterior = "full",
      weighted_ase = F,
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = b, r = 0.5
    ) %>% mutate(baseline = b, n_i = n_i)

    r_est_w <- sensitivityR(
      models = c("ase", "joint"),
      method = "sampling", return_posterior = "full",
      weighted_ase = T,
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = b, r = 0.5
    ) %>% mutate(model = paste0(model, "_weighted"), baseline = b, n_i = n_i)

    r_est_sen.s <- rbind(r_est_sen.s, r_est, r_est_w)
  }
}

g1 <- r_est_sen.s %>%
  mutate(model = factor(model, c("ase", "joint", "ase_weighted", "joint_weighted"))) %>%
  ggplot(aes(x = r_est, group = model, color = model)) +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  geom_vline(xintercept = 0.5, color = "black") +
  geom_density(size = 1) +
  labs(x = "R_est", y = "Posterior density", title = "Sensitivity: Posterior + weighted") +
  scale_color_manual(values = pal_col) +
  theme_bw()

# optimizing
r_est_sen.o <- data.frame()
for (b in c(1, 3, 5)) {
  for (n_i in c(100, 300, 500)) {
    r_est <- sensitivityR(
      models = c("ase", "joint"),
      method = "optimizing", weighted_ase = F,
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = b, r = 0.5
    ) %>% mutate(baseline = b, n_i = n_i)

    r_est_w <- sensitivityR(
      models = c("ase", "joint"),
      method = "optimizing", weighted_ase = T,
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = b, r = 0.5
    ) %>% mutate(model = paste0(model, "_weighted"), baseline = b, n_i = n_i)

    r_est_sen.o <- rbind(r_est_sen.o, r_est, r_est_w)
  }
}

g2 <- r_est_sen.o %>%
  mutate(model = factor(model, c("ase", "joint", "ase_weighted", "joint_weighted"))) %>%
  ggplot(aes(x = model, y = r_est, group = model, fill = model)) +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0.5) +
  labs(x = "R_true", y = "R_est", title = "Sensitivity: MLE + weighted") +
  ylim(c(-0.5, 1)) +
  scale_fill_manual(values = pal_col) +
  theme_bw()

cowplot::plot_grid(g1, g2)


## 4) False positivity

# optimizing
design_fp <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(b = c(1, 3, 5))
)

fp_weighted.o <- testFalsePositive(
  models = c("ase", "joint"), weighted_ase = T,
  design_fp = design_fp, stan_models = stan_models, n_rep = 50
)
df_fp_weighted.o <- fp_weighted.o %>%
  mutate(rep_idx = rep(rep(1:50, each = 200), 9)) %>%
  mutate(model = paste0(model, "_weighted")) %>%
  group_by(n_i, baseline, model, rep_idx) %>%
  summarise(fp = mean(p_val <= 0.05))

g1 <- rbind(df_fp.o, df_fp_weighted.o) %>%
  mutate(model = factor(model, c("ase", "trc", "joint", "ase_weighted", "joint_weighted"))) %>%
  ggplot(aes(x = model, y = fp, fill = model)) +
  geom_hline(yintercept = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  labs(title = "False Positivity: MLE + weighted", y = "False Positive", x = "Model") +
  scale_fill_manual(values = pal_col) +
  theme_bw()


## 5) Power

# optimizing
design_power <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(b = c(1, 3, 5))
) %>% merge(y = data.frame(r = seq(0.01, 1.01, 0.05)))

power_weighted.o <- testPower(
  models = c("ase", "joint"), weighted_ase = T,
  design_power = design_power, stan_models = stan_models
)
df_power_weighted.o <- power_weighted.o %>%
  mutate(model = paste0(model, "_weighted")) %>%
  group_by(n_i, baseline, r_true, model) %>%
  summarise(power = mean(p_val <= 0.05))

g2 <- rbind(df_power.o, df_power_weighted.o) %>%
  mutate(model = factor(model, c("ase", "trc", "joint", "ase_weighted", "joint_weighted"))) %>%
  ggplot(aes(x = r_true, y = power, color = model, group = model)) +
  geom_point() +
  geom_line() +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  labs(title = "Power: MLE + weighted", y = "Power", x = "R (log FC)") +
  scale_color_manual(values = pal_col) +
  theme_bw()

cowplot::plot_grid(g1, g2)


## 6) ROC

# optimizing
roc <- rbind(
  power.o,
  power_weighted.o %>% mutate(model = paste0(model, "_weighted")),
  fp.o %>% group_by(n_i, baseline, model) %>% sample_n(2100),
  fp_weighted.o %>% mutate(model = paste0(model, "_weighted")) %>% group_by(n_i, baseline, model) %>% sample_n(2100)
) %>%
  mutate(is_negative = r_true == 0)

df_roc <- lapply(seq(0, 1.01, 0.01), function(threshold) {
  roc %>%
    mutate(is_pred_positive = p_val < threshold) %>%
    group_by(n_i, baseline, model) %>%
    summarise(
      fp = mean(is_pred_positive[is_negative]),
      tp = mean(is_pred_positive[!is_negative])
    ) %>%
    mutate(p_thresh = threshold)
}) %>% do.call(what = "rbind")

df_roc %>%
  mutate(model = factor(model, c("ase", "trc", "joint", "ase_weighted", "joint_weighted"))) %>%
  ggplot(aes(x = fp, y = tp, group = model, color = model)) +
  geom_line() +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  scale_color_manual(values = pal_col) +
  labs(title = "ROC curve: MLE + weighted", y = "True Positive Rate", x = "False Positive Rate") +
  theme_bw()
