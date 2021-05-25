library(tidyverse)
library(ggsci)
library(ggrepel)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load("./scripts/benchmark/basic.RData")

stan_models <- load_models()

pal_col <- c("ase" = pal_aaas()(10)[1],
             "joint" = pal_aaas()(10)[2],
             "trc" = pal_aaas()(10)[3])

## 1) Estimation accuracy

# sampling
r_est_ac.s <- data.frame()
for (r in c(0, 0.3, 0.6, 0.9)) {
  r_est <- sensitivityR(
    method = "sampling", return_posterior = "full",
    n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
    phi = 3, theta = 30, baseline = 3, r = r
  )

  r_est <- r_est %>% mutate(r_true = r)
  r_est_ac.s <- rbind(r_est_ac.s, r_est)
}

g1 <- r_est_ac.s %>%
  mutate(model = factor(model, c("ase", "trc", "joint"))) %>%
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
  r_est <- sensitivityR(
    method = "optimizing",
    n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
    phi = 3, theta = 30, baseline = 3, r = r
  )

  r_est <- r_est %>% mutate(r_true = r)
  r_est_ac.o <- rbind(r_est_ac.o, r_est)
}

g2 <- r_est_ac.o %>%
  mutate(model = factor(model, c("ase", "trc", "joint"))) %>%
  ggplot(aes(x = r_true, y = r_est, group = interaction(model, r_true), fill = model)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(y = r_true), color = "black", shape = 17, size = 3) +
  labs(x = "R_true", y = "R_est", title = "Accuracy: MLE") +
  ylim(c(-0.5, 1.5)) +
  scale_fill_manual(values = pal_col) +
  theme_bw()

cowplot::plot_grid(g1, g2)


## 2) Sensitivity (on sample size and baseline expression)

# sampling
r_est_sen.s <- data.frame()
for (b in c(1, 3, 5)) {
  for (n_i in c(100, 300, 500)){
    r_est = sensitivityR(
      method = "sampling", return_posterior = "full",
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = b, r = 0.5
    )

    r_est <- r_est %>% mutate(baseline = b, n_i = n_i)
    r_est_sen.s <- rbind(r_est_sen.s, r_est)
  }
}

g1 <- r_est_sen.s %>%
  mutate(model = factor(model, c("ase", "trc", "joint"))) %>%
  ggplot(aes(x = r_est, group = model, color = model)) +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  geom_vline(xintercept = 0.5, color = "black") +
  geom_density(size = 1) +
  labs(x = "R_est", y = "Posterior density", title = "Sensitivity: Posterior") +
  scale_color_manual(values = pal_col) +
  theme_bw()

# optimizing
r_est_sen.o <- data.frame()
for (b in c(1, 3, 5)) {
  for (n_i in c(100, 300, 500)){
    r_est = sensitivityR(
      method = "optimizing",
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = b, r = 0.5
    )

    r_est <- r_est %>% mutate(baseline = b, n_i = n_i)
    r_est_sen.o <- rbind(r_est_sen.o, r_est)
  }
}

g2 <- r_est_sen.o %>%
  mutate(model = factor(model, c("ase", "trc", "joint"))) %>%
  ggplot(aes(x = model, y = r_est, group = model, fill = model)) +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0.5) +
  labs(x = "R_true", y = "R_est", title = "Sensitivity: MLE") +
  ylim(c(-0.5, 1)) +
  scale_fill_manual(values = pal_col) +
  theme_bw()

cowplot::plot_grid(g1, g2)


## 3) Run time

# sampling
runtime.s <- data.frame()
for (n_i in c(100, 400, 800, 1600, 3200, 6400)){
  runtime <- lapply(1:3, function(i){
    runtime_notest <- getRunTime(
      method = "sampling", significance_test = F,
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5
    )
    runtime_test <- getRunTime(
      method = "sampling", significance_test = T,
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5
    )
    rbind(runtime_notest, runtime_test) %>%
      mutate(sig_test = rep(c(0, 1), each = 3), n_i = n_i)
  }) %>% do.call(what = "rbind")

  runtime.s <- rbind(runtime.s, runtime)
}

g1 <- runtime.s %>%
  ggplot(aes(x = n_i, y = elapsed, group = interaction(model, sig_test), color = model, shape = factor(sig_test))) +
  geom_point() +
  stat_smooth(aes(linetype = factor(sig_test)), method = "lm", formula = y ~ x, se = F) +
  scale_color_manual(values = pal_col) +
  ylim(c(0, 50)) +
  labs(x = "Number of samples", y = "Total Run Time (seconds)",
       color = "Model", shape = "+ LRT", linetype = "+ LRT",
       title = "Runtime: Posterior") +
  theme_bw()


# optimizing
runtime.o <- data.frame()
for (n_i in c(100, 400, 800, 1600, 3200, 6400)){
  runtime <- lapply(1:100, function(i){
    runtime_notest <- getRunTime(
      method = "optimizing", significance_test = F,
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5
    )
    runtime_test <- getRunTime(
      method = "optimizing", significance_test = T,
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3, r = 0.5
    )
    rbind(runtime_notest, runtime_test) %>%
      mutate(sig_test = rep(c(0, 1), each = 3), n_i = n_i)
  }) %>% do.call(what = "rbind")

  runtime.o <- rbind(runtime.o, runtime)
}

g2 <- runtime.o %>%
  ggplot(aes(x = n_i, y = elapsed, group = interaction(model, sig_test), color = model, shape = factor(sig_test))) +
  geom_point() +
  stat_smooth(aes(linetype = factor(sig_test)), method = "lm", formula = y ~ x, se = F) +
  scale_color_manual(values = pal_col) +
  ylim(c(0, 0.3)) +
  labs(x = "Number of samples", y = "Total Run Time (seconds)",
       color = "Model", shape = "+ LRT", linetype = "+ LRT",
       title = "Runtime: MLE") +
  theme_bw()

cowplot::plot_grid(g1, g2)

## 4) False positivity

# optimizing
design_fp <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(b = c(1, 3, 5))
)

fp.o <- testFalsePositive(design_fp = design_fp, stan_models = stan_models, n_rep = 50)
df_fp.o <- fp.o %>%
  mutate(rep_idx = rep(rep(1:50, each = 300), 9), model = factor(model, c("ase", "trc", "joint"))) %>%
  group_by(n_i, baseline, model, rep_idx) %>%
  summarise(fp = mean(p_val <= 0.05))
g1 <- df_fp.o %>% ggplot(aes(x = model, y = fp, fill = model)) +
  geom_hline(yintercept = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  labs(title = "False Positivity: MLE", y = "False Positive", x = "Model") +
  scale_fill_manual(values = pal_col) +
  theme_bw()

## 5) Power

# optimizing
design_power <- merge(
  data.frame(n_i = c(100, 300, 500)),
  data.frame(b = c(1, 3, 5))
) %>% merge(y = data.frame(r = seq(0.01, 1.01, 0.05)))

power.o <- testPower(design_power = design_power, stan_models = stan_models)
df_power.o <- power.o %>%
  group_by(n_i, baseline, r_true, model) %>%
  summarise(power = mean(p_val <= 0.05)) %>%
  mutate(model = factor(model, c("ase", "trc", "joint")))
g2 <- df_power.o %>% ggplot(aes(x = r_true, y = power, color = model, group = model)) +
  geom_point() +
  geom_line() +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  labs(title = "Power: MLE", y = "Power", x = "R (log FC)") +
  scale_color_manual(values = pal_col) +
  theme_bw()

cowplot::plot_grid(g1, g2)

## 6) ROC

# optimizing
roc <- rbind(power.o, fp.o %>% group_by(n_i, baseline, model) %>% sample_n(2100)) %>%
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
  mutate(model = factor(model, c("ase", "trc", "joint"))) %>%
  ggplot(aes(x = fp, y = tp, group = model, color = model)) +
  geom_line() +
  facet_grid(n_i ~ baseline, labeller = label_both) +
  scale_color_manual(values = pal_col) +
  labs(title = "ROC curve: MLE", y = "True Positive Rate", x = "False Positive Rate") +
  theme_bw()
