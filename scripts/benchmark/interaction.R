library(tidyverse)
library(ggsci)
library(ggrepel)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load("./scripts/benchmark/basic.RData")

stan_models <- load_models()

pal_col <- c(
  "ase" = pal_aaas()(10)[1],
  "joint" = pal_aaas()(10)[2],
  "trc" = pal_aaas()(10)[3],
  "ase_int" = pal_npg()(10)[2],
  "joint_int" = pal_npg()(10)[5],
  "trc_int" = pal_npg()(10)[9]
)


## 1) Estimation accuracy

# simulate test data
n_i <- 1000
n_cond <- 10
cond <- sample(1:n_cond, n_i, replace = T)
r.cond.effect <- 1:n_cond * 0.1
Y <- simulateCisEffect.s2s(
  n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
  phi = 3, theta = 30, baseline = 3,
  cond = cond, r.cond.effect = r.cond.effect
)

# sampling
fit_trc_int <- sampling(
  stan_models$interaction$trc_int,
  data = Y$data,
  chain = 4, iter = 2000, refresh = 0,
  include = FALSE, pars = c("mu_t", "mu_a", "log_lik", "sum_log_lik")
)
fit_ase_int <- sampling(
  stan_models$interaction$ase_int,
  data = Y$data,
  chain = 4, iter = 2000, refresh = 0,
  include = FALSE, pars = c("mu_t", "mu_a", "log_lik", "sum_log_lik")
)
fit_joint_int <- sampling(
  stan_models$interaction$joint_int,
  data = Y$data,
  chain = 4, iter = 2000, refresh = 0,
  include = FALSE, pars = c("mu_t", "mu_a", "log_lik", "sum_log_lik")
)

r_est_ac.s <- data.frame(
  cond = rep(1:n_cond, each = 4000),
  r_true = rep(r.cond.effect, each = 4000),
  r_est = sapply(list(fit_trc_int, fit_ase_int, fit_joint_int), function(x) extract(x, pars = "R")[[1]] %>% as.vector()) %>% as.vector(),
  model = rep(c("trc_int", "ase_int", "joint_int"), each = 4000 * n_cond)
)

g1 <- r_est_ac.s %>%
  mutate(cond = factor(cond), model = factor(model, c("ase_int", "trc_int", "joint_int"))) %>%
  ggplot(aes(x = cond, y = r_est, group = interaction(model, cond), fill = model)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(y = r_true), color = "black", shape = 17, size = 3) +
  scale_fill_manual(values = pal_col) +
  labs(x = "Condition", y = "Regulatory effect", title = "Accuracy: Posterior") +
  ylim(c(-1, 2)) +
  theme_bw()

# optimizing
r_est_ac.o <- data.frame()
for (i in 1:200) {
  Y <- simulateCisEffect.s2s(
    n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
    phi = 3, theta = 30, baseline = 3,
    cond = cond, r.cond.effect = r.cond.effect
  )

  mle_trc_int <- optimizing(stan_models$interaction$trc_int, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
  mle_ase_int <- optimizing(stan_models$interaction$ase_int, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
  mle_joint_int <- optimizing(stan_models$interaction$joint_int, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)

  r_est_i <- data.frame(
    cond = 1:n_cond,
    r_true = r.cond.effect,
    r_est = sapply(list(mle_trc_int, mle_ase_int, mle_joint_int), function(x) x$par$R) %>% as.vector(),
    model = rep(c("trc_int", "ase_int", "joint_int"), each = n_cond)
  )

  r_est_ac.o <- rbind(r_est_ac.o, r_est_i)
}

g2 <- r_est_ac.o %>%
  mutate(cond = factor(cond), model = factor(model, c("ase_int", "trc_int", "joint_int"))) %>%
  ggplot(aes(x = cond, y = r_est, group = interaction(model, cond), fill = model)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(y = r_true), color = "black", shape = 17, size = 3) +
  scale_fill_manual(values = pal_col) +
  labs(x = "Condition", y = "Regulatory effect", title = "Accuracy: MLE") +
  ylim(c(-1, 2)) +
  theme_bw()

cowplot::plot_grid(g1, g2)


## 2) Sensitivity (on sample size)
n_cond <- 5
r.cond.effect <- 1:n_cond * 0.2

# sampling
r_est_sen.s <- lapply(c(500, 1000, 5000), function(n_i) {
  cond <- sample(1:n_cond, n_i, replace = T)
  Y <- simulateCisEffect.s2s(
    n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
    phi = 3, theta = 30, baseline = 3,
    cond = cond, r.cond.effect = r.cond.effect
  )

  fit_trc_int <- sampling(stan_models$interaction$trc_int,
    data = Y$data, chain = 4, iter = 2000, refresh = 0,
    include = FALSE, pars = c("mu_t", "mu_a", "log_lik", "sum_log_lik")
  )
  fit_ase_int <- sampling(stan_models$interaction$ase_int,
    data = Y$data, chain = 4, iter = 2000, refresh = 0,
    include = FALSE, pars = c("mu_t", "mu_a", "log_lik", "sum_log_lik")
  )
  fit_joint_int <- sampling(stan_models$interaction$joint_int,
    data = Y$data, chain = 4, iter = 2000, refresh = 0,
    include = FALSE, pars = c("mu_t", "mu_a", "log_lik", "sum_log_lik")
  )

  data.frame(
    cond = rep(1:n_cond, each = 4000),
    r_true = rep(r.cond.effect, each = 4000),
    n_i = n_i,
    r_est = sapply(list(fit_trc_int, fit_ase_int, fit_joint_int), function(x) extract(x, pars = "R")[[1]] %>% as.vector()) %>% as.vector(),
    model = rep(c("trc_int", "ase_int", "joint_int"), each = 4000 * n_cond)
  )
}) %>% do.call(what = "rbind")

g1 <- r_est_sen.s %>%
  mutate(cond = factor(cond), model = factor(model, c("ase_int", "trc_int", "joint_int"))) %>%
  ggplot(aes(x = cond, y = r_est, group = interaction(model, cond), fill = model)) +
  facet_wrap(~n_i, nrow = 1, labeller = label_both) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(y = r_true), color = "black", shape = 17, size = 3) +
  scale_fill_manual(values = pal_col) +
  labs(x = "Condition", y = "Regulatory effect", title = "Sensitivity: Posterior") +
  ylim(c(-1, 2)) +
  theme_bw()

# optimizing
r_est_sen.o <- lapply(c(500, 1000, 5000), function(n_i) {
  cond <- sample(1:n_cond, n_i, replace = T)

  r_est <- lapply(1:200, function(i) {
    Y <- simulateCisEffect.s2s(
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3,
      cond = cond, r.cond.effect = r.cond.effect
    )

    mle_trc_int <- optimizing(stan_models$interaction$trc_int, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
    mle_ase_int <- optimizing(stan_models$interaction$ase_int, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
    mle_joint_int <- optimizing(stan_models$interaction$joint_int, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)

    data.frame(
      cond = 1:n_cond,
      r_true = r.cond.effect,
      r_est = sapply(list(mle_trc_int, mle_ase_int, mle_joint_int), function(x) x$par$R) %>% as.vector(),
      model = rep(c("trc_int", "ase_int", "joint_int"), each = n_cond)
    )
  }) %>%
    do.call(what = "rbind") %>%
    mutate(n_i = n_i)
}) %>% do.call(what = "rbind")

g2 <- r_est_sen.o %>%
  mutate(cond = factor(cond), model = factor(model, c("ase_int", "trc_int", "joint_int"))) %>%
  ggplot(aes(x = cond, y = r_est, group = interaction(model, cond), fill = model)) +
  facet_wrap(~n_i, nrow = 1, labeller = label_both) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(y = r_true), color = "black", shape = 17, size = 3) +
  scale_fill_manual(values = pal_col) +
  labs(x = "Condition", y = "Regulatory effect", title = "Sensitivity: MLE") +
  ylim(c(-1, 2)) +
  theme_bw()

cowplot::plot_grid(g1, g2)


## 4) False positivity

# optimizing
design_fp <- merge(
  data.frame(n_i = c(500, 1000, 5000)),
  data.frame(n_cond = c(1, 5, 10))
)

fp_int.o <- apply(design_fp[rep(1:nrow(design_fp), each = 10), ], 1, function(x) {
  n_i <- x["n_i"]
  n_cond <- x["n_cond"]
  cond <- sample(1:n_cond, n_i, replace = T)
  r.cond.effect <- rep(0, n_cond)

  lapply(1:50, function(i){
    Y <- simulateCisEffect.s2s(
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3,
      cond = cond, r.cond.effect = r.cond.effect
    )

    mle_trc <- optimizing(stan_models$basic$trc, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
    mle_ase <- optimizing(stan_models$basic$ase, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
    mle_joint <- optimizing(stan_models$basic$joint, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)

    mle_trc_int <- optimizing(stan_models$interaction$trc_int, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
    mle_ase_int <- optimizing(stan_models$interaction$ase_int, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
    mle_joint_int <- optimizing(stan_models$interaction$joint_int, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)

    mle_trc_null <- optimizing(stan_models$null_models$trc_null, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
    mle_ase_null <- optimizing(stan_models$null_models$ase_null, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
    mle_joint_null <- optimizing(stan_models$null_models$joint_null, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)

    data.frame(
      r_est = c(
        mle_trc$par$r, mle_ase$par$r, mle_joint$par$r,
        mle_trc_int$par$R[1], mle_ase_int$par$R[1], mle_joint_int$par$R[1]
      ),
      model = c("trc", "ase", "joint", "trc_int", "ase_int", "joint_int")
    ) %>% mutate(
      lr = 2 * c(
        mle_trc$par$sum_log_lik - mle_trc_null$par$sum_log_lik,
        mle_ase$par$sum_log_lik - mle_ase_null$par$sum_log_lik,
        mle_joint$par$sum_log_lik - mle_joint_null$par$sum_log_lik,
        mle_trc_int$par$sum_log_lik - mle_trc_null$par$sum_log_lik,
        mle_ase_int$par$sum_log_lik - mle_ase_null$par$sum_log_lik,
        mle_joint_int$par$sum_log_lik - mle_joint_null$par$sum_log_lik
      ),
      df = c(1, 1, 1, n_cond, n_cond, n_cond),
      p_val = 1 - pchisq(q = lr, df = df)
    )
  }) %>% do.call(what = "rbind") %>% mutate(n_i = n_i, n_cond = n_cond)
}) %>% do.call(what = "rbind")

df_fp_int.o <- fp_int.o %>%
  mutate(rep_idx = rep(rep(1:10, each = 300), 9)) %>%
  group_by(n_i, n_cond, model, rep_idx) %>%
  summarise(fp = mean(p_val <= 0.05))

g1 <- df_fp_int.o %>%
  mutate(model = factor(model, levels = c("ase", "trc", "joint", "ase_int", "trc_int", "joint_int"))) %>%
  ggplot(aes(x = model, y = fp, fill = model)) +
  geom_hline(yintercept = 0.05) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(n_i ~ n_cond, labeller = label_both) +
  labs(title = "False Positivity: MLE", y = "False Positive", x = "Model") +
  scale_fill_manual(values = pal_col) +
  theme_bw()
cowplot::plot_grid(g1, g2)


## 5) Power

# optimizing
design_power <- merge(
  data.frame(n_i = c(500, 1000, 5000)),
  data.frame(n_cond = c(1, 5, 10))
) %>% merge(y = data.frame(r = seq(0.01, 1.01, 0.05)))

power_int.o <- apply(design_power, 1, function(x) {
  n_i <- x["n_i"]
  n_cond <- x["n_cond"]
  cond <- sample(1:n_cond, n_i, replace = T)
  r.cond.effect <- c(x["r"], rep(0, n_cond - 1))

  lapply(1:20, function(i){
    Y <- simulateCisEffect.s2s(
      n_i = n_i, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
      phi = 3, theta = 30, baseline = 3,
      cond = cond, r.cond.effect = r.cond.effect
    )

    mle_trc <- optimizing(stan_models$basic$trc, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
    mle_ase <- optimizing(stan_models$basic$ase, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
    mle_joint <- optimizing(stan_models$basic$joint, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)

    mle_trc_int <- optimizing(stan_models$interaction$trc_int, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
    mle_ase_int <- optimizing(stan_models$interaction$ase_int, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
    mle_joint_int <- optimizing(stan_models$interaction$joint_int, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)

    mle_trc_null <- optimizing(stan_models$null_models$trc_null, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
    mle_ase_null <- optimizing(stan_models$null_models$ase_null, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)
    mle_joint_null <- optimizing(stan_models$null_models$joint_null, data = Y$data, algorithm = "LBFGS", verbos = F, as_vector = F)

    data.frame(
      r_est = c(
        mle_trc$par$r, mle_ase$par$r, mle_joint$par$r,
        mle_trc_int$par$R[1], mle_ase_int$par$R[1], mle_joint_int$par$R[1]
      ),
      model = c("trc", "ase", "joint", "trc_int", "ase_int", "joint_int")
    ) %>% mutate(
      lr = 2 * c(
        mle_trc$par$sum_log_lik - mle_trc_null$par$sum_log_lik,
        mle_ase$par$sum_log_lik - mle_ase_null$par$sum_log_lik,
        mle_joint$par$sum_log_lik - mle_joint_null$par$sum_log_lik,
        mle_trc_int$par$sum_log_lik - mle_trc_null$par$sum_log_lik,
        mle_ase_int$par$sum_log_lik - mle_ase_null$par$sum_log_lik,
        mle_joint_int$par$sum_log_lik - mle_joint_null$par$sum_log_lik
      ),
      df = c(1, 1, 1, n_cond, n_cond, n_cond),
      p_val = 1 - pchisq(q = lr, df = df)
    )
  }) %>% do.call(what = "rbind") %>% mutate(n_i = n_i, n_cond = n_cond, r_true = x["r"])
}) %>% do.call(what = "rbind")

df_power_int.o <- power_int.o %>%
  group_by(n_i, n_cond, r_true, model) %>%
  summarise(power = mean(p_val <= 0.05))
g2 <- df_power_int.o %>%
  mutate(model = factor(model, levels = c("ase", "trc", "joint", "ase_int", "trc_int", "joint_int"))) %>%
  ggplot(aes(x = r_true, y = power, color = model, group = model)) +
  geom_point() +
  geom_line() +
  facet_grid(n_i ~ n_cond, labeller = label_both) +
  labs(title = "Power: MLE", y = "Power", x = "R (log FC)") +
  scale_color_manual(values = pal_col) +
  theme_bw()

cowplot::plot_grid(g1, g2)


## 6) ROC

# optimizing
roc <- rbind(power_int.o, fp_int.o %>% mutate(r_true = 0)) %>%
  mutate(is_negative = r_true == 0)

df_roc <- lapply(seq(0, 1.01, 0.01), function(threshold) {
  roc %>%
    mutate(is_pred_positive = p_val < threshold) %>%
    group_by(n_i, n_cond, model) %>%
    summarise(
      fp = mean(is_pred_positive[is_negative]),
      tp = mean(is_pred_positive[!is_negative])
    ) %>%
    mutate(p_thresh = threshold)
}) %>% do.call(what = "rbind")

df_roc %>%
  mutate(model = factor(model, levels = c("ase", "trc", "joint", "ase_int", "trc_int", "joint_int"))) %>%
  ggplot(aes(x = fp, y = tp, group = model, color = model)) +
  geom_line() +
  facet_grid(n_i ~ n_cond, labeller = label_both) +
  scale_color_manual(values = pal_col) +
  labs(title = "ROC curve: MLE + interaction", y = "True Positive Rate", x = "False Positive Rate") +
  theme_bw()

