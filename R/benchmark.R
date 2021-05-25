# helper function to test the sensitivity of the estimation of R
sensitivityR <- function(method = "sampling", return_posterior = "full",
                         significance_test = F, weighted_ase = FALSE, phasing_error = "null",
                         n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
                         phi = 3, theta = 30, baseline = 3, r = 0.5, ...) {
  n_rep <- ifelse(method == "sampling", 1, 200)

  out <- lapply(1:n_rep, function(i) {
    Y <- simulateCisEffect.s2s(
      n_i = n_i, maf = maf, prob_ref = prob_ref, prob_as = prob_as,
      phi = phi, theta = theta, baseline = baseline, r = r, ...
    )

    trc <- estimateCisRegEffects(
      data = Y$data, stan_models = stan_models,
      model = "trc", method = method, return_posterior = return_posterior,
      significance_test = significance_test,
      weighted_ase = weighted_ase, phasing_error = phasing_error
    )

    joint <- estimateCisRegEffects(
      data = Y$data, stan_models = stan_models,
      model = "joint", method = method, return_posterior = return_posterior,
      significance_test = significance_test,
      weighted_ase = weighted_ase, phasing_error = phasing_error
    )

    ase <- estimateCisRegEffects(
      data = Y$data, stan_models = stan_models,
      model = "ase", method = method, return_posterior = return_posterior,
      significance_test = significance_test,
      weighted_ase = weighted_ase, phasing_error = phasing_error
    )

    rbind(joint, trc, ase)
  }) %>% do.call(what = "rbind")

  return(out)
}

# helper function to test for the runtime
getRunTime <- function(method = "sampling", return_posterior = "full",
                       significance_test = F, weighted_ase = FALSE, phasing_error = "null",
                       n_i = 1000, maf = 0.1, prob_ref = 0.5, prob_as = 0.5,
                       phi = 3, theta = 30, baseline = 3, r = 0.5, ...) {
  out <- lapply(c("trc", "ase", "joint"), function(m) {
    Y <- simulateCisEffect.s2s(
      n_i = n_i, maf = maf, prob_ref = prob_ref, prob_as = prob_as,
      phi = phi, theta = theta, baseline = baseline, r = r, ...
    )
    system.time({
      estimateCisRegEffects(
        data = Y$data, stan_models = stan_models,
        model = m, method = method, return_posterior = return_posterior,
        significance_test = significance_test,
        weighted_ase = weighted_ase, phasing_error = phasing_error
      )
    })
  }) %>%
    do.call(what = "rbind") %>%
    as.data.frame() %>%
    mutate(model = c("trc", "ase", "joint"))

  return(out)
}

# helper function to test the fp of mle estimation
testFalsePositive <- function(design_fp, stan_models, n_rep = 50,
                              par = "baseline", weighted_ase = FALSE,
                              phasing_error = "null", confounders = "null") {
  out <- apply(design_fp[rep(1:nrow(design_fp), each = n_rep), ], 1, function(x) {
    sim_pars <- list(
      n_i = x["n_i"], r = 0, baseline = 3,
      maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
    )

    if (par == "baseline") {
      sim_pars$baseline <- x["b"]
    } else if (par == "p_error") {
      sim_pars$p_error <- x["p_error"]
    }

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
      mutate(n_i = x[1], !!par := x[2], r_true = 0)
  }) %>% do.call(what = "rbind")

  return(out)
}

# helper function to test the power of mle estimation
testPower <- function(design_power, stan_models, par = "baseline",
                      weighted_ase = FALSE, phasing_error = "null",
                      confounders = "null") {
  out <- apply(design_power, 1, function(x) {
    lapply(1:100, function(i) {
      sim_pars <- list(
        n_i = x["n_i"], r = x["r"], baseline = 3,
        maf = 0.1, prob_ref = 0.5, prob_as = 0.5, phi = 3, theta = 30
      )
      if (par == "baseline") {
        sim_pars$baseline <- x["b"]
      } else if (par == "p_error") {
        sim_pars$p_error <- x["p_error"]
      }

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
      mutate(n_i = x[1], !!par := x[2], r_true = x[3])
  }) %>% do.call(what = "rbind")

  return(out)
}
