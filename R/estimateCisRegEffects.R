#' Estimate cis-regulatory QTL effects
#'
#' @importFrom bridgesampling bridge_sampler bf
#' @importFrom rstan optimizing sampling extract
#' @export
estimateCisRegEffects <- function(data, stan_models,
                                  model = "joint", method = "optimizing", significance_test = TRUE,
                                  phasing_error = "null", weighted_ase = FALSE, confounders = "null",
                                  return_posterior = "full") {
  model <- match.arg(model, c("trc", "ase", "joint"))
  method <- match.arg(method, c("optimizing", "sampling"))
  phasing_error <- match.arg(phasing_error, c("null", "known", "unknown"))
  confounders <- match.arg(confounders, c("null", "fixed", "random"))
  return_posterior <- match.arg(return_posterior, c("full", "mean", "median"))

  run_models <- list()

  if (model == "joint") { ## joint model (default)
    if (confounders == "null") { # no confounders
      run_models <- ifelse(
        weighted_ase,
        c(stan_models$basic$joint_weighted, stan_models$null_models$joint_weighted_null) %>% list(),
        switch(phasing_error,
          "null" = c(stan_models$basic$joint, stan_models$null_models$joint_null) %>% list(),
          "known" = c(stan_models$phasing_error$joint_mixture_known, stan_models$null_models$joint_null) %>% list(),
          "unknown" = c(stan_models$phasing_error$joint_mixture_unknown, stan_models$null_models$joint_null) %>% list()
        )
      )
    } else if (confounders == "fixed") { # fixed-effect confounders
      run_models <- c(stan_models$confounders$joint_multi_fe, stan_models$null_models$joint_multi_fe_null) %>% list()
    } else { # random-effect confounders
      run_models <- c(stan_models$confounders$joint_indiv_re, stan_models$null_models$joint_indiv_re_null) %>% list()
    }
  } else if (model == "trc") { ## total read count component only
    run_models <- switch(confounders,
      "null" = c(stan_models$basic$trc, stan_models$null_models$trc_null) %>% list(),
      "fixed" = c(stan_models$confounders$trc_multi_fe, stan_models$null_models$trc_multi_fe_null) %>% list(),
      "random" = c(stan_models$confounders$trc_indiv_re, stan_models$null_models$trc_indiv_re_null) %>% list()
    )
  } else { ## allele-specific component only
    run_models <- ifelse(
      weighted_ase,
      c(stan_models$basic$ase_weighted, stan_models$null_models$ase_weighted_null) %>% list(),
      switch(phasing_error,
        "null" = c(stan_models$basic$ase, stan_models$null_models$ase_null) %>% list(),
        "known" = c(stan_models$phasing_error$ase_mixture_known, stan_models$null_models$ase_null) %>% list(),
        "unknown" = c(stan_models$phasing_error$ase_mixture_unknown, stan_models$null_models$ase_null) %>% list()
      )
    )
  }

  run_models <- run_models[[1]]

  if (method == "optimizing") { ## mle
    mle <- lapply(run_models, function(m) {
      rstan::optimizing(m, data = data, algorithm = "LBFGS", verbos = F, as_vector = F)
    })

    out <- data.frame(
      r_est = mle[[1]]$par$r,
      model = model
    )

    if (significance_test) { # likelihood ratio test
      out <- out %>% mutate(
        lr = -2 * (mle[[2]]$value - mle[[1]]$value),
        df = 1,
        p_val = 1 - pchisq(q = lr, df = df)
      )
    }
  } else { ## posterior
    if (!significance_test) {
      fit <- rstan::sampling(run_models[[1]], data = data, chain = 4, iter = 2000, refresh = 0)
      if (return_posterior == "full") {
        out <- data.frame(r_est = rstan::extract(fit, pars = "r")[[1]])
      } else {
        out <- data.frame(
          r_est = list(rstan::extract(fit, pars = "r")[[1]]) %>%
            do.call(what = return_posterior)
        )
      }
      out$model <- model
    } else { # bayes factor
      fit <- lapply(run_models, function(m) {
        rstan::sampling(m, data = data, chain = 4, iter = 2000, refresh = 0)
      })
      if (return_posterior == "full") {
        out <- data.frame(r_est = rstan::extract(fit[[1]], pars = "r")[[1]])
      } else {
        out <- data.frame(
          r_est = list(rstan::extract(fit[[1]], pars = "r")[[1]]) %>%
            do.call(what = return_posterior)
        )
      }

      out$model <- model
      out$bf <- bridgesampling::bf(
        x1 = bridgesampling::bridge_sampler(fit[[2]], silent = T),
        x2 = bridgesampling::bridge_sampler(fit[[1]], silent = T)
      )[["bf"]]
    }
  }

  return(out)
}
