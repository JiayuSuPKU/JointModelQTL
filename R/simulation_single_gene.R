#' Simulate read for one gene
#'
#' This function simulates read counts data from negative binomial (total read counts)
#' and beta-binomial (allele-specific read counts) distributions for one
#'
#' J=1, K=1, L=1
#'
#' \code{T ~ NB(mu, phi)}
#' \code{A_alt ~ BB(T * prob_as := N_as, prob_alt, theta)}
#'
#' @param n_i Number of samples
#' @param mu Expected number of total reads
#' @param phi Overdispersion term in negative-binomial
#' @param prob_as Proportion of allele-specific reads
#' @param prob_alt Proportion of alt-allele reads in allele-specific reads
#' @param theta Overdispersion term in beta-binomial
#'
#' @importFrom extraDistr rbbinom
#' @importFrom stats rbinom rnbinom runif
#'
#' @export
simulateReadCounts.s <- function(n_i, mu, phi, prob_as, prob_alt, theta) {

  # total read counts
  trc <- rnbinom(n = n_i, mu = mu, size = phi)

  # allele-specific read counts
  n_as <- round(trc * prob_as)

  # alternative allele read counts
  a_alt <- rbbinom(
    n = n_i, size = n_as,
    alpha = prob_alt * (theta - 1),
    beta = (1 - prob_alt) * (theta - 1)
  )

  return(list(TotalReadCounts = trc, RefCounts = n_as - a_alt, AltCounts = a_alt))
}


#' Simulate genotype for one gene-snp pair
#'
#' This function simulates genotype G and phasing P for one gene-snp pair
#'
#' J=1, K=1, L=1
#'
#' @param n_i Number of samples
#' @param maf Minor allele frequency of the test snp
#' @param prob_ref Probability of the test snp on the ref allele of the target
#' gene if the test snp is heterogeneous
#' @param origin A vector specifying sample origin. If \code{NULL} then no two
#' samples come from a same individual.
#'
#' @export
simulateGenotype.s2s <- function(n_i, maf, prob_ref, origin = NULL) {
  maf <- ifelse(maf <= 0.5, maf, 1 - maf)

  if (is.null(origin)) {
    # all samples come from different individuals
    origin <- 1:n_i
  }

  # simulate two loci for each individual
  n_indiv <- length(unique(origin))
  stopifnot(n_indiv > 1)

  origin <- as.numeric(factor(origin))

  h1 <- rbinom(n = n_indiv, size = 1, prob = maf)
  h2 <- rbinom(n = n_indiv, size = 1, prob = maf)
  # G: n_indiv by 1, genotype matrix
  G <- h1 + h2

  # P: n_indiv by 1, phasing matrix
  P <- ifelse(G == 1, 1, 0)
  is_cis <- (runif(n_indiv) <= prob_ref) * 2 - 1
  P <- P * is_cis

  # calculate G and P for each sample
  G_tilta <- G[origin]
  P_tilta <- P[origin]

  return(list(
    Genotype = G_tilta, Phasing = P_tilta,
    Origin = origin, N_indiv = n_indiv
  ))
}


softmax <- function(x) {
  return(exp(x) / (1 + exp(x)))
}


#' Simulate cis-regulatory effect for one gene-snp pair
#'
#' This function simulates genotype G, phasing P, and expression
#' (total and allele-specific) profiles for one gene-snp pair
#'
#' J=1, K=1, L=1
#'
#' @inheritParams simulateReadCounts.s
#' @inheritParams simulateGenotype.s2s
#' @param baseline Baseline expression (under ref/ref scenario)
#' @param r Cis-regulatory effect defined as the log read ratio under alt/alt
#' and ref/ref scenario
#' @param p_error A vector of phasing error rate of the gene pair in each
#' sample. If \code{length(p_error) == 1} then all samples share the same
#' phasing error rate
#' @param origin.effect A vector specifying the effect of sample origin with
#' length \code{n_indiv} (the number of samples). If \code{NULL} then no effect
#' will be added
#' @param lib.size A vector of the relative library size of each sample. If
#' \code{NULL} then no adjustment will be made
#' @param confounder A n_i by N matrix where N is the number of confounding factors
#' @param confounder.effect A vector specifying the effect of each confounding factor
#' If \code{NULL} then no effect will be added
#'
#' @export
simulateCisEffect.s2s <- function(n_i, maf, prob_ref, phi, prob_as,
                                  theta, baseline, r, p_error = NULL,
                                  origin = NULL, origin.effect = NULL,
                                  lib.size = NULL,
                                  confounder = NULL, confounder.effect = NULL) {

  # simulate genotype
  meta <- simulateGenotype.s2s(
    n_i = n_i, maf = maf, prob_ref = prob_ref, origin = origin
  )

  # simulate genetic effects on total read counts
  log_mu <- baseline + ifelse(meta$Genotype == 1, log(1 + exp(r)) - log(2), 0) +
    ifelse(meta$Genotype == 2, r, 0)

  # simulate individual effect on total read counts
  if (!is.null(origin.effect)) {
    log_mu <- log_mu + origin.effect[meta$Origin]
  }

  # simulate library size effect on total read counts
  if (!is.null(lib.size)) {
    log_mu <- log_mu + lib.size
  }

  # simulate confounding effect on total read counts
  if (!is.null(confounder.effect)) {
    # check confounder dimensions
    stopifnot(length(confounder.effect) == ncol(confounder))

    log_mu <- log_mu + confounder %*% confounder.effect
  }

  # simulate genetic effects on allelic imbalance
  logit_prob_alt <- meta$Phasing * r

  # simulate expression profile
  Y <- simulateReadCounts.s(
    n_i = n_i, mu = exp(log_mu), phi = phi, prob_as = prob_as,
    prob_alt = softmax(logit_prob_alt), theta = theta
  )

  # simulate phasing error for each sample
  if (is.null(p_error)) {
    p_error <- rep(0, n_i)
  } else {
    if (length(p_error) == 1) {
      p_error <- rep(p_error, n_i)
    }
    # reverse P if error occurs
    is_p_error <- runif(n_i) <= p_error
    meta$Phasing[is_p_error] <- (-1) * meta$Phasing[is_p_error]
  }

  sim <- list(
    genotype_pars = list(
      n_i = n_i,
      maf = maf,
      prob_ref = prob_ref
    ),
    gene_pars = list(
      prob_as = prob_as,
      phi = phi,
      theta = theta,
      baseline = baseline,
      r = r
    ),
    conf_pars = list(
      origin.effect = origin.effect,
      lib.size = lib.size,
      confounder.effect = confounder.effect
    ),
    data = list(
      I = n_i,
      N_indiv = meta$N_indiv,
      N = ncol(confounder),
      P = meta$Phasing,
      P_error = p_error,
      G = meta$Genotype,
      T = Y$TotalReadCounts,
      log1p_T = log1p(Y$TotalReadCounts),
      A_ref = Y$RefCounts,
      A_alt = Y$AltCounts,
      Is_ase_het = Y$RefCounts * Y$AltCounts != 0,
      Ori = meta$Origin,
      Cov_ori = outer(meta$Origin, meta$Origin, FUN = "==") * 1,
      Lib_size = lib.size,
      X = confounder
    )
  )

  sim$data$logit_pi_alt <- ifelse(
    !sim$data$Is_ase_het, 0, log(sim$data$A_alt / sim$data$A_ref)
  )

  return(sim)
}

#' Simulate genotype for one-to-multi gene-snp pairs
#'
#' This function simulates genotype G and phasing P for one-to-multi gene-snp pairs
#'
#' J=1, K>=1, L=1
#'
#' @param n_i Number of samples
#' @param n_k Number of regulatory snps
#' @param snp_pars A list of length \code{n_k} or \code{1} specifying
#' generative parameters for each regulatory snp. Each element should be a
#' list containing parameters used in \code{\link{simulateGenotype.s2s}}.
#' If \code{length(snp_pars) == 1} then all \code{n_k} snps share the same set
#' of parameters.
#' @param origin A vector specifying sample origin. If \code{NULL} then no two
#' samples come from a same individual.
#'
#' @export
simulateGenotype.s2m <- function(n_i, n_k, snp_pars, origin = NULL) {
  if (length(snp_pars) == 1) {
    # all snps share same pars
    snp_pars <- rep(snp_pars, n_k)
  }

  if (is.null(origin)) {
    # all samples come from different individuals
    origin <- 1:n_i
  }

  # simulate genotype for each individual
  n_indiv <- length(unique(origin))
  stopifnot(n_indiv > 1)
  origin <- as.numeric(factor(origin))

  # G: n_indiv by n_k, genotype matrix
  G <- lapply(snp_pars, function(x) {
    h1 <- rbinom(n = n_indiv, size = 1, prob = x$maf)
    h2 <- rbinom(n = n_indiv, size = 1, prob = x$maf)
    return(h1 + h2)
  }) %>% do.call(what = "cbind")

  # P: n_indiv by n_k, phasing matrix
  P <- ifelse(G == 1, 1, 0)
  is_cis <- lapply(snp_pars, function(x) {
    (runif(n_indiv) <= x$prob_ref) * 2 - 1
  }) %>% do.call(what = "cbind")
  P <- P * is_cis

  # calculate G and P for each sample
  G_tilta <- G[origin,]
  P_tilta <- P[origin,]

  return(list(
    Genotype = G_tilta, Phasing = P_tilta,
    Origin = origin, N_indiv = n_indiv
  ))
}


#' Simulate cis-regulatory effect for one-to-multi gene-snp pairs
#'
#' This function simulates genotype G, phasing P, and expression
#' (total and allele-specific) profiles for one-to-multi gene-snp pairs
#'
#' J=1, K>=1, L=1
#'
#' @inheritParams simulateReadCounts.s
#' @inheritParams simulateGenotype.s2m
#' @param baseline Baseline expression (under ref/ref scenario)
#' @param snp_pars A list of length \code{n_k} or \code{1} specifying
#' generative parameters for each regulatory snp. Each element should be a
#' list containing parameters used in \code{\link{simulateGenotype.s2s}},
#' \code{r}, the cis-regulatory effect (log read ratio), and \code{p_error},
#' a vector of phasing error rate of the gene pair in each sample.
#' If \code{length(p_error) == 1} then all samples share the same
#' phasing error rate.
#' If \code{length(snp_pars) == 1} then all \code{n_k} snps share the same set
#' of parameters.
#' @param origin.effect A vector specifying the effect of sample origin with
#' length \code{n_indiv} (the number of samples). If \code{NULL} then no effect
#' will be added
#' @param lib.size A vector of the relative library size of each sample. If
#' \code{NULL} then no adjustment will be made
#' @param confounder A n_i by N matrix where N is the number of confounding factors
#' @param confounder.effect A vector specifying the effect of each confounding factor
#' If \code{NULL} then no effect will be added
#'
#' @export
simulateCisEffect.s2m <- function(n_i, n_k, snp_pars,
                                  phi, prob_as, theta, baseline,
                                  origin = NULL, origin.effect = NULL,
                                  lib.size = NULL,
                                  confounder = NULL, confounder.effect = NULL) {

  if (length(snp_pars) == 1) {
    # all snps share same pars
    snp_pars <- rep(snp_pars, n_k)
  }

  # simulate genotype
  meta <- simulateGenotype.s2m(
    n_i = n_i, n_k = n_k, snp_pars = snp_pars, origin = origin
  )

  # extract regulatory effect
  R <- sapply(snp_pars, function(x) x$r)

  # simulate genetic effects on total read counts
  log_mu <- sapply(1:n_k, function(k){
    ifelse(meta$Genotype[,k] == 1, log(1 + exp(R[k])), 0) +
      ifelse(meta$Genotype[,k] == 2, R[k], 0)
  }) %>% rowSums() + baseline

  # simulate individual effect on total read counts
  if (!is.null(origin.effect)) {
    log_mu <- log_mu + origin.effect[meta$Origin]
  }

  # simulate library size effect on total read counts
  if (!is.null(lib.size)) {
    log_mu <- log_mu + lib.size
  }

  # simulate confounding effect on total read counts
  if (!is.null(confounder.effect)) {
    # check confounder dimensions
    stopifnot(length(confounder.effect) == ncol(confounder))

    log_mu <- log_mu + confounder %*% confounder.effect
  }

  # simulate genetic effects on allelic imbalance
  logit_prob_alt <- meta$Phasing %*% R

  # simulate expression profile
  Y <- simulateReadCounts.s(
    n_i = n_i, mu = exp(log_mu), phi = phi, prob_as = prob_as,
    prob_alt = softmax(logit_prob_alt), theta = theta
  )

  # simulate phasing error for each sample
  P_error <- lapply(snp_pars, function(x) x$p_error)
  for (k in 1:n_k){
    if (is.null(P_error[[k]])) {
      P_error[[k]] <- rep(0, n_i)
    } else {
      if (length(P_error[[k]]) == 1) {
        P_error[[k]] <- rep(P_error[[k]], n_i)
      }
      # reverse P if error occurs
      is_p_error <- runif(n_i) <= P_error[[k]]
      meta$Phasing[is_p_error, k] <- (-1) * meta$Phasing[is_p_error, k]
    }
  }


  sim <- list(
    genotype_pars = list(
      n_i = n_i
    ),
    snp_pars = snp_pars,
    gene_pars = list(
      prob_as = prob_as,
      phi = phi,
      theta = theta,
      baseline = baseline
    ),
    conf_pars = list(
      origin.effect = origin.effect,
      lib.size = lib.size,
      confounder.effect = confounder.effect
    ),
    data = list(
      I = n_i,
      K = n_k,
      N_indiv = meta$N_indiv,
      N = ncol(confounder),
      P = meta$Phasing,
      P_error = P_error,
      G = meta$Genotype,
      T = Y$TotalReadCounts,
      log1p_T = log1p(Y$TotalReadCounts),
      A_ref = Y$RefCounts,
      A_alt = Y$AltCounts,
      Is_ase_het = Y$RefCounts * Y$AltCounts != 0,
      Ori = meta$Origin,
      Cov_ori = outer(meta$Origin, meta$Origin, FUN = "==") * 1,
      Lib_size = lib.size,
      X = confounder
    )
  )

  sim$data$logit_pi_alt <- ifelse(
    !sim$data$Is_ase_het, 0, log(sim$data$A_alt / sim$data$A_ref)
  )

  return(sim)
}
