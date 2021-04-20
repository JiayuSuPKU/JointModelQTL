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
    Origin = origin, N_indiv = n_indiv))
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
#' will be added.
#'
#' @export
simulateCisEffect.s2s <- function(n_i, maf, prob_ref, phi, prob_as,
                                  theta, baseline, r, p_error = NULL,
                                  origin = NULL, origin.effect = NULL) {

  # simulate genotype
  meta <- simulateGenotype.s2s(
    n_i = n_i, maf = maf, prob_ref = prob_ref, origin = origin
  )

  # simulate genetic effects on total read counts
  log_mu <- baseline + ifelse(meta$Genotype == 1, log(1 + exp(r)) - log(2), 0) +
    ifelse(meta$Genotype == 2, r, 0)

  # simulate confounding factors on total read counts
  if (!is.null(origin.effect)) {
    log_mu <- log_mu + origin.effect[meta$Origin]
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
      origin.effect = origin.effect
    ),
    data = list(
      I = n_i,
      N_indiv = meta$N_indiv,
      P = meta$Phasing,
      P_error = p_error,
      G = meta$Genotype,
      T = Y$TotalReadCounts,
      log1p_T = log1p(Y$TotalReadCounts),
      A_ref = Y$RefCounts,
      A_alt = Y$AltCounts,
      Is_ase_het = Y$RefCounts * Y$AltCounts != 0,
      Ori = meta$Origin
    )
  )

  sim$data$logit_pi_alt <- ifelse(
    !sim$data$Is_ase_het, 0, log(sim$data$A_alt / sim$data$A_ref)
  )

  return(sim)
}
