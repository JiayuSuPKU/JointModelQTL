#' Simulate read for multiple genes
#'
#' This function simulates read counts data from negative binomial (total read counts)
#' and beta-binomial (allele-specific read counts) distributions for multiple genes
#'
#' J>=1, K=1, L=1
#'
#' \code{T ~ NB(mu, phi)}
#' \code{A_alt ~ BB(T * prob_as := N_as, prob_alt, theta)}
#'
#' @param n_i Numbet of samples
#' @param n_j Number of genes
#' @param gene_pars A list of length \code{n_j} or \code{1} specifying
#' gene-level statistics. Each element should be a list containing parameters
#' used in \code{\link{simulateReadCounts.s}}.
#' If \code{length(gene_pars) == 1} then all \code{n_j} genes share the same set
#' of parameters.
#'
#' @importFrom extraDistr rbbinom
#' @importFrom stats rbinom rnbinom runif
#'
#' @export
simulateReadCounts.m <- function(n_i, n_j, gene_pars) {
  if (length(gene_pars) == 1) {
    gene_pars <- rep(gene_pars, n_j)
  }

  # total read counts
  trc <- sapply(gene_pars, function(l) {
    return(rnbinom(n = n_i, mu = l[["mu"]], size = l[["phi"]]))
  })

  # allele-specific read counts
  prob_as <- sapply(gene_pars, function(l) {
    return(l[["prob_as"]])
  })
  n_as <- round(t(t(trc) * prob_as))

  # alternative allele read counts
  a_alt <- sapply(1:n_j, function(j) {
    prob_alt <- gene_pars[[j]][["prob_alt"]]
    theta <- gene_pars[[j]][["theta"]]
    return(
      rbbinom(
        n = n_i, size = n_as[, j],
        alpha = prob_alt * (theta - 1),
        beta = (1 - prob_alt) * (theta - 1)
      )
    )
  })

  return(list(TotalReadCounts = trc, RefCounts = n_as - a_alt, AltCounts = a_alt))
}


#' Simulate genotype for multi-to-one gene-snp pairs
#'
#' This function simulates genotype G and phasing P for multi-to-one gene-snp pairs
#'
#' J>=1, K=1, L=1
#'
#' @param n_i Number of samples
#' @param n_j Number of genes
#' @param maf Minor allele frequency of the test snp
#' @param prob_ref A vector specifying the probability of the test snp on the
#' ref allele of each gene if the test snp is heterogeneous.
#' If \code{length(prob_ref) == 1} then all genes share the same \code{prob_ref}
#'
#' @export
simulateGenotype.m2s <- function(n_i, n_j, maf, prob_ref) {
  if (length(prob_ref) == 1) {
    prob_ref <- rep(prob_ref, n_j)
  }

  maf <- ifelse(maf <= 0.5, maf, 1 - maf)

  # simulate two loci
  h1 <- rbinom(n = n_i, size = 1, prob = maf)
  h2 <- rbinom(n = n_i, size = 1, prob = maf)

  # G: n_i by 1, genotype matrix
  G <- h1 + h2

  # P: n_i by n_j, phasing matrix
  P <- sapply(prob_ref, function(x) {
    is_cis <- (runif(n_i) <= x) * 2 - 1
    return(ifelse(G == 1, 1, 0) * is_cis)
  })

  return(list(Genotype = G, Phasing = P))
}


#' Simulate cis-regulatory effect for multi-to-one gene-snp pairs
#'
#' This function simulates genotype G, phasing P, and expression
#' (total and allele-specific) profiles for multi-to-one gene-snp pairs
#'
#' J>=1, K=1, L=1
#'
#' @inheritParams simulateReadCounts.m
#' @inheritParams simulateGenotype.m2s
#' @param gene_pars A list of length \code{n_j} or \code{1} specifying
#' gene-level statistics. Each element should be a list containing parameters
#' used in \code{\link{simulateCisEffect.s2s}}, including phi, prob_as, theta,
#' baseline, and r. If \code{length(gene_pars) == 1} then all \code{n_j} genes
#' share the same set of parameters.
#'
#' @export
simulateCisEffect.m2s <- function(n_i, n_j, maf, prob_ref, gene_pars) {
  if (length(gene_pars) == 1) {
    gene_pars <- rep(gene_pars, n_j)
  }

  # simulate genotype
  meta <- simulateGenotype.m2s(n_i = n_i, n_j = n_j, maf = maf, prob_ref = prob_ref)

  gene_pars_mu <- list()

  for (j in 1:n_j) {
    phi <- gene_pars[[j]][["phi"]]
    theta <- gene_pars[[j]][["theta"]]
    prob_as <- gene_pars[[j]][["prob_as"]]
    baseline <- gene_pars[[j]][["baseline"]]
    r <- gene_pars[[j]][["r"]]

    # simulate genetic effects on total read counts
    log_mu <- baseline + ifelse(meta$Genotype == 1, log(1 + exp(r)) - log(2), 0) +
      ifelse(meta$Genotype == 2, r, 0)

    # simulate genetic effects on allelic imbalance
    logit_prob_alt <- meta$Phasing[, j] * r

    gene_pars_mu <- append(gene_pars_mu, list(list(
      mu = exp(log_mu), phi = phi, prob_as = prob_as,
      prob_alt = softmax(logit_prob_alt), theta = theta
    )))
  }

  # simulate expression profiles
  Y <- simulateReadCounts.m(n_i = n_i, n_j = n_j, gene_pars = gene_pars_mu)

  sim <- list(
    genotype_pars = list(
      n_i = n_i,
      maf = maf,
      prob_ref = prob_ref
    ),
    gene_pars = gene_pars,
    data = list(
      I = n_i,
      J = n_j,
      P = meta$Phasing,
      G = meta$Genotype,
      T = Y$TotalReadCounts,
      log1p_T = log1p(Y$TotalReadCounts),
      A_ref = Y$RefCounts,
      A_alt = Y$AltCounts,
      Is_ase_het = as.matrix(Y$RefCounts * Y$AltCounts != 0) * 1
    )
  )
  sim$data$logit_pi_alt <- ifelse(
    !sim$data$Is_ase_het, 0, log(sim$data$A_alt / sim$data$A_ref)
  )

  return(sim)
}
