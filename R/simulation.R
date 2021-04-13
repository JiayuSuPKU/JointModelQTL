#' Simulate read for one gene
#'
#' This function simulates read counts data from negative binomial (total read counts)
#' and beta-binomial (allele-specific read counts) distributions for one gene
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
simulateReadCountsSingle <- function(n_i, mu, phi, prob_as, prob_alt, theta) {

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
#'
#' @export
simulateGenotypeSingle <- function(n_i, maf, prob_ref) {

  maf <- ifelse(maf <= 0.5, maf, 1 - maf)

  # simulate two loci
  h1 <- rbinom(n = n_i, size = 1, prob = maf)
  h2 <- rbinom(n = n_i, size = 1, prob = maf)
  # G: n_i by 1, genotype matrix
  G <- h1 + h2

  # P: n_i by 1, phasing matrix
  P <- ifelse(G == 1, 1, 0)
  is_cis <- (runif(n_i) <= prob_ref) * 2 - 1
  P <- P * is_cis

  return(list(Genotype = G, Phasing = P))
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
#' @inheritParams simulateReadCountsSingle
#' @inheritParams simulateGenotypeSingle
#' @param baseline Baseline expression (under ref/ref scenario)
#' @param r Cis-regulatory effect defined as the log read ratio under alt/alt
#' and ref/ref scenario
#' @param p_error A vector of phasing error rate of the gene pair in each
#' sample. If \code{length(p_error) == 1} then all samples share the same
#' phasing error rate
#'
#' @export
simulateCisEffectSingle <- function(
  n_i, maf, prob_ref, phi, prob_as, theta, baseline, r, p_error = NULL) {
  ## n_i: number of samples
  ## maf: minor allele frequency of the test snp
  ## prob_ref: probability of the test snp on the ref allele of the target gene
  ##   if the test snp is heterogeneous
  ## T ~ NB(mu, phi)
  ## A_alt ~ BB(T * prob_as = N, prob_alt, theta)
  ## prob_as: proportion of allele-specific reads = N/T
  ## prob_alt: proportion of alt-allele reads in allele-specific reads =
  ##    A_alt / (A_ref + A_alt)

  meta <- simulateGenotypeSingle(n_i = n_i, maf = maf, prob_ref = prob_ref)

  # add genetic effects on total read counts
  log_mu <- baseline + ifelse(meta$Genotype == 1, log(1 + exp(r)) - log(2), 0) +
    ifelse(meta$Genotype == 2, r, 0)

  # add genetic effects on allelic imbalance
  logit_prob_alt <- meta$Phasing * r

  # simulate expression profile
  Y <- simulateReadCountsSingle(
    n_i = n_i, mu = exp(log_mu), phi = phi, prob_as = prob_as,
    prob_alt = softmax(logit_prob_alt), theta = theta
  )

  # consider phasing error for each sample
  if (is.null(p_error)){
    p_error <- rep(0, n_i)
  } else{
    if (length(p_error) == 1)
      p_error <- rep(p_error, n_i)
    # reverse P if error occurs
    is_p_error <- runif(n_i) <= p_error
    meta$Phasing[is_p_error] <- (-1) *meta$Phasing[is_p_error]
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
    data = list(
      I = n_i,
      P = meta$Phasing,
      P_error = p_error,
      G = meta$Genotype,
      T = Y$TotalReadCounts,
      log1p_T = log1p(Y$TotalReadCounts),
      A_ref = Y$RefCounts,
      A_alt = Y$AltCounts,
      Is_ase_het = Y$RefCounts * Y$AltCounts != 0
    )
  )

  sim$data$logit_pi_alt <- ifelse(
    !sim$data$Is_ase_het, 0, log(sim$data$A_alt / sim$data$A_ref)
  )

  return(sim)
}


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
#' used in \code{\link{simulateReadCountsSingle}}.
#' If \code{length(gene_pars) == 1} then all \code{n_j} genes share the same set
#' of parameters.
#'
#' @importFrom extraDistr rbbinom
#' @importFrom stats rbinom rnbinom runif
#'
#' @export
simulateReadCountsMulti <- function(n_i, n_j, gene_pars) {

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
        n = n_i, size = n_as[,j],
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
simulateGenotypeMulti <- function(n_i, n_j, maf, prob_ref) {

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
#' @inheritParams simulateReadCountsMulti
#' @inheritParams simulateGenotypeMulti
#' @param gene_pars A list of length \code{n_j} or \code{1} specifying
#' gene-level statistics. Each element should be a list containing parameters
#' used in \code{\link{simulateCisEffectSingle}}, including phi, prob_as, theta,
#' baseline, and r. If \code{length(gene_pars) == 1} then all \code{n_j} genes
#' share the same set of parameters.
#'
#' @export
simulateCisEffectMulti <- function(n_i, n_j, maf, prob_ref, gene_pars) {

  if (length(gene_pars) == 1) {
    gene_pars <- rep(gene_pars, n_j)
  }

  meta <- simulateGenotypeMulti(n_i = n_i, n_j = n_j, maf = maf, prob_ref = prob_ref)

  gene_pars_mu <- list()

  for (j in 1:n_j) {
    phi <- gene_pars[[j]][["phi"]]
    theta <- gene_pars[[j]][["theta"]]
    prob_as <- gene_pars[[j]][["prob_as"]]
    baseline <- gene_pars[[j]][["baseline"]]
    r <- gene_pars[[j]][["r"]]

    # add genetic effects on total read counts
    log_mu <- baseline + ifelse(meta$Genotype == 1, log(1 + exp(r)) - log(2), 0) +
      ifelse(meta$Genotype == 2, r, 0)

    # add genetic effects on allelic imbalance
    logit_prob_alt <- meta$Phasing[,j] * r

    gene_pars_mu <- append(gene_pars_mu, list(list(
      mu = exp(log_mu), phi = phi, prob_as = prob_as,
      prob_alt = softmax(logit_prob_alt), theta = theta
    )))
  }

  Y <- simulateReadCountsMulti(n_i = n_i, n_j = n_j, gene_pars = gene_pars_mu)

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
