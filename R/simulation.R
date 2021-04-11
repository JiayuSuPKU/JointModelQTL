library(extraDistr)

# J = 1, K = 1, L = 1

simulateReadCountsSingle <- function(n_i, mu, phi, prob_as, prob_alt, theta) {
  ## n_i: number of samples
  ## T ~ NB(mu, phi)
  ## A_alt ~ BB(T * prob_as = N, prob_alt, theta)
  ## prob_as: proportion of allele-specific reads = N/T
  ## prob_alt: proportion of alt-allele reads in allele-specific reads =
  ##    A_alt / (A_ref + A_alt)

  trc <- rnbinom(n = n_i, mu = mu, size = phi)
  n_as <- round(trc * prob_as)
  a_alt <- rbbinom(
    n = n_i, size = n_as,
    alpha = prob_alt * (1 / theta - 1),
    beta = (1 - prob_alt) * (1 / theta - 1)
  )

  return(list(TotalReadCounts = trc, RefCounts = n_as - a_alt, AltCounts = a_alt))
}

simulateGenotypeSingle <- function(n_i, maf, prob_ref) {
  ## n_i: number of samples
  ## maf: minor allele frequency of the test snp
  ## prob_ref: probability of the test snp on the ref allele of the target gene
  ##   if the test snp is heterogeneous

  maf <- ifelse(maf <= 0.5, maf, 1 - maf)
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

logistic <- function(x) {
  return(exp(x) / (1 + exp(x)))
}

simulateCisEffectSingle <- function(n_i, maf, prob_ref, phi, prob_as, theta, baseline, r) {
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
    prob_alt = logistic(logit_prob_alt), theta = theta
  )

  sim <- list(
    pars = list(
      n_i = n_i,
      maf = maf,
      prob_ref = prob_ref,
      prob_as = prob_as,
      phi = phi,
      theta = theta,
      baseline = baseline,
      r = r
    ),
    data = list(
      I = n_i,
      P = meta$Phasing,
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

simulateReadCountsMulti <- function(n_i, n_j, gene_pars) {
  ## n_i: number of samples
  ## n_j: number of genes
  ## gene_pars = list(list(mu, phi, prob_as, prob_alt, theta), ...)
  ## T ~ NB(mu, phi)
  ## A_alt ~ BB(T * prob_as = N, prob_alt, theta)
  ## prob_as: proportion of allele-specific reads = N/T
  ## prob_alt: proportion of alt-allele reads in allele-specific reads =
  ##    A_alt / (A_ref + A_alt)

  if (length(gene_pars) == 1) {
    gene_pars <- rep(gene_pars, n_j)
  }

  trc <- sapply(gene_pars, function(l) {
    return(rnbinom(n = n_i, mu = l[["mu"]], size = l[["phi"]]))
  })

  prob_as <- sapply(gene_pars, function(l) {
    return(l[["prob_as"]])
  })
  n_as <- round(t(t(trc) * prob_as))

  a_alt <- sapply(1:n_j, function(j) {
    prob_alt <- gene_pars[[j]][["prob_alt"]]
    theta <- gene_pars[[j]][["theta"]]
    return(
      rbbinom(
        n = n_i, size = n_as[,j],
        alpha = prob_alt * (1 / theta - 1),
        beta = (1 - prob_alt) * (1 / theta - 1)
      )
    )
  })


  return(list(TotalReadCounts = trc, RefCounts = n_as - a_alt, AltCounts = a_alt))
}

simulateGenotypeMulti <- function(n_i, n_j, maf, prob_ref) {
  ## n_i: number of samples
  ## n_j: number of genes
  ## maf: minor allele frequency of the test snp
  ## prob_ref: probability of the test snp on the ref allele of the target gene
  ##   if the test snp is heterogeneous

  if (length(prob_ref) == 1) {
    prob_ref <- rep(prob_ref, n_j)
  }

  maf <- ifelse(maf <= 0.5, maf, 1 - maf)
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

simulateCisEffectMulti <- function(n_i, n_j, maf, prob_ref, gene_pars) {
  ## n_i: number of samples
  ## maf: minor allele frequency of the test snp
  ## prob_ref: probability of the test snp on the ref allele of the target gene
  ##   if the test snp is heterogeneous
  ## gene_pars = list(list(phi, prob_as, theta, baseline, r), ...)
  ## T ~ NB(mu, phi)
  ## A_alt ~ BB(T * prob_as = N, prob_alt, theta)
  ## prob_as: proportion of allele-specific reads = N/T
  ## prob_alt: proportion of alt-allele reads in allele-specific reads =
  ##    A_alt / (A_ref + A_alt)

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
    logit_prob_alt <- meta$Phasing * r

    gene_pars_mu <- append(gene_pars_mu, list(list(
      mu = exp(log_mu), phi = phi, prob_as = prob_as,
      prob_alt = logistic(logit_prob_alt), theta = theta
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
      Is_ase_het = Y$RefCounts * Y$AltCounts != 0
    )
  )
  sim$data$logit_pi_alt <- ifelse(
    !sim$data$Is_ase_het, 0, log(sim$data$A_alt / sim$data$A_ref)
  )

  return(sim)
}
