#' Hardy-Weinberg Equilibrium Test (Multiallelic, Unified Interface)
#'
#' Unified Hardy-Weinberg equilibrium (HWE) testing procedure for multiallelic loci.
#' All input formats are internally converted to a single genotype count matrix,
#' ensuring identical inference across representations.
#'
#' @param data genotype data in one of three formats:
#'  - alleles two-column allele pairs
#'  - genotypes integer genotype IDs (triangular encoding)
#'  - counts symmetric genotype count matrix
#'
#' @param type input format: "alleles", "genotypes", or "counts"
#' @param verbose logical; if TRUE, prints full test output
#' @param yates.correct logical; if TRUE applies Yates continuity correction
#'   to Pearson chi-square statistic only
#' @param B number of Monte Carlo replicates for exact test
#' @param seed random seed for Monte Carlo exact test
#'
#' @details
#' Let allele frequencies be:
#' \deqn{p_i = \frac{c_i}{2n}}
#'
#' Under Hardy–Weinberg equilibrium:
#' \deqn{
#' P_{ii} = p_i^2,\quad P_{ij} = 2p_i p_j \ (i \ne j)
#' }
#'
#' Expected genotype counts:
#' \deqn{E_{ij} = n P_{ij}}
#'
#' All input formats are mapped to the same genotype matrix,
#' ensuring algebraic equivalence.
#'
#' \emph{Pearson chi-square}
#' \deqn{
#' X^2 = \sum_{i \le j} \frac{(O_{ij} - E_{ij})^2}{E_{ij}}
#' }
#'
#' With optional Yates correction:
#' \deqn{
#' X^2 = \sum \frac{(|O_{ij} - E_{ij}| - 0.5)^2}{E_{ij}}
#' }
#'
#' \emph{Likelihood ratio test}
#' \deqn{
#' G^2 = 2 \sum_{i \le j} O_{ij} \log(O_{ij}/E_{ij})
#' }
#'
#' with convention 0log(0) = 0.
#'
#' \emph{Inbreeding coefficient}
#' \deqn{
#' F = \frac{H_{obs} - H_{exp}}{1 - H_{exp}}
#' }
#' where:
#' \deqn{
#' H_{obs} = \sum_i O_{ii}/n,\quad H_{exp} = \sum_i p_i^2
#' }
#'
#' Under:
#' \deqn{X \sim \text{Multinomial}(n, P)}
#'
#' the p-value is:
#' \deqn{
#' p = \Pr(\ell(X^{sim}) \le \ell(X^{obs}))
#' }
#' where:
#' \deqn{
#' \ell(x) = \sum x_i \log p_i + \log \frac{n!}{\prod x_i!}
#' }
#'
#' \deqn{
#' \text{alleles} \equiv \text{genotype IDs} \equiv \text{count matrix}
#' \rightarrow M_{ij}
#' }
#'
#' Hence all statistics are identical up to Monte Carlo variation.
#'
#' @note Note that
#' - Zero-frequency alleles are removed automatically
#' - Exact test is Monte Carlo (fixed seed for reproducibility)
#' - All statistics are computed on the same standardized genotype matrix
#'
#' @return Named list containing:
#' - `source` — input representation used ("alleles", "genotypes", or "counts")
#' - `X2` — Pearson chi-square statistic
#' - `p_X2` — asymptotic p-value of Pearson chi-square test
#' - `LRT` — likelihood ratio statistic
#' - `p_LRT` — asymptotic p-value of likelihood ratio test
#' - `p_exact` — Monte Carlo exact test p-value
#' - `freq` — allele frequency vector
#' - `rho` — inbreeding coefficient
#'
#' @examples
#' \dontrun{
#' a1 <- c(1,1,1,1,2,2,2,3,3,1,2,3,1,2,3,1,2,3)
#' a2 <- c(1,2,3,1,2,3,2,3,1,2,3,1,1,1,2,3,2,3)
#' r1 <- hwe(cbind(a1,a2), "alleles", FALSE)
#' g <- a2g(a1,a2)
#' r2 <- hwe(g, "genotypes", FALSE)
#' g_tab <- table(g)
#' pairs <- g2a(as.integer(names(g_tab)))
#' k <- max(pairs)
#' M <- matrix(0, k, k)
#' for(i in seq_len(nrow(pairs))) {
#'   M[pairs[i,1], pairs[i,2]] <- g_tab[i]
#' }
#' r3 <- hwe(M, "counts", FALSE)
#' r <- lapply(list(r1, r2, r3), \(x) within(as.data.frame(x), {
#'      freq <- paste(round(as.numeric(freq), 3), collapse=";")})[1,])
#' do.call(rbind,r)
#' }
#' @export
#'
hwe <- function(data,
                type=c("alleles","genotypes","counts"),
                verbose=TRUE,
                yates.correct=FALSE,
                B=1e5,
                seed=123){
  type <- match.arg(type)
  # Monte Carlo exact test
  hwe_exact <- function(obs, probs, n){
    logMultinom <- function(x, p){
      lgamma(sum(x)+1) - sum(lgamma(x+1)) + sum(x*log(p))
    }
    set.seed(seed)
    obs_ll <- logMultinom(obs, probs)
    sim_ll <- numeric(B)
    for(b in seq_len(B)){
      sim <- as.vector(rmultinom(1,n,probs))
      sim_ll[b] <- logMultinom(sim, probs)
    }
    mean(sim_ll <= obs_ll)
  }
  # INPUT NORMALISATION
  if(type=="alleles"){
    if(ncol(data)!=2) stop("Allele input must have 2 columns")
    g <- a2g(data[,1], data[,2])
    gtab <- table(g)
    pairs <- g2a(as.integer(names(gtab)))
    k <- max(pairs)
    M <- matrix(0,k,k)
    for(i in seq_len(nrow(pairs)))
      M[pairs[i,1],pairs[i,2]] <- gtab[i]
  }
  if(type=="genotypes"){
    gtab <- table(data)
    pairs <- g2a(as.integer(names(gtab)))
    k <- max(pairs)
    M <- matrix(0,k,k)
    for(i in seq_len(nrow(pairs)))
      M[pairs[i,1],pairs[i,2]] <- gtab[i]
  }
  if(type=="counts"){
    M <- as.matrix(data)
    if(nrow(M)!=ncol(M))
      stop("Genotype count matrix must be square")
  }
  # CORE STATISTICS
  n <- sum(M)
  k <- nrow(M)
  allele.counts <- numeric(k)
  for(i in seq_len(k)){
    allele.counts[i] <- 2*M[i,i] +
      sum(M[i,-i]) +
      sum(M[-i,i])
  }
  p <- allele.counts/(2*n)
  keep <- p > 0
  M <- M[keep,keep,drop=FALSE]
  p <- p[keep]
  k <- length(p)
  exp <- matrix(0,k,k)
  prob_mat <- matrix(0,k,k)
  for(i in seq_len(k)){
    for(j in i:k){
      exp[i,j] <- ifelse(i==j, n*p[i]^2, 2*n*p[i]*p[j])
      prob_mat[i,j] <- ifelse(i==j, p[i]^2, 2*p[i]*p[j])
    }
  }
  obs <- M[upper.tri(M,diag=TRUE)]
  expv <- exp[upper.tri(exp,diag=TRUE)]
  probs <- prob_mat[upper.tri(prob_mat,diag=TRUE)]
  # Pearson chi-square
  if(yates.correct){
    X2 <- sum((abs(obs-expv)-0.5)^2/expv)
  } else {
    X2 <- sum((obs-expv)^2/expv)
  }
  df <- length(obs)-length(p)
  p_chisq <- pchisq(X2,df,lower.tail=FALSE)
  # Likelihood ratio test
  LRT <- 2*sum(ifelse(obs>0, obs*log(obs/expv), 0))
  p_LRT <- pchisq(LRT,df,lower.tail=FALSE)
  # Inbreeding coefficient
  Hobs <- sum(diag(M))/n
  Hexp <- sum(p^2)
  rho <- (Hobs-Hexp)/(1-Hexp)
  # Exact test
  p_exact <- hwe_exact(obs, probs, n)
  if(verbose){
    cat("Hardy-Weinberg equilibrium test (multiallelic)\n\n")
    cat("Pearson chi-square:\n")
    cat(" X2 =",X2," df =",df," p =",p_chisq,"\n\n")
    cat("Likelihood ratio:\n")
    cat(" LRT =",LRT," df =",df," p =",p_LRT,"\n\n")
    cat("Exact test (MC):\n")
    cat(" p =",p_exact,"\n\n")
    cat("Allele frequencies:\n")
    print(round(p,4))
    cat("\nrho =",rho,"\n")
  }
  return(list(
    source = type,
    X2 = X2,
    p_X2 = p_chisq,
    LRT = LRT,
    p_LRT = p_LRT,
    p_exact = p_exact,
    freq = p,
    rho = rho
  ))
}

# 08-02-2004 Working but miss.value needs to be added later
# 09-02-2004 Add using genotype counts
# 11-02-2004 Works ok in all three modes
# 17-02-2004 Add is.miss function
# 04-10-2004 tidy is.count/is.genotype with data.type and document allele.freq
# 06-11-2004 simply code
# 08-11-2004 debug, refer to code of August 2004 and done
# 09-11-2004 correct allele lables for allele frequencies
# 10-11-2004 fix lrt statistic and p.lrt
# 16-01-2005 fix n.genotype
# 16-05-2026 Rewrite via vibe coding
