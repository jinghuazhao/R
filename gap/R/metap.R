#' Meta-analysis of p-values with heterogeneity and random effects
#'
#' Performs GWAS-style meta-analysis by combining p-values across studies.
#' Implements Fisher’s method, fixed-effect weighted Stouffer Z,
#' and random-effects Z meta-analysis with heterogeneity statistics.
#'
#' @param data data.frame containing study results.
#' @param N number of studies.
#' @param sided `"two"` (default) or `"one"` for one-sided p-values.
#' @param verbose logical; print summary output.
#' @param prefixp prefix for p-value columns (default `"p"`).
#' @param prefixn prefix for sample-size columns (default `"n"`).
#' @param prefixbeta optional prefix for beta columns (used for direction).
#' @param prefixdir optional prefix for sign columns (+1/-1).
#'
#' @return data.frame with
#' \describe{
#'   \item{fisher_p}{Fisher combined p-value}
#'   \item{stouffer_FE}{Fixed-effect Z meta p-value}
#'   \item{stouffer_RE}{Random-effects Z meta p-value}
#'   \item{Q}{Cochran heterogeneity statistic}
#'   \item{I2}{Proportion heterogeneity}
#'   \item{tau2}{Between-study variance}
#' }
#'
#' @details
#'
#' This function implements the meta-analysis approach used in the
#' Genetic Investigation of ANThropometric Traits (GIANT) consortium.
#'
#' Missing studies are automatically excluded per row, so the number of
#' contributing studies may vary across variants.
#'
#' ## Fisher’s method
#'
#' \deqn{X^2 = -2 \sum_{i=1}^k \log(p_i) \sim \chi^2_{2k}}
#'
#' ## Fixed-effect weighted Stouffer Z
#'
#' Convert p-values to Z:
#' \deqn{z_i = \Phi^{-1}(1 - p_i/2)}
#'
#' Sample-size weights:
#' \deqn{w_i = \sqrt{n_i}}
#'
#' \deqn{Z_{FE} = \frac{\sum w_i z_i}{\sqrt{\sum w_i^2}}}
#'
#' ## Heterogeneity
#'
#' \deqn{Q = \sum w_i (z_i - \bar z)^2}
#' \deqn{I^2 = \max(0, (Q-(k-1))/Q)}
#'
#' ## Random-effects Z meta-analysis
#'
#' DerSimonian–Laird variance:
#' \deqn{\tau^2 = \max\left(0,\frac{Q-(k-1)}{\sum w_i - \sum w_i^2/\sum w_i}\right)}
#'
#' Random-effects weights:
#' \deqn{w_i^* = 1/(1/w_i + \tau^2)}
#'
#' \deqn{Z_{RE} = \frac{\sum w_i^* z_i}{\sqrt{\sum w_i^*}}}
#'
#' @examples
#' \dontrun{
#' ## ------------------------------------------------------------------
#' ## Classic GIANT consortium example (historical demo)
#' ## Speliotes, Elizabeth K., M.D. [ESPELIOTES@PARTNERS.ORG]
#' ## 22-2-2008 MRC-Epid JHZ
#' ## ------------------------------------------------------------------
#'
#' s <- data.frame(
#'   p1 = 0.1^rep(8:2, each = 7),
#'   n1 = rep(32000, 49),
#'   p2 = 0.1^rep(8:2, times = 7),
#'   n2 = rep(8000, 49)
#' )
#'
#' res <- metap(s, 2)
#' head(res)
#'
#' ## Visual comparison equal vs unequal N (original demo)
#'
#' np <- 7
#' p  <- 0.1^((np + 1):2)
#' z  <- qnorm(1 - p/2)
#' n  <- c(32000, 8000)
#'
#' grid <- expand.grid(i = seq_len(np), j = seq_len(np))
#' za <- z[grid$i]; zb <- z[grid$j]
#'
#' metaz_equalN   <- (sqrt(n[1])*za + sqrt(n[1])*zb)/sqrt(n[1]+n[1])
#' metaz_unequalN <- (sqrt(n[1])*za + sqrt(n[2])*zb)/sqrt(n[1]+n[2])
#'
#' q  <- -log10(sort(p, decreasing=TRUE))
#' t1 <- matrix(-log10(sort(2*pnorm(-abs(metaz_equalN))),   decreasing=TRUE), np, np)
#' t2 <- matrix(-log10(sort(2*pnorm(-abs(metaz_unequalN))), decreasing=TRUE), np, np)
#'
#' par(mfrow=c(1,2), bg="white", mar=c(4.2,3.8,0.2,0.2))
#' persp(q,q,t1, main="Equal sample sizes")
#' persp(q,q,t2, main="Unequal sample sizes")
#'
#' ## ------------------------------------------------------------------
#' ## New example: missing studies + heterogeneity
#' ## ------------------------------------------------------------------
#'
#' set.seed(1)
#' dat <- data.frame(
#'   p1 = runif(100,1e-6,0.05),
#'   n1 = sample(5000:20000,100,TRUE),
#'   p2 = runif(100,1e-6,0.05),
#'   n2 = sample(5000:20000,100,TRUE)
#' )
#'
#' ## simulate missing studies
#' dat$p2[sample(1:100,20)] <- NA
#'
#' res <- metap(dat,2)
#' summary(res$I2)
#' }
#'
#' @author ChatGPT (Jing Hua Zhao)
#' @export
#'
metap <- function(
  data, N,
  sided = c("two","one"),
  verbose = TRUE,
  prefixp="p", prefixn="n",
  prefixbeta=NULL, prefixdir=NULL
){
  sided <- match.arg(sided)

  stopifnot(is.data.frame(data), N >= 1)

  p_cols <- paste0(prefixp, seq_len(N))
  n_cols <- paste0(prefixn, seq_len(N))

  if (!all(p_cols %in% names(data)))
    stop("Missing p columns: ", paste(setdiff(p_cols,names(data)),collapse=", "))
  if (!all(n_cols %in% names(data)))
    stop("Missing n columns: ", paste(setdiff(n_cols,names(data)),collapse=", "))

  P    <- as.matrix(data[p_cols])
  Nmat <- as.matrix(data[n_cols])
  M    <- nrow(P)

  ## ---- optional direction ----
  sign_mat <- matrix(1, M, N)

  if (!is.null(prefixbeta))
    sign_mat <- sign(as.matrix(data[paste0(prefixbeta, seq_len(N))]))

  if (!is.null(prefixdir))
    sign_mat <- as.matrix(data[paste0(prefixdir, seq_len(N))])

  ## ---- convert p -> Z ----
  if (sided=="two") Z <- qnorm(1 - P/2)
  if (sided=="one") Z <- qnorm(1 - P)
  Z <- Z * sign_mat

  ## ---- outputs ----
  k_used <- fisher_p <- stouffer_FE <- stouffer_RE <- rep(NA, M)
  Q <- I2 <- tau2 <- rep(NA, M)

  for (j in seq_len(M)) {

    keep <- which(!is.na(P[j,]) & !is.na(Nmat[j,]) & P[j,] > 0 & P[j,] <= 1)
    k <- length(keep)
    if (k == 0) next

    pj <- P[j,keep]
    nj <- Nmat[j,keep]
    zj <- Z[j,keep]

    k_used[j] <- k

    ## =========================================================
    ## Fisher method (now handles missing studies correctly)
    ## =========================================================
    x2 <- -2 * sum(log(pj))
    fisher_p[j] <- pchisq(x2, df = 2*k, lower.tail = FALSE)

    ## =========================================================
    ## Fixed-effect Stouffer
    ## =========================================================
    w <- sqrt(nj)
    z_FE <- sum(w * zj) / sqrt(sum(w^2))
    stouffer_FE[j] <- 2 * pnorm(-abs(z_FE))

    ## =========================================================
    ## Heterogeneity + Random effects
    ## =========================================================
    if (k > 1) {

      wi <- nj                     # inverse-variance approx
      z_bar <- sum(wi*zj) / sum(wi)

      Q[j] <- sum(wi * (zj - z_bar)^2)

      denom <- sum(wi) - sum(wi^2)/sum(wi)
      tau2[j] <- max(0, (Q[j] - (k-1)) / denom)

      wi_star <- 1 / (1/wi + tau2[j])
      z_RE <- sum(wi_star*zj) / sqrt(sum(wi_star))
      stouffer_RE[j] <- 2 * pnorm(-abs(z_RE))

      I2[j] <- max(0, (Q[j] - (k-1)) / Q[j])

    } else {
      ## only one study → FE = RE
      stouffer_RE[j] <- stouffer_FE[j]
      Q[j] <- 0
      I2[j] <- 0
      tau2[j] <- 0
    }
  }

  out <- data.frame(
    k = k_used,
    fisher_p,
    stouffer_FE,
    stouffer_RE,
    Q, I2, tau2
  )

  if (verbose) {
    cat("\nMeta-analysis summary:\n")
    print(head(out))
  }

  return(out)
}
