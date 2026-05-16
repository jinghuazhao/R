#' Fixed and random effects meta-analysis (vectorised implementation)
#'
#' Performs inverse-variance weighted fixed- and random-effects
#' meta-analysis across multiple studies supplied in wide format.
#'
#' The function is designed for high-throughput settings (e.g. GWAS),
#' where each row represents an independent meta-analysis.
#'
#' @param data A data frame containing regression coefficients and
#'   standard errors in wide format.
#' @param N Integer. Number of studies included in the meta-analysis.
#' @param verbose Logical. If `TRUE`, prints a completion message.
#' @param prefixb Character. Prefix for regression coefficients.
#'   Default is `"b"` (columns `b1, b2, …, bN`).
#' @param prefixse Character. Prefix for standard errors.
#'   Default is `"se"` (columns `se1, se2, …, seN`).
#'
#' @details
#' The function accepts wide-format input with estimates
#' \eqn{b_1,\ldots,b_N} and standard errors \eqn{se_1,\ldots,se_N}.
#' Missing values are automatically ignored on a per-row basis.
#'
#' ## Fixed effects model
#' For \eqn{k} studies, inverse-variance weights are defined as
#' \deqn{w_i = 1 / se_i^2}
#'
#' The pooled estimate is
#' \deqn{\beta_f = \frac{\sum_i b_i w_i}{\sum_i w_i}}
#'
#' with standard error
#' \deqn{se_f = \sqrt{1 / \sum_i w_i}}
#'
#' and test statistic
#' \deqn{z_f = \beta_f / se_f}
#'
#' with p-value
#' \deqn{p_f = 2\,\Phi(-|z_f|)}.
#'
#' ## Random effects model (DerSimonian–Laird)
#' Cochran's Q statistic:
#' \deqn{Q = \sum_i w_i (b_i - \beta_f)^2}
#'
#' Between-study variance:
#' \deqn{\tau^2 = \max\left(0, \frac{Q - (k-1)}{\sum w_i - \sum w_i^2 / \sum w_i}\right)}
#'
#' Corrected weights:
#' \deqn{w_i^* = 1 / (1/w_i + \tau^2)}
#'
#' Random-effects pooled estimate:
#' \deqn{\beta_r = \frac{\sum_i b_i w_i^*}{\sum_i w_i^*}}
#'
#' with standard error
#' \deqn{se_r = \sqrt{1/\sum_i w_i^*}}
#'
#' and p-value
#' \deqn{p_r = 2\,\Phi(-|z_r|)}.
#'
#' Heterogeneity p-value:
#' \deqn{p_{heter} = P(\chi^2_{k-1} > Q)}.
#'
#' The heterogeneity statistic is reported as
#' \deqn{I^2 = \max(0, (Q-(k-1))/Q)}.
#'
#' @return A data frame with one row per meta-analysis containing:
#'
#' \describe{
#'   \item{beta_f}{Fixed-effects pooled estimate}
#'   \item{se_f}{Standard error of fixed-effects estimate}
#'   \item{z_f}{Z statistic (fixed effects)}
#'   \item{p_f}{P value (fixed effects)}
#'   \item{beta_r}{Random-effects pooled estimate}
#'   \item{se_r}{Standard error of random-effects estimate}
#'   \item{z_r}{Z statistic (random effects)}
#'   \item{p_r}{P value (random effects)}
#'   \item{p_heter}{Cochran's Q heterogeneity test p-value}
#'   \item{i2}{\eqn{I^2} heterogeneity statistic}
#'   \item{k}{Number of studies contributing to each row}
#'   \item{eps}{Smallest double-precision number used for stability}
#' }
#'
#' @references
#' \insertRef{higgins03}{gap}
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   b1 = 1,  se1 = 2,
#'   b2 = 2,  se2 = 6,
#'   b3 = 3,  se3 = 8
#' )
#' metareg(df, 3)
#'
#' df2 <- data.frame(
#'   b1 = c(1,2), se1 = c(2,4),
#'   b2 = c(2,3), se2 = c(4,6),
#'   b3 = c(3,4), se3 = c(6,8)
#' )
#' metareg(df2, 3)
#' }
#' @author ChatGPT
#'
#' @keywords models meta-analysis
#' @export
#'
metareg <- function(data, N, verbose = FALSE, prefixb = "b", prefixse = "se") {
  stopifnot(is.data.frame(data))
  stopifnot(N >= 2)
  eps <- .Machine$double.eps
  M <- nrow(data)
  # ---- Extract matrices safely ----
  bcols  <- paste0(prefixb,  seq_len(N))
  secols <- paste0(prefixse, seq_len(N))
  if (!all(bcols %in% names(data))) stop("Missing beta columns")
  if (!all(secols %in% names(data))) stop("Missing se columns")
  B  <- as.matrix(data[, bcols, drop = FALSE])
  SE <- as.matrix(data[, secols, drop = FALSE])
  # Valid observations
  OK <- !is.na(B) & !is.na(SE)
  SE[SE < eps] <- eps
  # Fixed effects
  W  <- 1 / SE^2
  W[!OK] <- 0
  SW  <- rowSums(W)
  BW  <- rowSums(W * B)
  SWW <- rowSums(W^2)
  K   <- rowSums(OK)
  beta_f <- BW / SW
  se_f   <- sqrt(1 / SW)
  z_f    <- beta_f / se_f
  p_f    <- 2 * pnorm(-abs(z_f))
  # Cochran Q
  QW <- rowSums(W * (B - beta_f)^2)
  p_heter <- pchisq(QW, K - 1, lower.tail = FALSE)
  # DerSimonian-Laird tau²
  DL <- pmax(0, (QW - (K - 1)) / (SW - SWW / SW))
  DL[!is.finite(DL)] <- 0
  # Random effects
  WC <- 1 / (1/W + DL)
  WC[!OK] <- 0
  SWC <- rowSums(WC)
  beta_r <- rowSums(WC * B) / SWC
  se_r   <- sqrt(1 / SWC)
  z_r    <- beta_r / se_r
  p_r    <- 2 * pnorm(-abs(z_r))
  # I²
  I2 <- pmax(0, (QW - (K - 1)) / QW)
  I2[!is.finite(I2)] <- 0
  if (isTRUE(verbose)) {
    message("Meta-analysis completed for ", M, " rows")
  }
  data.frame(
    beta_f, se_f, z_f, p_f,
    beta_r, se_r, z_r, p_r,
    p_heter, i2 = I2, k = K, eps
  )
}
