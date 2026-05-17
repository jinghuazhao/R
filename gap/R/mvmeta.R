#' Multivariate fixed-effects meta-analysis via generalized least squares
#'
#' Performs multivariate meta-analysis by pooling study-specific parameter
#' estimates using generalized least squares (GLS) under a fixed-effects model.
#'
#' The function accepts a matrix of parameter estimates and the corresponding
#' within-study covariance matrices (stored in upper-triangular vector form).
#'
#' This approach is appropriate when combining correlated effect estimates,
#' for example correlation coefficients of SNPs across studies.
#'
#' @param b Matrix of study estimates (studies × parameters).
#' @param V may be supplied with
#'  - a matrix of upper-triangular rows
#'  - a list of covariance matrices
#'  - a 3D array (p × p × k)
#'
#' @details
#' The function fits a multivariate **fixed-effects meta-analysis**
#' using generalized least squares (GLS).
#'
#' For study \eqn{i = 1,\dots,k}, let \eqn{d_i} be the vector of observed
#' parameter estimates and \eqn{\Psi_i} the corresponding within-study
#' covariance matrix:
#' \deqn{d_i \sim N(\beta,\; \Psi_i)}
#' where \eqn{\beta} is the vector of common (pooled) parameters.
#'
#' The study estimates are stacked into a single vector
#' \deqn{d = (d_1^T,\dots,d_k^T)^T}
#' with block-diagonal covariance matrix
#' \deqn{\Psi = \mathrm{blockdiag}(\Psi_1,\dots,\Psi_k).}
#'
#' The model can then be written in GLS regression form
#' \deqn{d = X\beta + \varepsilon, \qquad \varepsilon \sim N(0,\Psi)}
#' where \eqn{X} is a block design matrix that repeats an identity matrix
#' for each study (intercept-only multivariate meta-analysis). Missing
#' outcomes are automatically removed when constructing \eqn{d},
#' \eqn{\Psi} and \eqn{X}.
#'
#' The pooled estimator is the GLS estimator
#' \deqn{\hat{\beta} =
#' (X^T \Psi^{-1} X)^{-1} X^T \Psi^{-1} d.}
#'
#' Heterogeneity is assessed using the multivariate Cochran Q statistic
#' \deqn{Q = (d - X\hat{\beta})^T \Psi^{-1} (d - X\hat{\beta}),}
#' which is asymptotically \eqn{\chi^2_{N-p}}, where \eqn{N} is the number
#' of observed estimates and \eqn{p} the number of pooled parameters.
#'
#' This implementation corresponds to the multivariate fixed-effects
#' model described in Hartung et al. (2008, Example 11.3).
#'
#' @importFrom stats AIC BIC confint cov2cor nobs
#' @return An object of class `"mvmeta"` with the following elements:
#' - beta: pooled estimates
#' - vcov: covariance matrix of pooled estimates
#' - se: standard errors
#' - z: z statistics
#' - pval: p values
#' - ci: 95% confidence intervals
#' - X2, df, p: heterogeneity test
#' - logLik: model log-likelihood
#' - k: number of studies
#' - p_outcomes: number of pooled outcomes
#'
#' @references
#' \insertRef{hartung08}{gap}
#'
#' @seealso \code{\link{metareg}}
#'
#' @examples
#' \dontrun{
#' # Example 11.3 from Hartung et al.
#' b <- matrix(c(
#' 0.808, 1.308, 1.379, NA, NA,
#' NA, 1.266, 1.828, 1.962, NA,
#' NA, 1.835, NA, 2.568, NA,
#' NA, 1.272, NA, NA, 2.038,
#' 1.171, 2.024, 2.423, 3.159, NA,
#' 0.681, NA, NA, NA, NA), ncol=5, byrow=TRUE)
#'
#' psi1 <- psi2 <- psi3 <- psi4 <- psi5 <- psi6 <- matrix(0,5,5)
#' psi1[1,1] <- 0.0985; psi1[1,2] <- 0.0611; psi1[1,3] <- 0.0623
#' psi1[2,2] <- 0.1142; psi1[2,3] <- 0.0761; psi1[3,3] <- 0.1215
#'
#' psi2[2,2] <- 0.0713; psi2[2,3] <- 0.0539; psi2[2,4] <- 0.0561
#' psi2[3,3] <- 0.0938; psi2[3,4] <- 0.0698; psi2[4,4] <- 0.0981
#'
#' psi3[2,2] <- 0.1228; psi3[2,4] <- 0.1119; psi3[4,4] <- 0.1790
#' psi4[2,2] <- 0.0562; psi4[2,5] <- 0.0459; psi4[5,5] <- 0.0815
#'
#' psi5[1,1] <- 0.0895; psi5[1,2] <- 0.0729; psi5[1,3] <- 0.0806
#' psi5[1,4] <- 0.0950; psi5[2,2] <- 0.1350; psi5[2,3] <- 0.1151
#' psi5[2,4] <- 0.1394; psi5[3,3] <- 0.1669; psi5[3,4] <- 0.1609
#' psi5[4,4] <- 0.2381
#'
#' psi6[1,1] <- 0.0223
#'
#' V <- rbind(psi1[upper.tri(psi1,diag=TRUE)],
#'            psi2[upper.tri(psi2,diag=TRUE)],
#'            psi3[upper.tri(psi3,diag=TRUE)],
#'            psi4[upper.tri(psi4,diag=TRUE)],
#'            psi5[upper.tri(psi5,diag=TRUE)],
#'            psi6[upper.tri(psi6,diag=TRUE)])
#'
#' fit <- mvmeta(b, V)
#' summary(fit)
#' logLik(fit)
#' AIC(fit)
#' BIC(fit)
#' }
#'
#' @author Jing Hua Zhao
#' @export
#'
mvmeta <- function(b, V)
{
  cl <- match.call()
  if (!is.matrix(b)) stop("b must be a matrix")
  k <- nrow(b); p <- ncol(b)
  ## convert upper-triangle matrix → list of covariance matrices
  if (is.matrix(V)) {
    make_S <- function(vrow) {
      M <- matrix(0,p,p)
      M[upper.tri(M,diag=TRUE)] <- vrow
      M[lower.tri(M)] <- t(M)[lower.tri(M)]
      M
    }
    V <- lapply(seq_len(k), function(i) make_S(V[i,]))
  }
  ## convert 3D array → list
  if (is.array(V) && length(dim(V))==3)
    V <- lapply(seq_len(dim(V)[3]), function(i) V[,,i])

  if (!is.list(V))
    stop("V must be a list or 3D array of covariance matrices")

  ## ensure symmetry
  symmetrise <- function(M) {
    M[lower.tri(M)] <- t(M)[lower.tri(M)]
    M
  }
  V <- lapply(V, symmetrise)
  ## ---- stack outcomes across studies ----
  d <- c(); Psi <- NULL; X <- NULL
  for(i in seq_len(k))
  {
    bi <- b[i, ]
    Si_full <- V[[i]]
    ## observed outcomes = estimate present AND variance present
    obs <- !is.na(bi) & diag(Si_full) > 0
    if(!any(obs)) next
    bi <- bi[obs]
    Si <- Si_full[obs, obs, drop=FALSE]
    d <- c(d, bi)
    Psi <- if(is.null(Psi)) Si else magic::adiag(Psi, Si)
    ## design matrix mapping outcomes → pooled β
    Xi <- matrix(0, length(bi), p)
    Xi[cbind(seq_along(bi), which(obs))] <- 1
    X  <- if(is.null(X)) Xi else rbind(X, Xi)
  }
  ## safe inverse
  inv <- function(M) tryCatch(solve(M), error=function(e) MASS::ginv(M))
  Psi_inv <- inv(Psi)
  XtW  <- crossprod(X, Psi_inv)
  XtWX <- XtW %*% X
  cov.beta <- inv(XtWX)
  beta <- cov.beta %*% XtW %*% d
  ## heterogeneity statistic (multivariate Cochran Q)
  resid <- d - X %*% beta
  X2 <- as.numeric(crossprod(resid, Psi_inv %*% resid))
  df <- length(d) - length(beta)
  pval <- pchisq(X2, df, lower.tail=FALSE)
  if(!is.null(colnames(b))) {
    names(beta) <- colnames(b)
    rownames(cov.beta) <- colnames(cov.beta) <- colnames(b)
  }
  res <- list(
    call = cl, d = d, Psi = Psi, X = X,
    beta = beta, cov.beta = cov.beta,
    X2 = X2, df = df, p = pval
  )
  class(res) <- "mvmeta"
  res
}

#' @export
coef.mvmeta <- function(object, ...) drop(object$beta)

#' @export
vcov.mvmeta <- function(object, ...) object$cov.beta

#' @export
confint.mvmeta <- function(object, parm = NULL, level = 0.95, ...)
{
  se <- sqrt(diag(object$cov.beta))
  z  <- qnorm((1 + level) / 2)
  ci <- cbind(
    lower = coef(object) - z * se,
    upper = coef(object) + z * se
  )
  if (!is.null(names(coef(object))))
    rownames(ci) <- names(coef(object))
  ci
}

#' @export
print.mvmeta <- function(x, digits=4, ...)
{
  cat("Call:\n"); print(x$call)
  cat("\nMultivariate meta-analysis (GLS)\n")
  cat("Number of estimates:", length(x$d), "\n\n")
  cat("Pooled estimates:\n")
  print(round(drop(x$beta), digits))
  cat("\nHeterogeneity test:\n")
  cat("Chi^2 =", round(x$X2,digits),
      " df =", x$df,
      " p =", format.pval(x$p,digits=digits), "\n")
  invisible(x)
}

#' @export
summary.mvmeta <- function(object, ...)
{
  se <- sqrt(diag(object$cov.beta))
  z  <- coef(object)/se
  p  <- 2*pnorm(abs(z), lower.tail=FALSE)
  tab <- cbind(Estimate=coef(object), SE=se, z=z, p=p)
  rownames(tab) <- names(coef(object))

  res <- list(
    call = object$call,
    coefficients = tab,
    confint = confint(object),
    correlation = cov2cor(object$cov.beta),
    heterogeneity = c(X2=object$X2, df=object$df, p=object$p)
  )
  class(res) <- "summary.mvmeta"
  res
}

#' @export
print.summary.mvmeta <- function(x, digits=4, ...)
{
  cat("Call:\n"); print(x$call)
  cat("\nModel summary\n\n")
  print(round(x$coefficients, digits))
  cat("\n95% CI:\n"); print(round(x$confint, digits))
  cat("\nCorrelation of pooled estimates:\n")
  print(round(x$correlation, digits))
  cat("\nHeterogeneity test:\n")
  cat("Chi^2 =", round(x$heterogeneity["X2"],digits),
      " df =", x$heterogeneity["df"],
      " p =", format.pval(x$heterogeneity["p"],digits=digits), "\n")
  invisible(x)
}

#' @export
logLik.mvmeta <- function(object, ...)
{
  r <- object$d - object$X %*% object$beta
  Psi_inv <- tryCatch(solve(object$Psi), error=function(e) MASS::ginv(object$Psi))
  quad <- as.numeric(crossprod(r, Psi_inv %*% r))
  logdet <- as.numeric(determinant(object$Psi, TRUE)$modulus)
  ll <- -0.5*(length(r)*log(2*pi) + logdet + quad)
  attr(ll,"df") <- length(object$beta)
  class(ll) <- "logLik"
  ll
}

#' @export
AIC.mvmeta <- function(object, ...)
{
  ll <- logLik(object); k <- attr(ll,"df")
  -2*as.numeric(ll) + 2*k
}

#' @export
BIC.mvmeta <- function(object, ...)
{
  ll <- logLik(object); k <- attr(ll,"df")
  -2*as.numeric(ll) + log(nobs(object))*k
}

#' @export
nobs.mvmeta <- function(object, ...) length(object$d)
