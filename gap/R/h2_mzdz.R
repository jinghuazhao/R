#' ACE heritability estimation from MZ/DZ twin correlations
#'
#' Estimates additive genetic (A), shared environmental (C), and
#' unique environmental (E) variance components using Falconer's
#' equations from monozygotic (MZ) and dizygotic (DZ) twin data.
#'
#' @param mzDat Optional data.frame containing MZ twin data.
#' @param dzDat Optional data.frame containing DZ twin data.
#' @param rmz Optional numeric scalar. Correlation for MZ twins.
#' @param rdz Optional numeric scalar. Correlation for DZ twins.
#' @param nmz Optional integer. Sample size for MZ twins.
#' @param ndz Optional integer. Sample size for DZ twins.
#' @param selV Character vector of length 2 specifying twin variables.
#' @param use Correlation missing-data method passed to stats::cor().
#' Default is "complete.obs".
#' @param ci Logical; if TRUE, returns approximate 95% confidence intervals.
#' @param bounds Logical; if TRUE, constrain h2/c2/e2 estimates to \eqn{[0,1]}.
#' @param digits Integer controlling printed rounding.
#'
#' @details
#' Falconer's estimators are:
#'
#' \eqn{h^2 = 2(r_{MZ} - r_{DZ})},
#' \eqn{c^2 = 2r_{DZ} - r_{MZ}}, and
#' \eqn{e^2 = 1 - r_{MZ}}.
#'
#' Approximate variances are computed using
#' \deqn{Var(r) \approx \frac{(1-r^2)^2}{n-1}}
#' Given MZ/DZ data or their correlations and sample sizes, it obtains
#' heritability and variance estimates under an ACE model as in
#' \insertCite{gidziela23;textual}{gap} and \insertCite{keeping95;textual}{gap}.
#'
#' @return
#' A data.frame containing:
#'
#' - `h2`: additive genetic variance estimate.
#' - `c2`: shared environmental variance estimate.
#' - `e2`: unique environmental variance estimate.
#' - `vh2`: estimated variance of `h2`.
#' - `vc2`: estimated variance of `c2`.
#' - `ve2`: estimated variance of `e2`.
#' - `se_h2`: standard error of `h2`.
#' - `se_c2`: standard error of `c2`.
#' - `se_e2`: standard error of `e2`.
#' - `lcl_h2`: lower 95% confidence limit for `h2`.
#' - `ucl_h2`: upper 95% confidence limit for `h2`.
#' - `lcl_c2`: lower 95% confidence limit for `c2`.
#' - `ucl_c2`: upper 95% confidence limit for `c2`.
#' - `lcl_e2`: lower 95% confidence limit for `e2`.
#' - `ucl_e2`: upper 95% confidence limit for `e2`.
#'
#' @examples
#' \dontrun{
#' h2_mzdz(rmz = 0.72,
#'         rdz = 0.33,
#'         nmz = 384,
#'         ndz = 243)
#' library(mvtnorm)
#' set.seed(12345)
#' selVars <- c("bmi1", "bmi2")
#' mzm <- as.data.frame(
#'   mvtnorm::rmvnorm(
#'     195,
#'     c(22.75, 22.75),
#'     matrix(2.66^2 * c(1, 0.67, 0.67, 1), 2)
#'   )
#' )
#' dzm <- as.data.frame(
#'   mvtnorm::rmvnorm(
#'     130,
#'     c(23.44, 23.44),
#'     matrix(2.75^2 * c(1, 0.32, 0.32, 1), 2)
#'   )
#' )
#' names(mzm) <- selVars
#' names(dzm) <- selVars
#' res <- ACE_CI(
#'   mzData = mzm,
#'   dzData = dzm,
#'   selV = selVars,
#'   n.boot = 50,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' print(res)
#' }
#'
#' @export
#' @references
#' \insertRef{elks12}{gap}
#'
#' \insertAllCited{}
#'
h2_mzdz <- function(
    mzDat = NULL,
    dzDat = NULL,
    rmz = NULL,
    rdz = NULL,
    nmz = NULL,
    ndz = NULL,
    selV = NULL,
    use = "complete.obs",
    ci = TRUE,
    bounds = FALSE,
    digits = 3
) {

  #-----------------------------
  # Internal helper
  #-----------------------------
  get_corr_info <- function(dat, r, n, selV, use, group) {

    if (!is.null(dat)) {

      if (is.null(selV) || length(selV) != 2) {
        stop("selV must contain exactly two variable names.")
      }

      if (!all(selV %in% names(dat))) {
        missing_vars <- selV[!selV %in% names(dat)]
        stop(sprintf(
          "%s data missing variable(s): %s",
          group,
          paste(missing_vars, collapse = ", ")
        ))
      }

      complete_idx <- stats::complete.cases(dat[, selV])
      n_complete <- sum(complete_idx)

      if (n_complete < 4) {
        stop(sprintf(
          "%s data has fewer than 4 complete pairs.",
          group
        ))
      }

      r_est <- stats::cor(
        dat[[selV[1]]],
        dat[[selV[2]]],
        use = use
      )

      if (!is.finite(r_est)) {
        stop(sprintf("Unable to compute %s correlation.", group))
      }

      return(list(r = unname(r_est), n = n_complete))
    }

    # Correlation/sample-size input path
    if (is.null(r) || is.null(n)) {
      stop(sprintf(
        "%s requires either raw data or both correlation and sample size.",
        group
      ))
    }

    if (!is.numeric(r) || length(r) != 1 || abs(r) > 1) {
      stop(sprintf("%s correlation must be a numeric scalar in [-1, 1].", group))
    }

    if (!is.numeric(n) || length(n) != 1 || n < 4) {
      stop(sprintf("%s sample size must be >= 4.", group))
    }

    list(r = r, n = as.integer(n))
  }

  #-----------------------------
  # Extract correlations and n
  #-----------------------------
  mz <- get_corr_info(
    dat = mzDat,
    r = rmz,
    n = nmz,
    selV = selV,
    use = use,
    group = "MZ"
  )

  dz <- get_corr_info(
    dat = dzDat,
    r = rdz,
    n = ndz,
    selV = selV,
    use = use,
    group = "DZ"
  )

  r1 <- mz$r
  r2 <- dz$r
  n1 <- mz$n
  n2 <- dz$n

  #-----------------------------
  # Falconer ACE estimates
  #-----------------------------
  h2 <- 2 * (r1 - r2)
  c2 <- 2 * r2 - r1
  e2 <- 1 - r1

  if (bounds) {
    h2 <- max(0, min(1, h2))
    c2 <- max(0, min(1, c2))
    e2 <- max(0, min(1, e2))
  }

  #-----------------------------
  # Variance estimates
  #-----------------------------
  vmz <- ((1 - r1^2)^2) / (n1 - 1)
  vdz <- ((1 - r2^2)^2) / (n2 - 1)

  vh2 <- 4 * (vmz + vdz)
  vc2 <- vmz + 4 * vdz
  ve2 <- vmz

  se_h2 <- sqrt(vh2)
  se_c2 <- sqrt(vc2)
  se_e2 <- sqrt(ve2)

  out <- data.frame(
    rmz = r1,
    rdz = r2,
    nmz = n1,
    ndz = n2,
    h2 = h2,
    c2 = c2,
    e2 = e2,
    vh2 = vh2,
    vc2 = vc2,
    ve2 = ve2,
    se_h2 = se_h2,
    se_c2 = se_c2,
    se_e2 = se_e2,
    row.names = NULL
  )

  #-----------------------------
  # Confidence intervals
  #-----------------------------
  if (ci) {

    zcrit <- stats::qnorm(0.975)

    out$lcl_h2 <- h2 - zcrit * se_h2
    out$ucl_h2 <- h2 + zcrit * se_h2

    out$lcl_c2 <- c2 - zcrit * se_c2
    out$ucl_c2 <- c2 + zcrit * se_c2

    out$lcl_e2 <- e2 - zcrit * se_e2
    out$ucl_e2 <- e2 + zcrit * se_e2

    if (bounds) {
      ci_cols <- grep("^(lcl|ucl)_", names(out))
      out[ci_cols] <- lapply(out[ci_cols], function(x) {
        pmax(0, pmin(1, x))
      })
    }
  }

  #-----------------------------
  # Nicely rounded print method
  #-----------------------------
  numeric_cols <- vapply(out, is.numeric, logical(1))
  out[numeric_cols] <- lapply(
    out[numeric_cols],
    round,
    digits = digits
  )

  return(out)
}

#' Bootstrap ACE confidence intervals
#'
#' Estimates ACE components (additive genetic, shared environment,
#' unique environment) using nonparametric bootstrap resampling of
#' twin data (MZ and DZ pairs).
#'
#' @param mzData Data frame for monozygotic twins.
#' @param dzData Data frame for dizygotic twins.
#' @param selV Character vector of length 2 specifying twin variables.
#' @param n.boot Number of bootstrap replications.
#' @param conf.level Confidence level for percentile intervals.
#' @param use Missing-data handling method passed to cor().
#' @param bounds Logical; constrain ACE estimates to \eqn{[0,1]}.
#' @param seed Random seed for reproducibility.
#' @param verbose Logical; show progress messages.
#'
#' @return A list containing:
#' \itemize{
#'   \item observed: ACE estimates from original data
#'   \item bootstrap: bootstrap replicate estimates
#'   \item summary: bootstrap summary statistics
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' res <- ACE_CI(mzData = mzm, dzData = dzm, selV = c("bmi1","bmi2"))
#' print(res)
#' }
#'
#' @export
ACE_CI <- function(
    mzData,
    dzData,
    selV,
    n.boot = 1000,
    conf.level = 0.95,
    use = "complete.obs",
    bounds = FALSE,
    seed = NULL,
    verbose = interactive()
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  observed <- h2_mzdz(
    mzDat = mzData,
    dzDat = dzData,
    selV = selV,
    use = use,
    ci = TRUE,
    bounds = bounds
  )

  nmz <- nrow(mzData)
  ndz <- nrow(dzData)

  boot_res <- vector("list", n.boot)

  for (i in seq_len(n.boot)) {

    if (verbose && (i %% max(1, floor(n.boot / 20)) == 0)) {
      message(sprintf("Bootstrap %d / %d", i, n.boot))
    }

    mz_idx <- sample.int(nmz, replace = TRUE)
    dz_idx <- sample.int(ndz, replace = TRUE)

    mz_boot <- mzData[mz_idx, , drop = FALSE]
    dz_boot <- dzData[dz_idx, , drop = FALSE]

    boot_res[[i]] <- tryCatch(
      h2_mzdz(
        mzDat = mz_boot,
        dzDat = dz_boot,
        selV = selV,
        use = use,
        ci = FALSE,
        bounds = bounds,
        digits = 10
      ),
      error = function(e) NULL
    )
  }

  boot_df <- do.call(rbind, boot_res)

  if (is.null(boot_df) || nrow(boot_df) == 0) {
    stop("All bootstrap replicates failed.")
  }

  alpha <- 1 - conf.level
  params <- c("h2", "c2", "e2")

  summary_df <- do.call(
    rbind,
    lapply(params, function(p) {

      x <- boot_df[[p]]
      x <- x[is.finite(x)]

      data.frame(
        parameter = p,
        mean = mean(x),
        sd = stats::sd(x),
        median = stats::median(x),
        lcl = stats::quantile(x, alpha / 2),
        ucl = stats::quantile(x, 1 - alpha / 2),
        row.names = NULL
      )
    })
  )

  out <- list(
    observed = observed,
    bootstrap = boot_df,
    summary = summary_df
  )

  class(out) <- "ACE_CI"

  return(out)
}

#' @noRd
#' @export
print.ACE_CI <- function(x, digits = 3, ...) {

  cat("\nObserved ACE estimates\n")
  print(round(x$observed, digits))

  cat("\nBootstrap summary\n")

  sum_df <- x$summary
  num_cols <- vapply(sum_df, is.numeric, logical(1))
  sum_df[num_cols] <- lapply(sum_df[num_cols], round, digits)

  print(sum_df, row.names = FALSE)

  invisible(x)
}
