#' Test and Power Calculation for Mediating Effect
#'
#' This function computes the indirect (mediated) effect in a mediation model
#' and either tests its significance or calculates statistical power.
#' It uses regression coefficients from an independent variable to a mediator (`a`)
#' and from the mediator to an outcome (`b`) with their standard errors.
#' Variance can be estimated using Sobel, Aroian, or Goodman formulas.
#'
#' @param type Character; `"test"` to perform a significance test, `"power"` to estimate power.
#' @param n Numeric; baseline sample size for power calculation (default 25000).
#' @param a Numeric; regression coefficient from independent variable to mediator.
#' @param sa Numeric; standard error of `a`.
#' @param b Numeric; regression coefficient from mediator to outcome.
#' @param sb Numeric; standard error of `b`.
#' @param alpha Numeric; significance level (default 0.05) for power calculation.
#' @param fold Numeric; multiplicative factor for sample size in power calculation.
#' @param method Character; one of `"sobel"`, `"aroian"`, `"goodman"` specifying the variance formula.
#'
#' @details
#' The indirect effect is \eqn{ab}, representing the portion of the effect of \eqn{X} on \eqn{Y}
#' transmitted through mediator \eqn{M}.
#'
#' Standard errors of the indirect effect are approximated using the delta method:
#' - Sobel: \eqn{SE = \sqrt{b^2 s_a^2 + a^2 s_b^2}}
#' - Aroian: \eqn{SE = \sqrt{b^2 s_a^2 + a^2 s_b^2 + s_a^2 s_b^2}}
#' - Goodman: \eqn{SE = \sqrt{b^2 s_a^2 + a^2 s_b^2 - s_a^2 s_b^2}}
#'
#' The z-statistic is \eqn{z = ab / SE(ab)}, assumed to follow a standard normal
#' distribution under \eqn{H_0: ab = 0}.
#'
#' Two-sided p-value: \eqn{p = 2 \Phi(-|z|)}
#' 95% confidence interval: \eqn{ab \pm 1.96 \cdot SE(ab)}
#'
#' Power is based on a non-central chi-square approximation:
#' \eqn{\lambda = z^2 \cdot \text{fold}}
#' \eqn{Power = P(\chi^2_1(\lambda) > \chi^2_{1,1-\alpha})},
#' where \eqn{\chi^2_{1,1-\alpha}} is the critical value at significance level \eqn{\alpha}.
#'
#' These tests assume large samples. For small samples, bootstrap confidence intervals are recommended.
#'
#' @return
#' If `type = "test"`, a named numeric vector with:
#' - `effect`: indirect effect estimate (\eqn{ab})
#' - `se`: standard error
#' - `z`: z-statistic
#' - `p`: two-sided p-value
#' - `lcl`: lower 95% confidence limit
#' - `ucl`: upper 95% confidence limit
#'
#' If `type = "power"`, a named numeric vector with:
#' - `sample_size`: effective sample size (`fold * n`)
#' - `power`: estimated statistical power
#'
#' @references
#' \insertRef{freathy08}{gap}
#'
#' \insertRef{kline05}{gap}
#'
#' \insertRef{mackinnon08}{gap}
#'
#' \insertRef{preacher01}{gap}
#'
#' @examples
#' \dontrun{
#' # Test mediation effect
#' ab(a = 0.15, sa = 0.01,
#'    b = log(1.19), sb = 0.01,
#'    method = "aroian")
#'
#' # Power calculation for 10% increase in sample size
#' ab(type = "power", n = 25000,
#'    a = 0.15, sa = 0.01,
#'    b = log(1.19), sb = 0.01,
#'    fold = 1.1)
#' }
#'
ab <- function(type = c("power", "test"),
               n = 25000,
               a = 0.15,
               sa = 0.01,
               b = log(1.19),
               sb = 0.01,
               alpha = 0.05,
               fold = 1,
               method = c("sobel", "aroian", "goodman"))
{
  type <- match.arg(type)
  method <- match.arg(method)
  # Indirect effect
  effect <- a * b
  # Standard error according to method
  se <- switch(method,
    sobel   = sqrt(b^2 * sa^2 + a^2 * sb^2),
    aroian  = sqrt(b^2 * sa^2 + a^2 * sb^2 + sa^2 * sb^2),
    goodman = sqrt(b^2 * sa^2 + a^2 * sb^2 - sa^2 * sb^2)
  )
  z <- effect / se
  p <- 2 * pnorm(-abs(z))
  lcl <- effect - qnorm(0.975) * se
  ucl <- effect + qnorm(0.975) * se
  if (type == "power") {
    x2 <- z^2
    crit <- qchisq(alpha, df = 1, lower.tail = FALSE)
    power <- pchisq(crit, df = 1, ncp = x2 * fold, lower.tail = FALSE)
    stats <- c(sample_size = fold * n,
               power = power)
  } else {
    stats <- c(effect = effect,
               se = se,
               z = z,
               p = p,
               lcl = lcl,
               ucl = ucl)
  }
  invisible(stats)
}
