#' Power and sample size for case-cohort design
#'
#' @param n the total number of subjects in the cohort.
#' @param q the sampling fraction of the subcohort.
#' @param pD the proportion of the failures in the full cohort.
#' @param p1 proportions of the two groups (p2=1-p1).
#' @param theta log-hazard ratio for two groups.
#' @param alpha type I error -- significant level.
#' @param beta type II error.
#' @param power if specified, the power for which sample size is calculated.
#' @param verbose error messages are explicitly printed out.
#'
#' @details
#' The power of the test is according to 
#' \deqn{\Phi\left(Z_\alpha+m^{1/2}\theta\sqrt{\frac{p_1p_2p_D}{q+(1-q)p_D}}\right)}{Phi(Z_alpha+m^0.5*theta*sqrt(p_1p_2p_D/q+(1-q)p_D))}
#' where \eqn{\alpha}{alpha} is the significance level, \eqn{\theta}{theta} is the log-hazard ratio for two groups, \eqn{p_j}{p_j}, 
#' j=1, 2, are the proportion of the two groups in the population. \eqn{m} is the total number of subjects in the subcohort, 
#' \eqn{p_D} is the proportion of the failures in the full cohort, and \eqn{q} is the sampling fraction of the subcohort.
#'
#' Alternatively, the sample size required for the subcohort is \deqn{m=nBp_D/(n-B(1-p_D))}{m=nBp_D/(n-B(1-p_D))}
#' where \eqn{B=(Z_{1-\alpha}+Z_\beta)^2/(\theta^2p_1p_2p_D)}{B=(Z_{1-alpha}+Z_beta)^2/(theta^2p_1p_2p_D))}, and \eqn{n} is the size of cohort.
#'
#' When infeaisble configurations are specified, a sample size of -999 is returned.
#'
#' @export
#' @return
#' The returned value is a value indicating the power or required sample size.
#'
#' @references
#' \insertRef{cai04}{gap}
#'
#' @seealso [`pbsize`]
#'
#' @examples
#' \dontrun{
#' # Table 1 of Cai & Zeng (2004).
#' alpha <- 0.05
#' table1 <- rbind(
#'   transform(
#'     within(expand.grid(
#'       pD = c(0.10, 0.05),
#'       p1 = c(0.3, 0.5),
#'       theta = c(0.5, 1.0),
#'       q = c(0.1, 0.2)
#'     ), {
#'       n <- 1000
#'       power <- mapply(ccsize,
#'         n = n, q = q, pD = pD, p1 = p1, theta = theta,
#'         MoreArgs = list(alpha = alpha)
#'       )
#'     }),
#'     power = signif(power, 3)
#'   ),
#'   transform(
#'     within(expand.grid(
#'       pD = c(0.05, 0.01),
#'       p1 = c(0.3, 0.5),
#'       theta = c(0.5, 1.0),
#'       q = c(0.01, 0.02)
#'     ), {
#'       n <- 5000
#'       power <- mapply(ccsize,
#'         n = n, q = q, pD = pD, p1 = p1, theta = theta,
#'         MoreArgs = list(alpha = alpha)
#'       )
#'     }),
#'     power = signif(power, 3)
#'   )
#' )
#'
#' # ARIC study
#' aric <- within(
#'   data.frame(
#'     n = 15792,
#'     pD = 0.03,
#'     p1 = 0.25,
#'     hr = c(1.35, 1.40, 1.45),
#'     q = c(1463, 722, 468) / 15792
#'   ), {
#'     alpha <- 0.05
#'     beta <- 0.2
#'     power <- mapply(ccsize,
#'       n = n, q = q, pD = pD, p1 = p1, theta = log(hr),
#'       MoreArgs = list(alpha = alpha)
#'     )
#'     ssize <- mapply(ccsize,
#'       n = n, q = q, pD = pD, p1 = p1, theta = log(hr),
#'       MoreArgs = list(alpha = alpha, beta = beta, power = FALSE)
#'     )
#'     power <- signif(power, 3)
#'   }
#' )
#'
#' # EPIC study
#' epic <- within(
#'   expand.grid(
#'     pD = c(0.3, 0.2, 0.1, 0.05),
#'     p1 = seq(0.1, 0.5, by = 0.1),
#'     hr = seq(1.1, 1.4, by = 0.1)
#'   ), {
#'     n <- 25000
#'     q <- 0.1
#'     alpha <- 5e-8
#'     beta <- 0.2
#'     ssize <- mapply(ccsize,
#'       n = n, q = q, pD = pD, p1 = p1, theta = log(hr),
#'       MoreArgs = list(alpha = alpha, beta = beta, power = FALSE)
#'     )
#'   }
#' )
#' epic <- subset(epic, !is.na(ssize) & ssize > 0)
#'
#' # exhaustive search
#' search <- within(
#'   expand.grid(
#'     pD = c(0.3, 0.2, 0.1, 0.05),
#'     p1 = seq(0.1, 0.5, by = 0.1),
#'     hr = seq(1.1, 1.4, by = 0.1),
#'     q  = seq(0.01, 0.5, by = 0.01)
#'   ), {
#'     n <- 25000
#'     alpha <- 5e-8
#'     power <- mapply(ccsize,
#'       n = n, q = q, pD = pD, p1 = p1, theta = log(hr),
#'       MoreArgs = list(alpha = alpha)
#'     )
#'     nq <- n * q
#'   }
#' )
#' }
#' @author Jing Hua Zhao
#' @note Programmed for EPIC study.
#' keywords misc

ccsize <- function(n, q, pD, p1, theta, alpha,
                   beta = 0.2, power = TRUE, verbose = FALSE)
{
  p2 <- 1 - p1
  if (power)
  {
    # Equation (5): power
    z_alpha <- qnorm(alpha)
    m <- n * q
    z <- z_alpha +
         sqrt(m) * theta *
         sqrt(p1 * p2 * pD / (q + (1 - q) * pD))
    return(pnorm(z))
  } else {
    # Equation (6): sample size
    z_alpha <- qnorm(alpha, lower.tail = FALSE)
    z_beta  <- qnorm(beta,  lower.tail = FALSE)
    B <- (z_alpha + z_beta)^2 / (theta^2 * p1 * p2 * pD)
    nb <- ceiling(n * B * pD / (n - B * (1 - pD)))
    nb[nb > n] <- NA
    if (any(nb <= 0) && verbose)
      cat("Infeasible configuration\n")
    return(nb)
  }
}
