#' Test/Power calculation for mediating effect
#'
#' @param type string option: "test", "power".
#' @param n default sample size to be used for power calculation.
#' @param a regression coefficient from indepdendent variable to mediator.
#' @param sa SE(a).
#' @param b regression coefficient from mediator variable to outcome.
#' @param sb SE(b).
#' @param alpha size of siginficance test for power calculation.
#' @param fold fold change for power calculation, as appropriate for a range of sample sizes.
#'
#' @details
#' This function tests for or obtains power of mediating effect based on estimates of
#' two regression coefficients and their standard errors. Note that for binary outcome
#' or mediator, one should use log-odds ratio and its standard error.
#'
#' @export
#' @return The returned value are z-test and significance level for significant testing or sample size/power for a
#' given fold change of the default sample size.
#'
#' @references
#' \insertRef{freathy08}{gap}
#'
#' Kline RB. Principles and practice of structural equation modeling, Second Edition. The Guilford Press 2005.
#'
#' MacKinnon DP. Introduction to Statistical Mediation Analysis. Taylor & Francis Group 2008.
#'
#' Preacher KJ, Leonardelli GJ. Calculation for the Sobel Test-An interactive calculation tool for mediation tests
#' https://quantpsy.org/sobel/sobel.htm
#'
#' @seealso [`ccsize`]
#'
#' @examples
#' \dontrun{
#' ab()
#' n <- power <- vector()
#' for (j in 1:10)
#' {
#'    z <- ab(fold=j*0.01)
#'    n[j] <- z[1]
#'    power[j] <- z[2]
#' }
#' plot(n,power,xlab="Sample size",ylab="Power")
#' title("SNP-BMI-T2D association in EPIC-Norfolk study")
#' }
#'
#' @author Jing Hua Zhao
#' @keywords htest

ab <- function(type="power",n=25000,a=0.15,sa=0.01,b=log(1.19),sb=0.01,alpha=0.05,fold=1)
{
   ab <- a*b
   s <- sqrt(a^2*sb^2+b^2*sa^2)
   z <- ab/s
   if (type=="power") 
   {
      x2 <- z^2
      x <- qchisq(alpha,1,lower.tail=FALSE)
      power <- pchisq(x,1,ncp=x2*fold,lower.tail=FALSE)
      cat(fold*n, ",", power, "\n")
      stats <- c(fold*n,power)
   } else if (type=="test") stats <- c(z,2*pnorm(-abs(z)))
   else stop("Invalid option")
   invisible(stats)
}

# 10-11-2009 Modified from Stata code
