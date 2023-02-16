#' Bayesian false-discovery probability
#'
#' @param a parameter value at which the power is to be evaluated.
#' @param b the variance for a, or the uppoer point (\eqn{RR_{hi}}{RRhi}) of a 95%CI if logscale=FALSE.
#' @param pi1 the prior probabiility of a non-null association.
#' @param W the prior variance.
#' @param logscale FALSE=the orginal scale, TRUE=the log scale.
#'
#' @details
#' This function calculates BFDP, the approximate \eqn{P(H_0|\hat\theta)}{Pr( H0 | thetahat )},
#' given an estiamte of the log relative risk, \eqn{\hat\theta}{thetahat}, the variance of
#' this estimate, \eqn{V}, the prior variance, \eqn{W}, and the prior probability of
#' a non-null association. When logscale=TRUE, the function accepts an estimate of the relative
#' risk, \eqn{\hat{RR}}{RRhat}, and the upper point of a 95% confidence interval \eqn{RR_{hi}}{RRhi}.
#'
#' @export
#' @return
#' The returned value is a list with the following components:
#' PH0. probability given a,b).
#' PH1. probability given a,b,W).
#' BF. Bayes factor, \eqn{P_{H_0}/P_{H_1}}{PH0/PH1}.
#' BFDP. Bayesian false-discovery probability.
#' ABF. approxmiate Bayes factor.
#' ABFDP. approximate Bayesian false-discovery probability.
#'
#' @references
#' \insertRef{wakefield07}{gap}
#'
#' @seealso [`FPRP`]
#'
#' @examples
#' \dontrun{
#' # Example from BDFP.xls by Jon Wakefield and Stephanie Monnier
#' # Step 1 - Pre-set an BFDP-level threshold for noteworthiness: BFDP values below this
#' #          threshold are noteworthy
#' # The threshold is given by R/(1+R) where R is the ratio of the cost of a false
#' # non-discovery to the cost of a false discovery
#'
#' T <- 0.8
#'
#' # Step 2 - Enter up values for the prior that there is an association
#'
#' pi0 <- c(0.7,0.5,0.01,0.001,0.00001,0.6)
#'
#' # Step 3 - Enter the value of the OR that is the 97.5% point of the prior, for example
#' #          if we pick the value 1.5 we believe that the prior probability that the
#' #          odds ratio is bigger than 1.5 is 0.025.
#'
#' ORhi <- 3
#'
#' W <- (log(ORhi)/1.96)^2
#' W
#'
#' # Step 4 - Enter OR estimate and 95% confidence interval (CI) to obtain BFDP
#'
#' OR <- 1.316
#' OR_L <- 1.10
#' OR_U <- 2.50
#' logOR <- log(OR)
#' selogOR <- (log(OR_U)-log(OR))/1.96
#' r <- W/(W+selogOR^2)
#' r
#' z <- logOR/selogOR
#' z
#' ABF <- exp(-z^2*r/2)/sqrt(1-r)
#' ABF
#' FF <- (1-pi0)/pi0
#' FF
#' BFDPex <- FF*ABF/(FF*ABF+1)
#' BFDPex
#' pi0[BFDPex>T]
#'
#' ## now turn to BFDP
#'
#' pi0 <- c(0.7,0.5,0.01,0.001,0.00001,0.6)
#' ORhi <- 3
#' OR <- 1.316
#' OR_U <- 2.50
#' W <- (log(ORhi)/1.96)^2
#' z <- BFDP(OR,OR_U,pi0,W)
#' z
#' }
#' @author Jon Wakefield, Jing Hua Zhao
#' @note Adapted from BFDP functions by Jon Wakefield on 17th April, 2007.
#' @keywords models

BFDP <- function(a,b,pi1,W,logscale=FALSE)
{
  if (logscale) {
     thetahat <- a
     V <- b
     PO <- (1-pi1)/pi1
  } else {
     thetahat <- log(a)
     V <- ((log(b)-log(a))/1.96)^2
     pi0 <- 1-pi1
     PO <- pi0/(1-pi0)
  }
  postvar <- V + W
  r <- W / postvar
  z <- thetahat / sqrt(V)
  ABF <- exp(-r*z^2/2) / sqrt(1-r)
  ABFDP <- ABF*PO/(ABF*PO+1)
  pH0 <- dnorm(thetahat,mean=0,sd=sqrt(V))
  pH1 <- dnorm(thetahat,mean=0,sd=sqrt(postvar))
  BF <- pH0/pH1
  BFDP <- BF*PO/(BF*PO+1)
  invisible(list(pH0=pH0,pH1=pH1,BF=BF,BFDP=BFDP,ABF=ABF,ABFDP=ABFDP))
}

# 31-7-2007
