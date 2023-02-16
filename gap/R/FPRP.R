#' False-positive report probability
#'
#' @param a parameter value at which the power is to be evaluated.
#' @param b the variance for a, or the uppoer point of a 95%CI if logscale=FALSE.
#' @param pi0 the prior probabiility that \eqn{H_0}{H0} is true.
#' @param ORlist a vector of ORs that is most likely.
#' @param logscale FALSE=a,b in orginal scale, TRUE=a, b in log scale.
#'
#' @details
#' The function calculates the false positive report probability (FPRP), the probability of no true
#' association beteween a genetic variant and disease given a statistically significant finding,
#' which depends not only on the observed P value but also on both the prior probability that the
#' assocition is real and the statistical power of the test. An associate result is the
#' false negative reported probability (FNRP).  See example for the recommended steps.
#'
#' The FPRP and FNRP are derived as follows. Let \eqn{H_0}=null hypothesis (no association),
#' \eqn{H_A}=alternative hypothesis (association). Since classic frequentist theory considers 
#' they are fixed, one has to resort to Bayesian framework by introduing prior,
#' \eqn{\pi=P(H_0=TRUE)=P(association)}. Let \eqn{T}=test statistic, and \eqn{P(T>z_\alpha|H_0=TRUE)=P(rejecting\ 
#' H_0|H_0=TRUE)=\alpha}, \eqn{P(T>z_\alpha|H_0=FALSE)=P(rejecting\ H_0|H_A=TRUE)=1-\beta}. The joint
#' probability of test and truth of hypothesis can be expressed by \eqn{\alpha}, \eqn{\beta} and \eqn{\pi}.
#'
#' Joint probability of significance of test and truth of hypothesis
#' \tabular{llll}{
#' Truth of \eqn{H_A} \tab significant \tab nonsignificant \tab Total\cr 
#' TRUE \tab \eqn{(1-\beta)\pi} \tab \eqn{\beta\pi} \tab \eqn{\pi}\cr
#' FALSE \tab \eqn{\alpha (1-\pi)} \tab \eqn{(1-\alpha)(1-\pi)} \tab \eqn{1-\pi}\cr
#' Total \tab \eqn{(1-\beta)\pi+\alpha (1-\pi)} \tab \eqn{\beta\pi+(1-\alpha)(1-\pi)} \tab 1\cr
#' }
#'
#' We have \eqn{FPRP=P(H_0=TRUE|T>z_\alpha)= 
#' \alpha(1-\pi)/[\alpha(1-\pi)+(1-\beta)\pi]=\{1+\pi/(1-\pi)][(1-\beta)/\alpha]\}^{-1}}
#' and similarly \eqn{FNRP=\{1+[(1-\alpha)/\beta][(1-\pi)/\pi]\}^{-1}}.
#'
#' @export
#' @return
#' The returned value is a list with compoents,
#' p p value corresponding to a,b.
#' power the power corresponding to the vector of ORs.
#' FPRP False-positive report probability.
#' FNRP False-negative report probability.
#' 
#' @references
#' \insertRef{wacholder04}{gap}
#'
#' @seealso [`BFDP`]
#'
#' @examples
#' \dontrun{
#' # Example by Laure El ghormli & Sholom Wacholder on 25-Feb-2004
#' # Step 1 - Pre-set an FPRP-level criterion for noteworthiness
#'
#' T <- 0.2
#'
#' # Step 2 - Enter values for the prior that there is an association
#'
#' pi0 <- c(0.25,0.1,0.01,0.001,0.0001,0.00001)
#'
#' # Step 3 - Enter values of odds ratios (OR) that are most likely, assuming that
#' #          there is a non-null association
#'
#' ORlist <- c(1.2,1.5,2.0)
#'
#' # Step 4 - Enter OR estimate and 95% confidence interval (CI) to obtain FPRP 												
#'
#' OR <- 1.316
#' ORlo <- 1.08
#' ORhi <- 1.60
#'
#' logOR <- log(OR)
#' selogOR <- abs(logOR-log(ORhi))/1.96
#' p <- ifelse(logOR>0,2*(1-pnorm(logOR/selogOR)),2*pnorm(logOR/selogOR))
#' p
#' q <- qnorm(1-p/2)
#' POWER <- ifelse(log(ORlist)>0,1-pnorm(q-log(ORlist)/selogOR),
#'                 pnorm(-q-log(ORlist)/selogOR))
#' POWER
#' FPRPex <- t(p*(1-pi0)/(p*(1-pi0)+POWER\%o\%pi0))
#' row.names(FPRPex) <- pi0
#' colnames(FPRPex) <- ORlist
#' FPRPex
#' FPRPex>T
#'
#' ## now turn to FPRP
#' OR <- 1.316
#' ORhi <- 1.60
#' ORlist <- c(1.2,1.5,2.0)
#' pi0 <- c(0.25,0.1,0.01,0.001,0.0001,0.00001)
#' z <- FPRP(OR,ORhi,pi0,ORlist,logscale=FALSE)
#' z
#' }
#'
#' @author Jing Hua Zhao
#' @keywords models

FPRP <- function (a, b, pi0, ORlist, logscale = FALSE)
{
    if (logscale) {
        thetahat <- a
        V <- b
    }
    else {
        thetahat <- log(a)
        V <- ((log(b) - log(a))/1.96)^2
    }
    s <- sqrt(V)
    if (thetahat > 0) {
        p <- 2 * (1 - pnorm(thetahat/s))
        q <- qnorm(1 - p/2)
        b1 <- 1 - pnorm(q - log(ORlist)/s)
    }
    else {
        p <- 2 * pnorm(thetahat/s)
        q <- qnorm(1 - p/2)
        b1 <- pnorm(-q - log(ORlist)/s)
    }
    FPRP <- t(p * (1 - pi0)/(p * (1 - pi0) + b1 %o% pi0))
    row.names(FPRP) <- pi0
    colnames(FPRP) <- ORlist
    FNRP <- 1/(1 + (1 - pi0)/pi0 %o% ((1 - p)/(1 - b1)))
    row.names(FNRP) <- pi0
    colnames(FNRP) <- ORlist
    invisible(list(p = p, power = b1, FPRP = FPRP, FNRP = FNRP))
}

# 2-8-2007
