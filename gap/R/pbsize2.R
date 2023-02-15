#' A simple chi-squared test of two proportions
#' @noRd

x2 <- function(p1,p2,n1,n2)
# 3-3-2008 MRC-Epid JHZ
{
   m <- p1-p2
   v <- p1*(1-p1)/n1+p2*(1-p2)/n2
   invisible(m*m/v)
}

#' Power for case-control association design
#'
#' @param N The sample size.
#' @param fc The proportion of cases in the sample.
#' @param alpha Type I error rate.
#' @param gamma The genotype relative risk (GRR).
#' @param p Frequency of the risk allele.
#' @param kp The prevalence of disease in the population.
#' @param model Disease model, i.e., "multiplicative","additive","dominant","recessive","overdominant".
#'
#' @details
#' This extends [`pbsize`] from a multiplicative model for a case-control
#' design under a range of disease models. Essentially, for given sample sizes(s), a
#' proportion of which (fc) being cases, the function calculates power estimate for a
#' given type I error (alpha), genotype relative risk (gamma), frequency of the risk
#' allele (p), the prevalence of disease in the population (kp) and optionally a disease
#' model (model). A major difference would be the consideration of case/control
#' ascertainment in [`pbsize`].
#'
#' Internally, the function obtains a baseline risk to make the disease model consistent
#' with Kp as in [`tscc`] and should produce accurate power estimate. It
#' provides power estimates for given sample size(s) only.
#'
#' @export
#' @return
#' The returned value is the power for the specified design.
#'
#' @seealso The design follows that of [`pbsize`].
#'
#' @examples
#' \dontrun{
#' # single calculation
#' m <- c("multiplicative","recessive","dominant","additive","overdominant")
#' for(i in 1:5) print(pbsize2(N=50,alpha=5e-2,gamma=1.1,p=0.1,kp=0.1, model=m[i]))
#'
#' # a range of sample sizes
#' pbsize2(p=0.1, N=c(25,50,100,200,500), gamma=1.2, kp=.1, alpha=5e-2, model='r')
#'   
#' # a power table
#' m <- sapply(seq(0.1,0.9, by=0.1),
#'             function(x) pbsize2(p=x, N=seq(100,1000,by=100),
#'                         gamma=1.2, kp=.1, alpha=5e-2, model='recessive'))
#' colnames(m) <- seq(0.1,0.9, by=0.1)
#' rownames(m) <- seq(100,1000,by=100)
#' print(round(m,2))
#' }
#'
#' @keywords misc

pbsize2 <- function(N,fc=0.5,alpha=0.05,gamma=4.5,p=0.15,kp=0.1,model="additive")
{
   z <- KCC(model,gamma,p,kp)
   pp <- function(ssize)
   {
      x <- x2(z$pprime,z$p,fc*ssize,(1-fc)*ssize)
      q <- qchisq(1-alpha,1)
      power <- pchisq(q,1,x,lower.tail=FALSE)
   }
   sapply(N,pp)
}
