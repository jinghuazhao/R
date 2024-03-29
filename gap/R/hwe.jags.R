#' Hardy-Weinberg equlibrium test for a multiallelic marker using JAGS
#'
#' @param k number of alleles.
#' @param n a vector of k(k+1)/2 genotype counts.
#' @param delta initial parameter values.
#' @param lambda initial parameter values.
#' @param lambdamu initial parameter values.
#' @param lambdasd initial parameter values.
#' @param parms monitored parmeters.
#' @param ... parameters passed to jags, e.g., n.chains, n.burnin, n.iter.
#'
#' @details
#' Hardy-Weinberg equilibrium test.
#'
#' This function performs Bayesian Hardy-Weinberg equilibrium test, which mirrors [hwe.hardy],
#' another implementation for highly polymorphic markers when asymptotic results do not hold.
#'
#' @export
#' @return The returned value is a fitted model from jags().
#'
#' @references
#' \insertRef{wakefield10}{gap}
#'
#' @seealso [`hwe.hardy`]
#'
#' @examples
#' \dontrun{
#' ex1 <- hwe.jags(4,c(5,6,1,7,11,2,8,19,26,15))
#' print(ex1)
#' ex2 <- hwe.jags(2,c(49,45,6))
#' print(ex2)
#' ex3 <- hwe.jags(4,c(0,3,1,5,18,1,3,7,5,2),lambda=0.5,lambdamu=-2.95,lambdasd=1.07)
#' print(ex3)
#' ex4 <- hwe.jags(9,c(1236,120,3,18,0,0,982,55,7,249,32,1,0,12,0,2582,132,20,1162,29,
#'                     1312,6,0,0,4,0,4,0,2,0,0,0,0,0,0,0,115,5,2,53,1,149,0,0,4),
#'                 delta=c(1,1,1,1,1,1,1,1,1),lambdamu=-4.65,lambdasd=0.21)
#' print(ex4)
#' ex5 <- hwe.jags(8,n=c(
#'          3,
#'          4, 2,
#'          2, 2, 2,
#'          3, 3, 2, 1,
#'          0, 1, 0, 0, 0,
#'          0, 0, 0, 0, 0, 1,
#'          0, 0, 1, 0, 0, 0, 0,
#'          0, 0, 0, 2, 1, 0, 0, 0))
#' print(ex5)
#'
#' # Data and code accordining to the following URL,
#' # http://darwin.eeb.uconn.edu/eeb348-notes/testing-hardy-weinberg.pdf
#' hwe.jags.ABO <- function(n,...)
#' {
#'   hwe <- function() {
#'      # likelihood
#'      pi[1] <- p.a*p.a + 2*p.a*p.o
#'      pi[2] <- 2*p.a*p.b
#'      pi[3] <- p.b*p.b + 2*p.b*p.o
#'      pi[4] <- p.o*p.o
#'      n[1:4] ~ dmulti(pi[],N)
#'      # priors
#'      a1 ~ dexp(1)
#'      b1 ~ dexp(1)
#'      o1 ~ dexp(1)
#'      p.a <- a1/(a1 + b1 + o1)
#'      p.b <- b1/(a1 + b1 + o1)
#'      p.o <- o1/(a1 + b1 + o1)
#'   }
#'   hwd <- function() {
#'      # likelihood
#'      pi[1] <- p.a*p.a + f*p.a*(1-p.a) + 2*p.a*p.o*(1-f)
#'      pi[2] <- 2*p.a*p.b*(1-f)
#'      pi[3] <- p.b*p.b + f*p.b*(1-p.b) + 2*p.b*p.o*(1-f)
#'      pi[4] <- p.o*p.o + f*p.o*(1-p.o)
#'      n[1:4] ~ dmulti(pi[],N)
#'      # priors
#'      a1 ~ dexp(1)
#'      b1 ~ dexp(1)
#'      o1 ~ dexp(1)
#'      p.a <- a1/(a1 + b1 + o1)
#'      p.b <- b1/(a1 + b1 + o1)
#'      p.o <- o1/(a1 + b1 + o1)
#'      f ~ dunif(0,1)
#'   }
#'   N <- sum(n)
#'   ABO.hwe <- R2jags::jags(list(n=n,N=N),,c("pi","p.a","p.b","p.o"),hwe,...)
#'   ABO.hwd <- R2jags::jags(list(n=n,N=N),,c("pi","p.a","p.b","p.o","f"),hwd,...)
#'   invisible(list(hwe=ABO.hwe,hwd=ABO.hwd))
#' }

#' hwe.jags.ABO.results <- hwe.jags.ABO(n=c(862, 131, 365, 702))
#' hwe.jags.ABO.results
#' }
#'
#' @author Jing Hua Zhao, Jon Wakefield
#' @keywords htest

hwe.jags <- function(k,n,delta=rep(1/k,k),lambda=0,lambdamu=-1,lambdasd=1,
                     parms=c("p","f","q","theta","lambda"), ...)
{
  ncell <- k*(k+1)/2
  N <- sum(n)
  data <- list(k=k,n=n,ncell=ncell,N=N,lambdamu=lambdamu,lambdasd=lambdasd)
  inits <- function() list(delta=delta,lambda=lambda)
  modelfile <- system.file(package="gap","JAGS","hwe.jags")
  jagsfit <- R2jags::jags(data, inits, parms, modelfile, ...)
}
