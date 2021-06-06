#' Heritability under the liability threshold model
#'
#' @md
#' @param K Disease prevalence.
#' @param P Phenotypeic variance.
#' @param h2 Heritability estimate.
#' @param se Standard error.
#' @param verbose Detailed output.
#' @export
#' @return A list of the input heritability estimate/standard error and their counterpart under liability threshold model, the normal deviate..

h2l <- function(K=0.05,P=0.5,h2,se,verbose=TRUE)
{
  x <- qnorm(1-K)
  z <- dnorm(x)
  1/sqrt(2*pi)*exp(-x^2/2)
  fK <- (K*(1-K)/z)^2
  fP <- P*(1-P)
  f <- fK/fP
  h2l <- f*h2
  sel <- f*se
  z2 <- K^2*(1-K)^2/(f*fP)
  x2 <- -log(2*pi*z2)
  if (verbose) {
     cat("K = ", K, "P = ", P, "\n")
     cat("h2 =",h2,"SE =",se,"h2l =",h2l,"SE =",sel,"\n")
  }
  invisible(list(h2=h2,se=se,h2l=h2l,sel=sel,cc=f,z=sqrt(x2)))
}

