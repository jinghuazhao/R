#' genomic control
#'
#' @param data the data matrix.
#' @param zeta program constant with default value 1000.
#' @param kappa multiplier in prior for mean with default value 4.
#' @param tau2 multiplier in prior for variance with default value 1.
#' @param epsilon prior probability of marker association with default value 0.01.
#' @param ngib number of Gibbs steps, with default value 500.
#' @param burn number of burn-ins with default value 50.
#' @param idum seed for pseudorandom number sequence.
#'
#' @details
#' The Bayesian genomic control statistics with the following parameters,
#'
#' \tabular{ll}{
#' n \tab number of loci under consideration \cr
#' lambdahat \tab median(of the n trend statistics)/0.46 \cr
#'   \tab Prior for noncentrality parameter Ai is \cr
#'   \tab Normal(sqrt(lambdahat)kappa,lambdahat*tau2) \cr
#' kappa \tab multiplier in prior above, set at 1.6 * sqrt(log(n)) \cr
#' tau2  \tab multiplier in prior above \cr
#' epsilon \tab prior probability a marker is associated, set at 10/n \cr
#' ngib \tab number of cycles for the Gibbs sampler after burn in \cr
#' burn \tab number of cycles for the Gibbs sampler to burn in
#' }
#'
#' Armitage's trend test along with the posterior probability that each marker is
#' associated with the disorder is given. The latter is not a p-value but any value
#' greater than 0.5 (pout) suggests association. 
#'
#' @source \url{https://www.cmu.edu/dietrich/statistics-datascience/index.html}
#'
#' @return The returned value is a list containing:
#' - deltot the probability of being an outlier.
#' - x2 the \eqn{\chi^2}{chi-squared} statistic.
#' - A the A vector.
#'
#' @references
#' \insertRef{devlin99}{gap}
#'
#' @examples
#' \dontrun{
#' test<-c(1,2,3,4,5,6,  1,2,1,23,1,2, 100,1,2,12,1,1, 
#'         1,2,3,4,5,61, 1,2,11,23,1,2, 10,11,2,12,1,11)
#' test<-matrix(test,nrow=6,byrow=T)
#' gcontrol(test)
#' }
#' @author Bobby Jones, Jing Hua Zhao
#'
#' @note Adapted from gcontrol by Bobby Jones and Kathryn Roeder, 
#' use -Dexecutable for standalone program, function getnum in the original 
#' code needs \%*s to skip id string
#' @export
#' @keywords models

gcontrol<-function(data,zeta=1000,kappa=4,tau2=1,epsilon=0.01,ngib=500,burn=50,idum=2348)
{
  nkdata <- dim(data)[1]
  deltot <- rep(0,nkdata)
  kdata <- as.matrix(data)
  x <- a <- rep(0,nkdata)
  z<-.C("gcontrol_c",kdata=as.double(kdata),nkdata=as.integer(nkdata),
        zeta=as.double(zeta),kappa=as.double(kappa),tau2=as.double(tau2),
        epsilon=as.double(epsilon),ngib=as.integer(ngib),burn=as.integer(burn),
        idumR=as.integer(idum),deltot=as.double(deltot),x=as.double(array(x)),A=as.double(a),PACKAGE="gap")

  nkdata6 <- nkdata/6
  list(deltot=z$deltot[1:nkdata6],x2=z$x[1:nkdata6],A=z$A[1:nkdata6])
}
