#' Heritability estimation based on genomic relationship matrix using JAGS
#'
#' Heritability estimation based on genomic relationship matrix using JAGS.
#'
#' @param y outcome vector.
#' @param x covariate matrix.
#' @param G genomic relationship matrix.
#' @param eps a positive diagonal perturbation to G.
#' @param sigma.p initial parameter values.
#' @param sigma.r initial parameter values.
#' @param parms monitored parmeters.
#' @param ... parameters passed to jags, e.g., n.chains, n.burnin, n.iter.
#'
#' @details
#' This function performs Bayesian heritability estimation using genomic relationship matrix.
#'
#' @export
#' @return The returned value is a fitted model from jags().
#'
#' @references
#' Zhao JH, Luan JA, Congdon P (2018). Bayesian linear mixed models with polygenic effects. J Stat Soft 85(6):1-27, \doi{10.18637/jss.v085.i06}.
#'
#' @examples
#' \dontrun{
#' require(gap.datasets)
#' set.seed(1234567)
#' meyer <- within(meyer,{
#'     y[is.na(y)] <- rnorm(length(y[is.na(y)]),mean(y,na.rm=TRUE),sd(y,na.rm=TRUE))
#'     g1 <- ifelse(generation==1,1,0)
#'     g2 <- ifelse(generation==2,1,0)
#'     id <- animal
#'     animal <- ifelse(!is.na(animal),animal,0)
#'     dam <- ifelse(!is.na(dam),dam,0)
#'     sire <- ifelse(!is.na(sire),sire,0)
#' })
#' G <- kin.morgan(meyer)$kin.matrix*2
#' library(regress)
#' r <- regress(y~-1+g1+g2,~G,data=meyer)
#' r
#' with(r,h2G(sigma,sigma.cov))
#' eps <- 0.001
#' y <- with(meyer,y)
#' x <- with(meyer,cbind(g1,g2))
#' ex <- h2.jags(y,x,G,sigma.p=0.03,sigma.r=0.014)
#' print(ex)
#' }
#'
#' @author Jing Hua Zhao
#' keywords htest

h2.jags <- function(y,x,G,eps=0.0001,sigma.p=0,sigma.r=1,parms=c("b","p","r","h2"),...)
{
  M <- ncol(x)
  N <- length(y)
  data=list(M=M,N=N,eps=eps,y=y,x=x,CG=t(chol(G+diag(eps,N))))
  inits=function() list(b=rep(0,ncol(x)),sigma.p=sigma.p,sigma.r=sigma.r)
  modelfile <- system.file(package="gap","JAGS","h2.jags")
  jagsfit <- R2jags::jags(data, inits, parms, model.file=modelfile, ...)
}
