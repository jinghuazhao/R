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
