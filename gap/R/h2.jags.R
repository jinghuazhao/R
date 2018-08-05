h2.jags <- function(y,x,G,eps=0.0001,sigma.p=0,sigma.r=1,parms=c("b","p","r","h2"),...)
{
  M <- ncol(x)
  N <- length(y)
  data=list(M=M,N=N,eps=eps,y=y,x=x,CG=t(chol(G+diag(eps,N))))
  inits=function() list(b=rep(0,ncol(x)),sigma.p=sigma.p,sigma.r=sigma.r)
  modelfile <- system.file(package="gap","JAGS","h2.jags")
  jagsfit <- R2jags::jags(data, inits, parms, model.file=modelfile, ...)
}
