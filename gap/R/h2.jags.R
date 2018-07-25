h2.jags <- function(y,x,G,eps=0.0001,sigma.p=0,sigma.r=1,parms=c("b","p","r","h2"),...)
{
  N <- length(y)
  data=list(N=N,eps=eps,y=y,x=x,CG=t(chol(G+diag(eps,N))))
  inits=function() list(b=rep(0,ncol(x)),sigma.p=sigma.p,sigma.r=sigma.r)
  modelfile <- function()
  {
    for (i in 1:2)
    {
        b[i] ~ dnorm(0, 0.001)
    }
    sigma.p ~ dunif(0,1000)
    sigma.r ~ dunif(0,1000)
    p <- pow(sigma.p, 2)
    r <- pow(sigma.r, 2)
    tau <- pow(sigma.r, -2)
    g[1:N] <- sigma.p*CG[,]%*%z[]
    for (i in 1:N)
    {
        z[i] ~ dnorm(0,1)
    }
    for (i in 1:N)
    {
        xb[i] <- inprod(b[], x[i,])
        y[i] ~ dnorm(xb[i] + g[i], tau)
    }
    h2 <- p / (p * (1 + eps) + r)
  }
  jagsfit <- R2jags::jags(data, inits, parms, model.file=modelfile, ...)
}
