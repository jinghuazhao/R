#
# Functions written by Jon Wakefield 17th April, 2007.
#
# For more details see "A Bayesian Measure of the Probability of False
# Discovery in Genetic Epidemiology Studies", under revision for American
# Journal of Human Genetics.
#
# This function calculates BFDP, the approximate Pr( H0 | thetahat ),
# given an estiamte of the log relative risk, thetahat, the variance of
# this estimate, V, the prior variance, W, and the prior probability of
# a non-null association.
#
BFDPfunV <- function(thetahat,V,W,pi1){
  pH0 <- dnorm(thetahat,m=0,s=sqrt(V))
  postvar <- V + W
  pH1 <- dnorm(thetahat,m=0,s=sqrt(postvar))
  BF <- pH0/pH1; PO <- (1-pi1)/pi1
  BFDP <- BF*PO/(BF*PO+1)
  list(BF=BF,pH0=pH0,pH1=pH1,BFDP=BFDP)
}
#
# This function calculates BFDP, the approximate Pr( H0 | thetahat ),
# given an estiamte of the relative risk, RRhat, and the upper point of
# a 95% confidence interval RRhi, the prior variance, W, and the prior
# probability of a non-null association.
#
BFDPfunCI <- function(RRhat,RRhi,W,pi1){
  thetahat <- log(RRhat)
  V <- ((log(RRhi)-log(RRhat))/1.96)^2
  pH0 <- dnorm(thetahat,m=0,s=sqrt(V))
  postvar <- V + W
  pH1 <- dnorm(thetahat,m=0,s=sqrt(postvar))
  BF <- pH0/pH1; pi0 <- 1-pi1; PO <- pi0/(1-pi0)
  BFDP <- BF*PO/(BF*PO+1)
  list(BF=BF,pH0=pH0,pH1=pH1,BFDP=BFDP)
}
BFDPfunV2 <- function(thetahat1,thetahat2,V1,V2,W,pi1){
  BF1 <-dnorm(thetahat1,m=0,s=sqrt(V1))/dnorm(thetahat1,m=0,s=sqrt(V1+W))
  BF2 <-dnorm(thetahat2,m=0,s=sqrt(V2))/dnorm(thetahat2,m=0,s=sqrt(V2+W))
  Rfac <- W/(V1*W + V2*W + V1*V2)
  term1 <- thetahat1^2/V1^2
  term2 <- 2*thetahat1*thetahat2/(V1*V2)
  term3 <- thetahat2^2/V2^2
  BF12 <- sqrt(W/(Rfac*V1*V2))*exp(-.5*(term1+term2+term3)*V1*V2*Rfac)
  BF2given1 <- BF12/BF1
  PO <- (1-pi1)/pi1
  BFDP1 <- BF1*PO/(BF1*PO+1)
  BFDP2 <- BF2*PO/(BF2*PO+1)
  BFDP12 <- BF12*PO/(BF12*PO+1)
  smallr1 <- W/(V1+W)
  smallr2 <- W/(V2+W)
# pm1 and pv1 are the posterior median and posterior variance from the first study only
  pm1 <- smallr1*thetahat1; pv1 <- smallr1*V1
# pm2 and pv2 are the posterior median and posterior variance from the second study only
  pm2 <- smallr2*thetahat2; pv2 <- smallr2*V2
# pm12 and pv12 are the posterior median and posterior variance from the first two studies combined
  pm12 <- thetahat1*V2*Rfac + thetahat2*V1*Rfac; pv12 <- Rfac*V1*V2
  list(BF1=BF1,BF2=BF2,BF12=BF12,BF2given1=BF2given1,BFDP1=BFDP1,BFDP2=BFDP2,
       pm1=pm1,pv1=pv1,pm2=pm2,pv2=pv2,pm12=pm12,pv12=pv12)
}
BFDPfunCI2 <- function(RRhat1,RRhat2,RRhi1,RRhi2,W,pi1){
  thetahat1 <- log(RRhat1); thetahat2 <- log(RRhat2) 
  V1 <- ((log(RRhi1)-log(RRhat1))/1.96)^2
  V2 <- ((log(RRhi2)-log(RRhat2))/1.96)^2
  BF1 <-dnorm(thetahat1,m=0,s=sqrt(V1))/dnorm(thetahat1,m=0,s=sqrt(V1+W))
  BF2 <-dnorm(thetahat2,m=0,s=sqrt(V2))/dnorm(thetahat2,m=0,s=sqrt(V2+W))
  Rfac <- W/(V1*W + V2*W + V1*V2)
  term1 <- thetahat1^2/V1^2
  term2 <- 2*thetahat1*thetahat1/(V1*V2)
  term3 <- thetahat2^2/V2^2
  BF12 <- sqrt(W/(Rfac*V1*V2))*exp(-.5*(term1+term2+term3)*V1*V2*Rfac)
  BF2given1 <- BF12/BF1
  PO <- (1-pi1)/pi1
  BFDP1 <- BF1*PO/(BF1*PO+1)
  BFDP2 <- BF2*PO/(BF2*PO+1)
  BFDP12 <- BF12*PO/(BF12*PO+1)
  smallr1 <- W/(V1+W)
  smallr2 <- W/(V2+W)
# pm1 and pv1 are the posterior median and posterior variance from the first study only
  pm1 <- smallr1*thetahat1; pv1 <- smallr1*V1
# pm2 and pv2 are the posterior median and posterior variance from the second study only
  pm2 <- smallr2*thetahat2; pv2 <- smallr2*V2
# pm12 and pv12 are the posterior median and posterior variance from the first two studies combined
  pm12 <- thetahat1*V2*Rfac + thetahat2*V1*Rfac; pv12 <- Rfac*V1*V2
  list(BF1=BF1,BF2=BF2,BF12=BF12,BF2given1=BF2given1,BFDP1=BFDP1,BFDP2=BFDP2,BFDP12=BFDP12,
       pm1=pm1,pv1=pv1,pm2=pm2,pv2=pv2,pm12=pm12,pv12=pv12,V1=V1,V2=V2)
}

