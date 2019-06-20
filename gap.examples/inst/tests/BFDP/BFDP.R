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

