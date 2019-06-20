# 
# Data:
#	   CC   CT   TT
# Cases     6    8   75
# Controls 10   66  163
#
#  Likelihood first
# 
x <- c(0,1,2)
y <- c(6,8,75)
z <- c(10,66,163)
logitmod <- glm(cbind(y,z)~x,family="binomial")
thetahat <- logitmod$coeff[2]		# Log odds ratio
V <- vcov(logitmod)[2,2]		# se^2
exp(thetahat) 
exp(thetahat-1.96*sqrt(V))
exp(thetahat+1.96*sqrt(V))
summary(logitmod)
# 
# Bayesian
#
source("http://faculty.washington.edu/jonno/BFDP.R")
pi1 <- c(1/100,1/1000,1/10000,1/100000) # Prior on alternative
Upper975 <- 1.5
W <- (log(Upper975)/1.96)^2 # 97.5 point of prior is log(1.5) so that we 
                            # believe with prior prob 0.95 that the odds
                            # corresponding to one T allele more,
			    # lies in (2/3,1.5)
BFcall <- BFDPfunV(thetahat,V,W,pi1)
r <- W/(V+W)
exp(r*thetahat)
exp(r*thetahat-1.96*sqrt(r*V))
exp(r*thetahat+1.96*sqrt(r*V))
cat("log odds, standard error and ratio: ",thetahat,sqrt(V),
     thetahat/sqrt(V),"\n")
# Posterior distribution is asymptotically lognormal with mean r x thetahat
# and variance r x thetahat
RRseq <- seq(.1,4,.01)
postseq <- dlnorm(RRseq,mean=r*thetahat,sd=sqrt(r*V))
plot(postseq~RRseq,type="n",xlab="Odds Ratio",ylab="Posterior")
priorseq <- dlnorm(RRseq,mean=0,sd=sqrt(W))
lines(priorseq~RRseq,lty=2)
lines(postseq~RRseq)
legend("topright",legend=c("Prior Distribution","Posterior Distribution"),
lty=2:1,bty="n")
cat("Bayes Factor of H0 over H1: ",BFcall$BF,"\n")
cat("Posterior probs of the null, given priors on the null of: ",1-pi1," are ",BFcall$BFDP,"\n")
