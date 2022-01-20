for(p in c("MASS")) {
   if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
      if (!requireNamespace(p, quietly = TRUE))
      warning(paste("hwe.cc needs package `", p, "' to be fully functional; please install", sep=""))
   }
}

#######################################################
# Cox procedure
#######################################################
 
## alpha is constrained by beta, gamma, q, and k
alpha.f <-function(b,r,q,k) { k/( (1-q)^2+2*q*(1-q)*b+q^2*r )}
 
 
pfix<-function(parm,low,high)

# adjust parm to lie between low and high
   {
      if(parm<low)return(low)
      if(parm>high)return(high)
      return(parm)  # otherwise parm is OK
   }

###################################################################
invlogit <- function (x=0){
      exp(x)/(1+exp(x))
    } 
logit <- function (p=0.5){
   log(p/(1-p)) 
   } 
se.invlogit <- function(logit.p, se.logit){
   f1 <- invlogit(logit.p)
   se.inv <- (f1-f1^2)*se.logit 
   se.inv 
   }

cov.invlogit <- function(logit.p1, logit.p2,  cov.logit){
   f1 <- invlogit(logit.p1)
   f2 <- invlogit(logit.p2)

   cov.inv <-  (f1-f1^2)*(f2-f2^2)*cov.logit  
   cov.inv 
   }

###################################################################
   
se.exp <- function(p, se.p){
   f1 <- exp(p)
   se.delta <- sqrt( (f1)^2*se.p^2 )
   se.delta 
   }

################################################################

## Cox T statistics 

# case and control data in order of AA,Aa,aa
Cox.T <- function(parms, case, control, k)
{
    b<-parms[1]
    r<-parms[2]
    q<-parms[3]

    b<-pfix(b,0.001,19.99)
    r<-pfix(r,0.001,19.99)
    q<-pfix(q,0.001,0.999)
        
    a<-k/( (1-q)^2+2*q*(1-q)*b+q^2*r )
    
    nt<-sum(case)
    nc<-sum(control)

    gAA<-(1-q)^2*a/k
    gAa<-2*q*(1-q)*a*b/k
    gaa<-q^2*a*r/k
    ecase<-c(gAA,gAa,gaa)*nt
    
    t<-sum((ecase-case)^2/ecase)
    
    gAA<-(1-q)^2*(1-a)/(1-k)
    gAa<-2*q*(1-q)*(1-a*b)/(1-k)
    gaa<-q^2*(1-a*r)/(1-k)
    econtrol<-c(gAA,gAa,gaa)*nc

    t<-t+sum((econtrol-control)^2/econtrol)
    
    return(t)
}

##  estimate functions 

Cox.est <- function(case,ctl,k0,initial) {

   fitout<-nlm(Cox.T,initial,case=case,control=ctl,k=k0,hessian=TRUE)

#### beta
          beta <- fitout$estimate[1]
          se.b <- sqrt(1/fitout$hessian[1,1])
### gamme
          gamma<-fitout$estimate[2]
          se.g <- sqrt(1/fitout$hessian[2,2])
### q
          q<-fitout$estimate[3]
          se.q<- sqrt(1/fitout$hessian[3,3])
### alpha
          alpha<- alpha.f(fitout$estimate[1], fitout$estimate[2],
                                 fitout$estimate[3],k0)
### Deviance
          deviance<-fitout$minimum[1]

    list(alpha=alpha, beta=c(beta,se.b), gamma=c(gamma,se.g),
         q=c(q,se.q),deviance=c(deviance, 1-pchisq(deviance,1)) ) 
}

###########################################################################


#######################################################
# recessive models 
#######################################################

######################################################################
## deviance, H0 and 1 parameters: (beta=1), gamma  
# recessive model
# case and control data in order of AA,Aa,aa
# parms=[p,gamma], logit and log transformed 


DevH0recessive <- function(parms, case, control, k)
{
    b<-1

    p<-parms[1]
    r<-parms[2]

    p<-pfix(p,-15,15)
    p <- exp(p)/(1+exp(p))

    r<-pfix(r,-10,10)
    r <- exp(r)

    q<-1-p
        
    a<-k/( (1-q)^2+2*q*(1-q)*b+q^2*r )
    
    nt<-sum(case)
    nc<-sum(control)

    gAA<-(1-q)^2*a/k
    gAa<-2*q*(1-q)*a*b/k
    gaa<-q^2*a*r/k
    ecase<-c(gAA,gAa,gaa)*nt

    if ( length(ecase[ecase<=0])>0 )
      {
      #   cat(lcase, "- -", ecase, "       1    \n")
       ll <- 1.0e10
       return(2*ll)
       }

    ll <-  sum(case * log(case/ecase))
    
    gAA<-(1-q)^2*(1-a)/(1-k)
    gAa<-2*q*(1-q)*(1-a*b)/(1-k)
    gaa<-q^2*(1-a*r)/(1-k)
    econtrol<-c(gAA,gAa,gaa)*nc

    if (length(econtrol[econtrol<=0])>0)
      {
      #   cat("b=", b, ", r=",  r,", q=", q, ", a=", a, "       2   \n")
      ll <- 1.0e10
      return(2*ll)
      }

    ll <- ll + sum(control * log(control/econtrol))

    return(2*ll)
}

DevH0recessive.est <- function(case,ctl,k0,initial) {

   fitout<-nlm(DevH0recessive,initial,case=case,control=ctl,k=k0,hessian=TRUE)

         hessinv<-MASS::ginv(fitout$hessian)     # generalized inverse of Hessian

#### beta fixed at 1
          beta <- 1 
          
### p
          p <- invlogit(fitout$estimate[1])
          se.p <- se.invlogit(fitout$estimate[1],sqrt(hessinv[1,1]))

          q <- 1-p
          se.q <- se.p

### gamma
          gamma <- exp(fitout$estimate[2])
          se.g <- se.exp(fitout$estimate[2],sqrt(hessinv[2,2]))

### alpha
          alpha<- alpha.f(beta, gamma, q, k0)

### Deviance
          deviance<-fitout$minimum[1]

    list(alpha=alpha, beta=c(beta,0), gamma=c(gamma,se.g),
         q=c(q,se.q),deviance=c(deviance, 1-pchisq(deviance,2)) ) 
}


##########################################################################

## deviance, H_a (DHW) and 1 parameter: gamma, and
## beta=1 fixed, recessive model  
# case and control data in order of AA,Aa,aa
# parms=[p0,p1,gamma], logit or log transformed

DevHaGrecessive <- function(parms, case, control, k)
{
    b<-1

    p0 <- parms[1]
    p1 <- parms[2]
    r <- parms[3]

    p0 <- pfix(p0,-15,15)
    p0 <- exp(p0)/(1+exp(p0))

    p1 <- pfix(p1,-15,15)
    p1 <- exp(p1)/(1+exp(p1))

    r<- pfix(r,-10,10)
    r <- exp(r)
 
    pAA <- p0
    pAa <- p1
    paa <- 1-p0-p1

    a <- k/( pAA+ pAa*b+ paa*r )
    
    nt<-sum(case)
    nc<-sum(control)

    gAA <- pAA*a/k
    gAa <- pAa*a*b/k
    gaa <- paa*a*r/k
    ecase<-c(gAA,gAa,gaa)*nt

    if ( length(ecase[ecase<=0])>0 )
      {
      #   cat(lcase, "- -", ecase, "       1    \n")
       ll <- 1.0e10
       return(2*ll)
       }

    ll <-  sum(case * log(case/ecase))
    
    gAA <- pAA*(1-a)/(1-k)
    gAa <- pAa*(1-a*b)/(1-k)
    gaa <- paa*(1-a*r)/(1-k)
    econtrol<-c(gAA,gAa,gaa)*nc

    if (length(econtrol[econtrol<=0])>0)
      {
      #   cat("b=", b, ", r=",  r,", q=", q, ", a=", a, "       2   \n")
      ll <- 1.0e10
      return(2*ll)
      }

    ll <- ll + sum(control * log(control/econtrol))

    return(2*ll)
}

## deviance, H_a (DHW) and 1 parameter: gamma, and
## beta=1 fixed, recessive model  
# case and control data in order of AA,Aa,aa
# parms=[p0,p1,gamma], logit or log transformed

DevHaGrecessive.est <- function(case,ctl,k0,initial) {

   fitout<-nlm(DevHaGrecessive,initial,case=case,control=ctl,k=k0,hessian=TRUE)

         hessinv<-MASS::ginv(fitout$hessian)     # generalized inverse of Hessian

#### beta fixed at 1
          beta <- 1 
          
### p0
          p0 <- invlogit(fitout$estimate[1])
          se.p0 <- se.invlogit(fitout$estimate[1],sqrt(hessinv[1,1]))

### p1

          p1 <- invlogit(fitout$estimate[2])
          se.p1 <- se.invlogit(fitout$estimate[2],sqrt(hessinv[2,2]))

          p <- p0+p1/2
           ## covariance of p0 and p1
          a12 <- cov.invlogit (fitout$estimate[1], fitout$estimate[2],
                                  hessinv[1,2])
          se.p <- sqrt(se.p0^2+se.p1^2/4+a12)

          q <- 1-p
          se.q <- se.p

### gamma
          gamma <- exp(fitout$estimate[3])
          se.g <- se.exp(fitout$estimate[3], sqrt(hessinv[3,3]))

### imbreed coefficient f  estimate
   #      f <- fitout$estimate[3]
   #      f <- exp(f)/(1+exp(f)) 
   #      se.f <- sqrt(hessinv[3,3])

### alpha
          pAA <- p0
          pAa <- p1
          paa <- 1-p0-p1

          alpha <- k0/( pAA + pAa*beta + paa*gamma )

### Deviance
          deviance<-fitout$minimum[1]

    list(alpha=alpha, beta=c(beta,0), gamma=c(gamma,se.g),
         q=c(q,se.q),p0=c(p0,se.p0), p1=c(p1,se.p1),
         deviance=c(deviance, 1-pchisq(deviance,1)) ) 
}


########################################################################


#######################################################
# dominant models 
#######################################################


######################################################################
## deviance, H0 and 1 parameters: beta=gamma  
# dominant model
# case and control data in order of AA,Aa,aa
# parms=[p,gamma], logit and log transformed 


DevH0dominant <- function(parms, case, control, k)
{
    p<-parms[1]
    r<-parms[2]

    p<-pfix(p,-15,15)
    p <- exp(p)/(1+exp(p))

    r<-pfix(r,-10,10)
    r <- exp(r)

    b <- r     # beta=gamma

    q<-1-p
        
    a<-k/( (1-q)^2+2*q*(1-q)*b+q^2*r )
    
    nt<-sum(case)
    nc<-sum(control)

    gAA<-(1-q)^2*a/k
    gAa<-2*q*(1-q)*a*b/k
    gaa<-q^2*a*r/k
    ecase<-c(gAA,gAa,gaa)*nt

    if ( length(ecase[ecase<=0])>0 )
      {
      #   cat(lcase, "- -", ecase, "       1    \n")
       ll <- 1.0e10
       return(2*ll)
       }

    ll <-  sum(case * log(case/ecase))
    
    gAA<-(1-q)^2*(1-a)/(1-k)
    gAa<-2*q*(1-q)*(1-a*b)/(1-k)
    gaa<-q^2*(1-a*r)/(1-k)
    econtrol<-c(gAA,gAa,gaa)*nc

    if (length(econtrol[econtrol<=0])>0)
      {
      #   cat("b=", b, ", r=",  r,", q=", q, ", a=", a, "       2   \n")
      ll <- 1.0e10
      return(2*ll)
      }

    ll <- ll + sum(control * log(control/econtrol))

    return(2*ll)
}

DevH0dominant.est <- function(case,ctl,k0,initial) {

   fitout<-nlm(DevH0dominant,initial,case=case,control=ctl,k=k0,hessian=TRUE)

         hessinv<-MASS::ginv(fitout$hessian)     # generalized inverse of Hessian

### p
          p <- invlogit(fitout$estimate[1])
          se.p <- se.invlogit(fitout$estimate[1],sqrt(hessinv[1,1]))

          q <- 1-p
          se.q <- se.p

### gamma
          gamma <- exp(fitout$estimate[2])
          se.g <- se.exp(fitout$estimate[2],sqrt(hessinv[2,2]))

#### beta=gamma
          beta <- gamma
          se.b <- se.g 

### alpha
          alpha<- alpha.f(beta, gamma, q, k0)

### Deviance
          deviance<-fitout$minimum[1]

    list(alpha=alpha, beta=c(beta,se.g), gamma=c(gamma,se.g),
         q=c(q,se.q),deviance=c(deviance, 1-pchisq(deviance,2)) ) 
}


##########################################################################

## deviance, H_a (DHW) and 1 parameter: beta=gamma, 
## dominant model  
# case and control data in order of AA,Aa,aa
# parms=[p0,p1,gamma], logit or log transformed

DevHaGdominant <- function(parms, case, control, k)
{
    p0 <- parms[1]
    p1 <- parms[2]
    r <- parms[3]

    p0 <- pfix(p0,-15,15)
    p0 <- exp(p0)/(1+exp(p0))

    p1 <- pfix(p1,-15,15)
    p1 <- exp(p1)/(1+exp(p1))

    r<- pfix(r,-10,10)
    r <- exp(r)

    b <- r         # beta=gamma
 
    pAA <- p0
    pAa <- p1
    paa <- 1-p0-p1

    a <- k/( pAA+ pAa*b+ paa*r )
    
    nt<-sum(case)
    nc<-sum(control)

    gAA <- pAA*a/k
    gAa <- pAa*a*b/k
    gaa <- paa*a*r/k
    ecase<-c(gAA,gAa,gaa)*nt

    if ( length(ecase[ecase<=0])>0 )
      {
      #   cat(lcase, "- -", ecase, "       1    \n")
       ll <- 1.0e10
       return(2*ll)
       }

    ll <-  sum(case * log(case/ecase))
    
    gAA <- pAA*(1-a)/(1-k)
    gAa <- pAa*(1-a*b)/(1-k)
    gaa <- paa*(1-a*r)/(1-k)
    econtrol<-c(gAA,gAa,gaa)*nc

    if (length(econtrol[econtrol<=0])>0)
      {
      #   cat("b=", b, ", r=",  r,", q=", q, ", a=", a, "       2   \n")
      ll <- 1.0e10
      return(2*ll)
      }

    ll <- ll + sum(control * log(control/econtrol))

    return(2*ll)
}

## deviance, H_a (DHW) and 1 parameter: beta=gamma 
## dominant model  
# case and control data in order of AA,Aa,aa
# parms=[p0,p1,gamma], logit or log transformed

DevHaGdominant.est <- function(case,ctl,k0,initial) {

   fitout<-nlm(DevHaGdominant,initial,case=case,control=ctl,k=k0,hessian=TRUE)

         hessinv<-MASS::ginv(fitout$hessian)     # generalized inverse of Hessian

### p0
          p0 <- invlogit(fitout$estimate[1])
          se.p0 <- se.invlogit(fitout$estimate[1],sqrt(hessinv[1,1]))

### p1

          p1 <- invlogit(fitout$estimate[2])
          se.p1 <- se.invlogit(fitout$estimate[2],sqrt(hessinv[2,2]))

          p <- p0+p1/2
           ## covariance of p0 and p1
          a12 <- cov.invlogit (fitout$estimate[1], fitout$estimate[2],
                                  hessinv[1,2])
          se.p <- sqrt(se.p0^2+se.p1^2/4+a12)

          q <- 1-p
          se.q <- se.p

### gamma
          gamma <- exp(fitout$estimate[3])
          se.g <- se.exp(fitout$estimate[3], sqrt(hessinv[3,3]))

#### beta=gamma 
          beta <- gamma  
          se.b <- se.g

### alpha
          pAA <- p0
          pAa <- p1
          paa <- 1-p0-p1

          alpha <- k0/( pAA + pAa*beta + paa*gamma )

### Deviance
          deviance<-fitout$minimum[1]

    list(alpha=alpha, beta=c(beta,se.b), gamma=c(gamma,se.g),
         q=c(q,se.q),p0=c(p0,se.p0), p1=c(p1,se.p1),
         deviance=c(deviance, 1-pchisq(deviance,1)) ) 
}

#' A likelihood ratio test of population Hardy-Weinberg equilibrium for case-control studies
#'
#' A likelihood ratio test of population Hardy-Weinberg equilibrium for case-control studies
#'
#' @param model model specification, dominant, recessive.
#' @param case a vector of genotype counts in cases.
#' @param ctrl a vector of genotype counts in controls.
#' @param k0 prevalence of disease in the population.
#' @param initial1 initial values for beta, gamma, and q.
#' @param initial2 initial values for logit(p) and log(gamma).
#'
#' @details
#' This is a collection of utility functions. The null hypothesis declares that the proportions of
#' genotypes are according to Hardy-Weinberg law, while under the alternative hypothesis, the expected
#' genotype counts are according to the probabilities that particular genotypes are obtained conditional
#' on the prevalence of disease in the population. In so doing, Hardy-Weinberg equilibrium is considered
#' using both case and control samples but pending on the disease model such that 2-parameter multiplicative
#' model is built on baseline genotype \eqn{\alpha}{alpha}, \eqn{\alpha\beta}{alpha beta} and
#' \eqn{\alpha\gamma}{alpha gamma}.
#'
#' @export
#' @return
#' The returned value is a list with the following components.
#' \describe{
#' \item{Cox}{statistics under a general model}
#' \item{t2par}{under the null hypothesis}
#' \item{t3par}{under the alternative hypothesis}
#' \item{lrt.stat}{the log-likelihood ratio statistic}
#' \item{pval}{the corresponding p value}
#' }
#'
#' @references
#' Yu C, Zhang S, Zhou C, Sile S. A likelihood ratio test of population Hardy-Weinberg equilibrium
#' for case-control studies. Genetic Epidemiology 33:275-280, 2009
#'
#' @seealso \code{\link[gap]{hwe}}
#'
#' @examples
#' \dontrun{
#'
#' ### Saba Sile, email of Jan 26, 2007, data always in order of GG AG AA, p=Pr(G),
#' ### q=1-p=Pr(A)
#' case=c(155,27,4)
#' ctrl=c(408,55,15)
#' k0=.2
#' initial1=c(1.0,0.94,0.0904)
#' initial2=c(logit(1-0.0904),log(0.94))
#' hwe.cc("recessive",case,ctrl,k0, initial1, initial2)
#'
#' ### John Phillips III, TGFb1 data codon 10: TT CT CC, CC is abnormal and increasing
#' ### TGFb1 activity
#' case=c(29,78,13)
#' ctrl=c(17,28,6)
#' k0 <- 1e-5
#' initial1 <- c(2.45,2.45,0.34)
#' initial2 <- c(logit(1-0.34),log(2.45))
#' hwe.cc("dominant",case,ctrl,k0,initial1,initial2)
#' }
#'
#' @author Chang Yu, Li Wang, Jing Hua Zhao
#' @keywords htest

hwe.cc <- function(model, case, ctrl, k0, initial1, initial2)
{
model.int <- charmatch(model,c("dominant","recessive"))
if(is.na(model.int)) stop("Invalide model specification")
if(model.int==0) stop("Ambiguous model specification")

Cox1 <- Cox.est(case,ctrl,k0,initial1) 
if (model.int==1) {
   t2par <- DevH0dominant.est(case,ctrl,k0,initial2)
   p0i <- (case[1]+ctrl[1])/(sum(case)+sum(ctrl))
   p1i <- (case[2]+ctrl[2])/(sum(case)+sum(ctrl))
   initial <- c(logit(p0i), logit(p1i), log(1.0753569)) # logit(p0), logit(p1), log(gamma)
   t3par <- DevHaGdominant.est(case,ctrl,k0,initial)
}
if (model.int==2) {
   t2par <- DevH0recessive.est(case,ctrl,k0,initial2)
   p0i <- (case[1]+ctrl[1])/(sum(case)+sum(ctrl))
   p1i <- (case[2]+ctrl[2])/(sum(case)+sum(ctrl))
   initial <- c(logit(p0i), logit(p1i), log(1.0753569)) # logit(p0), logit(p1), log(gamma) 
   t3par <- DevHaGrecessive.est(case,ctrl,k0,initial)
}
lrt1.stat <- t2par$deviance[1]-t3par$deviance[1]

list(Cox=Cox1, t2par=t2par, t3par=t3par, lrt.stat=lrt1.stat, pval=pchisq(lrt1.stat,1,lower.tail=FALSE))
}
