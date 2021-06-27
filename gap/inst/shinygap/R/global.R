library(dplyr)
library(ggplot2)
library(htmltools)
library(rmarkdown)
library(plotly)
library(vroom)
library(shiny)
library(shinydashboard)

fbsize <- function (gamma,p,alpha=1e-4,beta=0.2,error=FALSE)
# Family-based sample sizes
# Jing Hua Zhao 30-12-98, 19-8-2009, 13-6-2021
# Risch & Merikangas 1996
# Science 273: 1516-17 13SEP1996
# Science 275: 1327-30 28FEB1997
{
  sn <- function (all,alpha,beta,op)
  # m=0,v=1 under the null hypotheses
  {
    m <- all[[1]]
    v <- all[[2]]
    z1beta <- qnorm(beta)                # -0.84, 1-beta=0.8 (-.84162123)
    zalpha <- -qnorm(alpha)              #  3.72, alpha=1E-4 (3.7190165); 5.33, alpha=5E-8 (5.3267239)
    s <- ((zalpha-sqrt(v)*z1beta)/m)^2/2 # shared/transmitted for each parent
    if(op==3) s <- s/2                   # the sample size is halved
    s
  }

  q <- 1-p
  k <- (p*gamma+q)^2
  va <- 2*p*q*((gamma-1)*(p*gamma+q))^2
  vd <- (p*q*(gamma-1)^2)^2

  w <- p*q*(gamma-1)^2/(p*gamma+q)^2
  y <- (1+w)/(2+w)
  lambdas <- (1+0.5*w)^2
  lambdao <- 1+w
  h <- h1 <- p*q*(gamma+1)/(p*gamma+q)
  pA <- gamma/(gamma+1)

# ASP
  nl.m <- 0
  nl.v <- 1
  aa.m <- 2*y-1
  if (error) aa.v <- 0
  else aa.v <- 4*y*(1-y)
  aa <- list(aa.m,aa.v)

  n1 <- sn(aa,alpha,beta,1)

# TDT
  aa.m <- sqrt(h)*(gamma-1)/(gamma+1)
  aa.v <- 1-h*((gamma-1)/(gamma+1))^2
  aa <- list(aa.m,aa.v)

  n2 <- sn(aa,alpha,beta,2)

# ASP-TDT
  h <- h2 <- p*q*(gamma+1)^2/(2*(p*gamma+q)^2+p*q*(gamma-1)^2)
  aa.m <- sqrt(h)*(gamma-1)/(gamma+1)
  aa.v <- 1-h*((gamma-1)/(gamma+1))^2
  aa <- list(aa.m,aa.v)

  n3 <- sn(aa,alpha,beta,3)

  list(gamma=gamma,p=p,y=y,n1=ceiling(n1),pA=pA,h1=h1,n2=ceiling(n2),h2=h2,n3=ceiling(n3),
      lambdao=lambdao,lambdas=lambdas)

}

#MSmodel <- function (lambdas, dom=0)
# Ebers G et al (1996) on multiple sclerosis Nat Genet 13:472
# for modest lambdas, two models are the same
# lambdas=sibling risk ratio associated with the locus
#{
#  if (dom!=0) {
#     # model with moderate amount of dominance variance
#     y <- 1-0.5/sqrt(lambdas)
#     z2 <- y*y
#     z1 <- 2*y*(1-y)
#     z0 <- (1-y)^2
#  }
#  else {
#    # additive model, dominance variance eq 0#
#    z2 <- 0.5-0.25/lambdas
#    z1 <- 0.5
#    z0 <- 0.25/lambdas
#  }
#  z <- c(z2,z1,z0)
#  z
#}

pbsize <- function (kp, g=4.5, p=0.15, alpha=5e-8, beta=0.2)
# population-based sample size
# alpha=5e-8, beta=0.8, alpha would give 5% genome-wide significance level
# x2alpha = 29.72 (Q=29.7168)
#
# lambda is the NCP from the marginal table
# pi is the pr(Affected|aa)
{
  z1alpha <- qnorm(1-alpha/2)
  zbeta <- qnorm(beta)
  q <- 1-p
  pi <- kp/(g*p+q)^2
  lambda <- pi*p*q*(g-1)^2/(1-kp)
  n <- (z1alpha-zbeta)^2/lambda
  n
}

ccsize <- function(n,q,pD,p1,theta,alpha,beta=0.2,power=FALSE,verbose=FALSE)
{
   p2 <- 1 - p1
   if (power)
   {
   # equation (5)
     z_alpha <- qnorm(alpha)
     z <- z_alpha + sqrt(n) * theta * sqrt(p1 * p2 / (1 / pD + (1 / q - 1)))
     invisible(pnorm(z))
   } else {
   # equation (6)
     z_alpha <- qnorm(alpha, lower.tail=FALSE)
     z_beta <- qnorm(beta, lower.tail=FALSE)
     theta_lon <- (z_alpha + z_beta) / sqrt(p1 * p2 * pD)
     d <- (theta / theta_lon)^2 - (1 - pD) / n
     nb <- ceiling(pD / d)
     nb [nb > n] <- -999
     if (any(d <= 0) & verbose) cat("bad hazard ratio =", exp(theta), "\n")
     else if (any(nb > n) & verbose) cat("bad subcohort size", nb, "\n")
     invisible(nb)
   }
}

# 1-3-2008, MRC-Epid, JHZ

KCC <- function(model,GRR,p1,K)
# 6-6-2018
{
   model.idx <- charmatch(model,c("multiplicative","additive","recessive","dominant","overdominant"))
   if(is.na(model.idx)) stop("Invalid model type")
   if(model.idx == 0) stop("Ambiguous model type")
   multiplicative <- c(1,GRR,GRR*GRR)
   additive <- c(1,GRR,2*GRR-1)
   recessive <- c(1,1,GRR)
   dominant <- c(1,GRR,GRR)
   overdominant <- c(GRR,1,GRR)
   f <- switch(model.idx,multiplicative,additive,recessive,dominant,overdominant)
   scale <- K/(f[1]*(1-p1)^2+f[2]*2*p1*(1-p1)+f[3]*p1^2)
   f <- f*scale
#  if(f[3]>1) stop("misspecified model")
   pprime <- (f[3]*p1^2+f[2]*p1*(1-p1))/K
   p <- ((1-f[3])*p1^2+(1-f[2])*p1*(1-p1))/(1-K)
   invisible(list(pprime=pprime,p=p))
}

z <- function(p1,p2,n1,n2,r)
{
   z.mean <- p1-p2
   z.var <- p1*(1-p1)/n1+p2*(1-p2)/n2
   invisible(z.mean/sqrt(z.var/(2*r)))
}

solve_skol <- function(rootfun,target,lo,hi,e)
{
   if(rootfun(lo)>rootfun(hi)) {
      temp <- lo
      lo <- hi
      hi <- temp
   }
   e <- e*target+1e-20
   while(1) {
      d <- hi-lo
      point <- lo+d/2
      fpoint <- rootfun(point)
      if (fpoint<target) {
         d <- lo-point
         lo <- point
      } else {
         d <- point-hi
         hi <- point
      }
      if(abs(d)<e|fpoint==target) break
   }
   point
}

tscc <- function(model,GRR,p1,n1,n2,M,alpha.genome,pi.samples,pi.markers,K)
{
   ce <- environment()
   l <- c(p1,pi.samples,pi.markers,K)
   if(any(sapply(l,">",1))|any(sapply(c(l,GRR,n1,n2,M),"<=",0))) stop("invalid input")
   x <- KCC(model,GRR,p1,K)
   pprime <- x$pprime
   p <- x$p
   alpha.marker <- alpha.genome/M
   C <- qnorm(1-alpha.marker/2)
   m <- z(pprime,p,n1,n2,1)
   power <- 1-pnorm(C-m)+pnorm(-m-C)
   m1 <- z(pprime,p,n1,n2,pi.samples)
   C1 <- qnorm(1-pi.markers/2)
   A <- 1-pnorm(C1-m1)
   B <- pnorm(-C1-m1)
   P1 <- A+B
   power1 <- A
   m2 <- z(pprime,p,n1,n2,(1-pi.samples))
   C2 <- qnorm(1-alpha.marker/pi.markers)
   P2 <- (1-pnorm(C2-m2))*A/(A+B)+pnorm(-C2-m2)*B/(A+B)
   power2 <- P1*P2
   u <- function(z1) {
      if(rootfinding) m1 <- m2 <- 0
      m <- m2*sqrt(1-pi.samples)+sqrt(pi.samples)*z1
      v <- 1-pi.samples
      u1 <- (-Cj-m)/sqrt(v)
      u2 <- ( Cj-m)/sqrt(v)
      (pnorm(u1)+1-pnorm(u2))/sqrt(2*pi)*exp(-0.5*(z1-m1)^2)
   }
   rootfun <- function(x) {
       assign("Cj",x,envir=ce)
       integrate(u,-Inf,-C1)$value+integrate(u,C1,Inf)$value
   }
   Cj <- 0
   rootfinding <- TRUE
   Cj <- solve_skol(rootfun,alpha.marker,C2,C,1e-6)
   rootfinding <- FALSE
   powerj <- integrate(u,-Inf,-C1)$value+integrate(u,C1,Inf)$value
   invisible(list(model=model,GRR=GRR,p1=p1,pprime=pprime,p=p,n1=n1,n2=n2,M=M,
         pi.samples=pi.samples,pi.markers=pi.markers,alpha.genome=alpha.genome,
         C=c(C,C1,C2,Cj),power=c(power,power1,power2,powerj),K=K))
}
