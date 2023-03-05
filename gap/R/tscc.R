#' A normal z-test of two proportions
#' @noRd

z <- function(p1,p2,n1,n2,r)
{
   z.mean <- p1-p2
   z.var <- p1*(1-p1)/n1+p2*(1-p2)/n2
   invisible(z.mean/sqrt(z.var/(2*r)))
}

#' A function used by tscc
#' @noRd

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

#' Power calculation for two-stage case-control design
#'
#' @param model any in c("multiplicative","additive","dominant","recessive").
#' @param GRR genotype relative risk.
#' @param p1 the estimated risk allele frequency in cases.
#' @param n1 total number of cases.
#' @param n2 total number of controls.
#' @param M total number of markers.
#' @param alpha.genome false positive rate at genome level.
#' @param pi.samples sample% to be genotyped at stage 1.
#' @param pi.markers markers% to be selected (also used as the false positive rate at stage 1).
#' @param K the population prevalence.
#'
#' @details
#' This function gives power estimates for two-stage case-control design for genetic association.
#'
#' The false positive rates are calculated as follows,
#'
#' \deqn{P(|z1|>C1)P(|z2|>C2,sign(z1)=sign(z2))} and \deqn{P(|z1|>C1)P(|zj|>Cj||z1|>C1)}
#' for replication-based and joint analyses, respectively; where C1, C2, and Cj
#' are threshoulds at stages 1, 2 replication and joint analysis, 
#'
#' \deqn{z1 = z(p1,p2,n1,n2,pi.samples)}
#' \deqn{z2 = z(p1,p2,n1,n2,1-pi.samples)}
#' \deqn{zj = sqrt(pi.samples)*z1+sqrt(1-pi.samples)*z2}
#'
#' @export
#' @return
#' The returned value is a list containing a copy of the input plus output as follows,
#' - model any in c("multiplicative","additive","dominant","recessive").
#' - GRR genotype relative risk.
#' - p1 the estimated risk allele frequency in cases.
#' - pprime expected risk allele frequency in cases.
#' - p expected risk allele frequency in controls.
#' - n1 total number of cases.
#' - n2 total number of controls.
#' - M total number of markers.
#' - alpha.genome false positive rate at genome level.
#' - pi.samples sample% to be genotyped at stage 1.
#' - pi.markers markers% to be selected (also used as the false positive rate at stage 1).
#' - K the population prevalence.
#' - C threshoulds for no stage, stage 1, stage 2, joint analysis.
#' - power power corresponding to C.
#'
#' @references
#' \insertRef{skol06}{gap}
#'
#' @examples
#' \dontrun{
#' K <- 0.1
#' p1 <- 0.4
#' n1 <- 1000
#' n2 <- 1000 
#' M <- 300000
#' alpha.genome <- 0.05
#' GRR <- 1.4
#' p1 <- 0.4
#' pi.samples <- 0.2
#' pi.markers <- 0.1
#'
#' options(echo=FALSE)
#' cat("sample%,marker%,GRR,(thresholds x 4)(power estimates x 4)","\n")
#' for(GRR in c(1.3,1.35,1.40))
#' {
#'    cat("\n")
#'    for(pi.samples in c(1.0,0.5,0.4,0.3,0.2))
#'    {
#'       if(pi.samples==1.0) s <- 1.0
#'       else s <- c(0.1,0.05,0.01)
#'       for(pi.markers in s)
#'       {
#'         x <- tscc("multiplicative",GRR,p1,n1,n2,M,alpha.genome,
#'                   pi.samples,pi.markers,K)
#'         l <- c(pi.samples,pi.markers,GRR,x$C,x$power)
#'         l <- sprintf("%.2f %.2f %.2f, %.2f %.2f %.2f %.2f, %.2f %.2f %.2f %.2f",
#'                      l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11])
#'         cat(l,"\n")
#'       }
#'       cat("\n")
#'    }
#' }
#' options(echo=TRUE)
#' }
#' @author Jing Hua Zhao
#' @note `solve.skol` is adapted from CaTS.
#' @keywords misc

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
