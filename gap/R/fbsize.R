#' Sample size for family-based linkage and association design
#'
#' This function implements Risch and Merikangas (1996) statistics 
#' evaluating power for family-based linkage (affected sib pairs, ASP) and
#' association design. They are potentially useful in the prospect of
#' genome-wide association studies.
#'
#' The function calls auxiliary functions sn() and strlen; sn() 
#' contains the necessary thresholds for power calculation while
#' strlen() evaluates length of a string (generic).
#'
#' @param gamma genotype relative risk assuming multiplicative model.
#' @param p frequency of disease allele.
#' @param alpha Type I error rates for ASP linkage, TDT and ASP-TDT.
#' @param beta Type II error rate.
#' @param debug verbose output.
#' @param error 0=use the correct formula,1=the original paper.
#'
#' @return The returned value is a list containing:
#'
#' \describe{
#'   \item{gamma}{input gamma.}
#'   \item{p}{input p.}
#'   \item{n1}{sample size for ASP.}
#'   \item{n2}{sample size for TDT.}
#'   \item{n3}{sample size for ASP-TDT.}
#'   \item{lambdao}{lambda o.}
#'   \item{lambdas}{lambda s.}
#' }
#'
#' @references
#'
#' Risch, N. and K. Merikangas (1996). The future of genetic studies of
#' complex human diseases. Science 273(September): 1516-1517.
#'
#' Risch, N. and K. Merikangas (1997). Reply to Scott el al. Science
#' 275(February): 1329-1330.
#'
#' Scott, W. K., M. A. Pericak-Vance, et al. (1997). Genetic analysis of 
#' complex diseases. Science 275: 1327.
#'
#' @seealso \code{\link[gap]{pbsize}}
#'
#' @examples
#' models <- matrix(c(
#'    4.0, 0.01,
#'    4.0, 0.10,
#'    4.0, 0.50, 
#'    4.0, 0.80,
#'    2.0, 0.01,
#'    2.0, 0.10,
#'    2.0, 0.50,
#'    2.0, 0.80,
#'    1.5, 0.01,    
#'    1.5, 0.10,
#'    1.5, 0.50,
#'    1.5, 0.80), ncol=2, byrow=TRUE)
#' outfile <- "fbsize.txt"
#' cat("gamma","p","Y","N_asp","P_A","H1","N_tdt","H2","N_asp/tdt","L_o","L_s\n",
#'     file=outfile,sep="\t")
#' for(i in 1:12) {
#'   g <- models[i,1]
#'   p <- models[i,2]
#'   z <- fbsize(g,p)
#'   cat(z$gamma,z$p,z$y,z$n1,z$pA,z$h1,z$n2,z$h2,z$n3,z$lambdao,z$lambdas,file=outfile,
#'       append=TRUE,sep="\t")
#'   cat("\n",file=outfile,append=TRUE)
#' }
#' table1 <- read.table(outfile,header=TRUE,sep="\t")
#' nc <- c(4,7,9)
#' table1[,nc] <- ceiling(table1[,nc])
#' dc <- c(3,5,6,8,10,11)
#' table1[,dc] <- round(table1[,dc],2)
#' unlink(outfile)
#' # APOE-4, Scott WK, Pericak-Vance, MA & Haines JL
#' # Genetic analysis of complex diseases 1327
#' g <- 4.5
#' p <- 0.15
#' cat("\nAlzheimer's:\n\n")
#' fbsize(g,p)
#' # note to replicate the Table we need set alpha=9.961139e-05,4.910638e-08 and
#' # beta=0.2004542 or reset the quantiles in fbsize.R
#'
#' @author Jing Hua Zhao
#' @note extracted from rm.c.
#' @keywords misc

fbsize <- function (gamma,p,alpha=c(1e-4,1e-8,1e-8),beta=0.2,debug=0,error=0)
# Family-based sample sizes
# Jing Hua Zhao 30-12-98, 19-8-2009
# Risch & Merikangas 1996
# Science 273: 1516-17 13SEP1996
# Science 275: 1327-30 28FEB1997
{
  sn <- function (all,alpha,beta,op)
  # m=0,v=1 under the null hypotheses
  # to be used by fbsize()
  {
    m <- all[1]
    v <- all[2]
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
  if (error==1) aa.v <- 0
  else aa.v <- 4*y*(1-y)
  aa <- c(aa.m,aa.v)

  n1 <- sn(aa,alpha[1],beta,1)

# TDT
  aa.m <- sqrt(h)*(gamma-1)/(gamma+1)
  aa.v <- 1-h*((gamma-1)/(gamma+1))^2
  aa <- c(aa.m,aa.v)

  n2 <- sn(aa,alpha[2],beta,2)

# ASP-TDT
  h <- h2 <- p*q*(gamma+1)^2/(2*(p*gamma+q)^2+p*q*(gamma-1)^2)
  aa.m <- sqrt(h)*(gamma-1)/(gamma+1)
  aa.v <- 1-h*((gamma-1)/(gamma+1))^2
  aa <- c(aa.m,aa.v)

  n3 <- sn(aa,alpha[3],beta,3)

  if (debug==1)
  {
     cat("K=",k, "VA=",va, "VD=",vd,"\n")
     cat(format(gamma,width=4,nsmall=2),
         format(p,width=5,nsmall=2),
         format(round(y,digits=3),nsmall=3),
         rep("",10-nchar(ceiling(n1))),ceiling(n1),
         format(round(pA,digits=3),nsmall=3),
         format(round(h1,digits=3),nsmall=3),
         rep("",8-nchar(ceiling(n2))),ceiling(n2),
         format(round(h2,digits=3),nsmall=3),
         rep("",8-nchar(ceiling(n3))),ceiling(n3),
         format(round(lambdao,digits=2),nsmall=2),
         format(round(lambdas,digits=2),nsmall=2),"\n")
  }
  list(gamma=gamma,p=p,y=y,n1=n1,pA=pA,h1=h1,n2=n2,h2=h2,n3=n3,
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
