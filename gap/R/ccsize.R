#' Power and sample size for case-cohort design
#'
#' @param n the total number of subjects in the cohort.
#' @param q the sampling fraction of the subcohort.
#' @param pD the proportion of the failures in the full cohort.
#' @param p1 proportions of the two groups (p2=1-p1).
#' @param theta log-hazard ratio for two groups.
#' @param alpha type I error -- significant level.
#' @param beta type II error.
#' @param power if specified, the power for which sample size is calculated.
#' @param verbose error messages are explicitly printed out.
#'
#' @details
#' The power of the test is according to 
#' \deqn{\Phi\left(Z_\alpha+m^{1/2}\theta\sqrt{\frac{p_1p_2p_D}{q+(1-q)p_D}}\right)}{Phi(Z_alpha+m^0.5*theta*sqrt(p_1p_2p_D/q+(1-q)p_D))}
#' where \eqn{\alpha}{alpha} is the significance level, \eqn{\theta}{theta} is the log-hazard ratio for two groups, \eqn{p_j}{p_j}, 
#' j=1, 2, are the proportion of the two groups in the population. \eqn{m} is the total number of subjects in the subcohort, 
#' \eqn{p_D} is the proportion of the failures in the full cohort, and \eqn{q} is the sampling fraction of the subcohort.
#'
#' Alternatively, the sample size required for the subcohort is \deqn{m=nBp_D/(n-B(1-p_D))}{m=nBp_D/(n-B(1-p_D))}
#' where \eqn{B=(Z_{1-\alpha}+Z_\beta)^2/(\theta^2p_1p_2p_D)}{B=(Z_{1-alpha}+Z_beta)^2/(theta^2p_1p_2p_D))}, and \eqn{n} is the size of cohort.
#'
#' When infeaisble configurations are specified, a sample size of -999 is returned.
#'
#' @export
#' @return
#' The returned value is a value indicating the power or required sample size.
#'
#' @references
#' \insertRef{cai04}{gap}
#'
#' @seealso [`pbsize`]
#'
#' @examples
#' \dontrun{
#' # Table 1 of Cai & Zeng (2004).
#' outfile <- "table1.txt"
#' cat("n","pD","p1","theta","q","power\n",file=outfile,sep="\t")
#' alpha <- 0.05
#' n <- 1000
#' for(pD in c(0.10,0.05))
#' {
#'    for(p1 in c(0.3,0.5))
#'    {
#'       for(theta in c(0.5,1.0))
#'       {
#'          for(q in c(0.1,0.2))
#'          {
#'             power <- ccsize(n,q,pD,p1,alpha,theta)
#'             cat(n,"\t",pD,"\t",p1,"\t",theta,"\t",q,"\t",signif(power,3),"\n",
#'                 file=outfile,append=TRUE)
#'          }
#'       }
#'    }
#' }
#' n <- 5000
#' for(pD in c(0.05,0.01))
#' {
#'    for(p1 in c(0.3,0.5))
#'    {
#'       for(theta in c(0.5,1.0))
#'       {
#'          for(q in c(0.01,0.02))
#'          {
#'             power <- ccsize(n,q,pD,p1,alpha,theta)
#'             cat(n,"\t",pD,"\t",p1,"\t",theta,"\t",q,"\t",signif(power,3),"\n",
#'                 file=outfile,append=TRUE)
#'          }
#'       }
#'    }
#' }
#' table1<-read.table(outfile,header=TRUE,sep="\t")
#' unlink(outfile)
#' # ARIC study
#' outfile <- "aric.txt"
#' n <- 15792
#' pD <- 0.03
#' p1 <- 0.25
#' alpha <- 0.05
#' theta <- c(1.35,1.40,1.45)
#' beta1 <- 0.8
#' s_nb <- c(1463,722,468)
#' cat("n","pD","p1","hr","q","power","ssize\n",file=outfile,sep="\t")
#' for(i in 1:3)
#' {
#'   q <- s_nb[i]/n
#'   power <- ccsize(n,q,pD,p1,alpha,log(theta[i]))
#'   ssize <- ccsize(n,q,pD,p1,alpha,log(theta[i]),beta1)
#'   cat(n,"\t",pD,"\t",p1,"\t",theta[i],"\t",q,"\t",signif(power,3),"\t",ssize,"\n",
#'       file=outfile,append=TRUE)
#' }
#' aric<-read.table(outfile,header=TRUE,sep="\t")
#' unlink(outfile)
#' # EPIC study
#' outfile <- "epic.txt"
#' n <- 25000
#' alpha <- 0.00000005
#' power <- 0.8
#' s_pD <- c(0.3,0.2,0.1,0.05)
#' s_p1 <- seq(0.1,0.5,by=0.1)
#' s_hr <- seq(1.1,1.4,by=0.1)
#' cat("n","pD","p1","hr","alpha","ssize\n",file=outfile,sep="\t")
#' # direct calculation
#' for(pD in s_pD)
#' {
#'    for(p1 in s_p1)
#'    {
#'       for(hr in s_hr)
#'       {
#'          ssize <- ccsize(n,q,pD,p1,alpha,log(hr),power)
#'          if (ssize>0) cat(n,"\t",pD,"\t",p1,"\t",hr,"\t",alpha,"\t",ssize,"\n",
#'                           file=outfile,append=TRUE)
#'       }
#'    }
#' }
#' epic<-read.table(outfile,header=TRUE,sep="\t")
#' unlink(outfile)
#' # exhaustive search
#' outfile <- "search.txt"
#' s_q <- seq(0.01,0.5,by=0.01)
#' cat("n","pD","p1","hr","nq","alpha","power\n",file=outfile,sep="\t")
#' for(pD in s_pD)
#' {
#'    for(p1 in s_p1)
#'    {
#'       for(hr in s_hr)
#'       {
#'          for(q in s_q)
#'          {
#'             power <- ccsize(n,q,pD,p1,alpha,log(hr))
#'             cat(n,"\t",pD,"\t",p1,"\t",hr,"\t",q*n,"\t",alpha,"\t",power,"\n",
#'                 file=outfile,append=TRUE)
#'          }
#'       }
#'    }
#' }
#' search<-read.table(outfile,header=TRUE,sep="\t")
#' unlink(outfile)
#' }
#' @author Jing Hua Zhao
#' @note Programmed for EPIC study.
#' keywords misc

ccsize <- function(n,q,pD,p1,theta,alpha,beta=0.2,power=TRUE,verbose=FALSE)
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
