# 7-3-2008 MRC-Epid JHZ

#' Meta-analysis of p values
#'
#' @param data data frame.
#' @param N Number of studies.
#' @param verbose Control of detailed output.
#' @param prefixp Prefix of p value, with default value "p".
#' @param prefixn Preifx of sample size, with default value "n".
#'
#' @details
#' This function is the method of meta-analysis used in the Genetic Investigation
#' of ANThropometric Traits (GIANT) consortium, which is based on normal approximation
#' of p values and weighted by sample sizes from individual studies.
#'
#' @export
#' @return
#' - x2 Fisher's chi-squared statistics.
#' - p P values from Fisher's method according to chi-squared distribution with 2*N degree(s) of freedom.
#' - z Combined z value.
#' - p1 One-sided p value.
#' - p2 Two-sided p value.
#'
#' @examples
#' \dontrun{
#' s <- data.frame(p1=0.1^rep(8:2,each=7,times=1),n1=rep(32000,49),
#'                 p2=0.1^rep(8:2,each=1,times=7),n2=rep(8000,49))
#' cbind(s,metap(s,2))
#'
#' # Speliotes, Elizabeth K., M.D. [ESPELIOTES@PARTNERS.ORG]
#' # 22-2-2008 MRC-Epid JHZ
#'
#' np <- 7
#' p <- 0.1^((np+1):2)
#' z <- qnorm(1-p/2)
#' n <- c(32000,8000)
#' n1 <- n[1]
#'
#' s1 <- s2 <- vector("numeric")
#'
#' for (i in 1:np)
#' {
#'    a <- z[i]
#'    for (j in 1:np)
#'    {
#'        b <- z[j]
#'        metaz1 <- (sqrt(n1)*a+sqrt(n[1])*b)/sqrt(n1+n[1])
#'        metap1 <- pnorm(-abs(metaz1))
#'        metaz2 <- (sqrt(n1)*a+sqrt(n[2])*b)/sqrt(n1+n[2])
#'        metap2 <- pnorm(-abs(metaz2))
#'        k <- (i-1)*np+j
#'        cat(k,"\t",p[i],"\t",p[j],"\t",metap1,metaz1,"\t",metap2,metaz2,"\n")
#'        s1[k] <- metap1
#'        s2[k] <- metap2
#'   }
#' }
#'
#' q <- -log10(sort(p,decreasing=TRUE))
#' t1 <- matrix(-log10(sort(s1,decreasing=TRUE)),np,np)
#' t2 <- matrix(-log10(sort(s2,decreasing=TRUE)),np,np)
#'
#' par(mfrow=c(1,2),bg="white",mar=c(4.2,3.8,0.2,0.2))
#' persp(q,q,t1)
#' persp(q,q,t2)
#' }
#'
#' @seealso [`metareg`]
#'
#' @author Jing Hua Zhao
#' @keywords htest

metap <- function(data,N,verbose="Y",prefixp="p",prefixn="n")
{
   M <- dim(data)[1]
   rawp <- rawn <- w <- matrix("numeric",M,N)
   x2 <- p <- metaz <- vector("numeric",M)
   for (j in 1:M)
   {
       for (i in 1:N)
       {
           rawp[j,i] <- data[paste(prefixp,i,sep="")][j,]
           rawn[j,i] <- data[paste(prefixn,i,sep="")][j,]
       }
       s1 <- rawn[j,]
       s2 <- as.numeric(s1)
       w[j,] <- sqrt(s2)/sqrt(sum(s2))
       t1 <- rawp[j,]
       t2 <- as.numeric(t1)
       x2[j] <- 2*sum(-log(t2))
       metaz[j] <- sum(as.numeric(w[j,])*qnorm(1-t2/2))
   }
   p <- pchisq(x2,2*N,lower.tail=FALSE)
   metap1 <- pnorm(-abs(metaz))
   metap2 <- 2*metap1
   if (toupper(verbose)=="Y")
   cat("\n\nFisher's method x2=",x2,"\ndf=",2*N,"\ntwo-sided p=",p,"\n\n")
   cat("\n\nCombined z=",metaz,"\none-sided p=",metap1,"\ntwo-sided p=",metap2,"\n\n")
   invisible(data.frame(x2=x2,p=p,z=metaz,p1=metap1,p2=metap2))
}
