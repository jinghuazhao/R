#' LD statistics for two diallelic markers
#'
#' @param h a vector of haplotype frequencies.
#' @param n number of haplotypes.
#'
#' @details
#' It is possible to perform permutation test of \eqn{r^2} by re-ordering the genotype through
#' R's sample function, obtaining the haplotype frequencies by [`gc.em`]
#' or [`genecounting`], supplying the estimated haplotype frequencies to
#' the current function and record x2, and comparing the observed x2 and that from the
#' replicates.
#'
#' @export
#' @return
#' The returned value is a list containing:
#' - h the original haplotype frequency vector.
#' - n the number of haplotypes.
#' - D the linkage disequilibrium parameter.
#' - VarD the variance of D.
#' - Dmax the maximum of D.
#' - VarDmax the variance of Dmax.
#' - Dprime the scaled disequilibrium parameter.
#' - VarDprime the variance of Dprime.
#' - x2 the Chi-squared statistic.
#' - lor the log(OR) statistic.
#' - vlor the var(log(OR)) statistic.
#'
#' @references
#' \insertRef{zabetian03}{gap}
#'
#' \insertRef{zapata97}{gap}
#'
#' @seealso [`LDkl`]
#'
#' @examples
#' \dontrun{
#' h <- c(0.442356,0.291532,0.245794,0.020319)
#' n <- 481*2
#' t <- LD22(h,n)
#' }
#'
#' @author Jing Hua Zhao
#' @note extracted from 2ld.c, worked 28/6/03, tables are symmetric do not fix, see kbyl below
#' @keywords models

LD22<-function(h,n) 
{
   D<-VarD<-Dmax<-VarDmax<-Dprime<-VarDprime<-x2<-lor<-vlor<-00

   z<-.C("tbyt",h=as.vector(h), haplotypes=as.double(n),
          D=as.double(D), VarD=as.double(VarD),
          Dmax=as.double(Dmax), VarDmax=as.double(VarDmax),
          Dprime=as.double(Dprime), VarDprime=as.double(VarDprime),
          x2=as.double(x2),lor=as.double(lor),vlor=as.double(vlor),PACKAGE="gap")

    invisible(list(h=h,n=n,D=z$D,VarD=z$VarD,
         Dmax=z$Dmax,VarDmax=z$VarDmax,Dprime=z$Dprime,
         VarDprime=z$VarDprime,x2=z$x2,lor=z$lor,vlor=z$vlor))
}

# refine on 17/4/2005
# verbose, default values, etc.

#' LD statistics for two multiallelic markers
#'
#' LD statistics for two multiallelic loci. For two diallelic makers,
#' the familiar \eqn{r^2}{r^2} has standard error seX2.
#'
#' @param n1 number of alleles at marker 1.
#' @param n2 number of alleles at marker 2.
#' @param h a vector of haplotype frequencies.
#' @param n number of haplotypes.
#' @param optrho type of contingency table association, 0=Pearson, 1=Tschuprow, 2=Cramer (default).
#' @param verbose detailed output of individual statistics.
#'
#' @export
#' @return
#' The returned value is a list containing:
#' - n1 the number of alleles at marker 1.
#' - n2 the number of alleles at marker 2.
#' - h the haplotype frequency vector.
#' - n the number of haplotypes.
#' - Dp D'.
#' - VarDp variance of D'.
#' - Dijtable table of Dij.
#' - VarDijtable table of variances for Dij.
#' - Dmaxtable table of Dmax.
#' - Dijptable table of Dij'.
#' - VarDijptable table of variances for Dij'.
#' - X2table table of Chi-squares (based on Dij).
#' - ptable table of p values.
#' - x2 the Chi-squared statistic.
#' - seX2 the standard error of x2/n.
#' - rho the measure of association.
#' - seR the standard error of rho.
#' - optrho the method for calculating rho.
#' - klinfo the Kullback-Leibler information.
#'
#' @references
#' \insertRef{bishop75}{gap}
#'
#' \insertRef{cramer46}{gap}
#'
#' \insertRef{zapata01}{gap}
#'
#' \insertRef{zhao04}{gap}
#'
#' @seealso [`LD22`]
#'
#' @examples
#' \dontrun{
#' # two examples in the C program 2LD:
#' # two SNPs as in 2by2.dat
#' # this can be compared with output from LD22
#'
#' h <- c(0.442356,0.291532,0.245794,0.020319)
#' n <- 481*2
#' t <- LDkl(2,2,h,n)
#' t
#'
#' # two multiallelic markers as in kbyl.dat
#' # the two-locus haplotype vector is in file "kbyl.dat"
#' # The data is now available from 2ld in Haplotype-Analysis
#'
#' filespec <- system.file("kbyl.dat")
#' h <- scan(filespec,skip=1)
#' t <- LDkl(9,5,h,213*2,verbose=TRUE)
#' }
#'
#' @author Jing Hua Zhao
#' @note adapted from 2ld.c.
#' @keywords models

LDkl<-function(n1=2,n2=2,h,n,optrho=2,verbose=FALSE)
{
   Dp<-x2<-seX2<-rho<-seR<-klinfo<-0
   VarDp<-0
   Dijtable<-Dmaxtable<-Dijptable<-VarDijtable<-X2table<-VarDijptable<-matrix(rep(0,n1*n2),nrow=n1)

   z<-.C("kbyl",nalleles1=as.integer(n1), nalleles2=as.integer(n2),
          h=as.double(h), haplotypes=as.double(n),
          Dp=as.double(Dp),VarDp=as.double(VarDp),
          Dijtable=matrix(Dijtable,nrow=n1),
          VarDijtable=matrix(VarDijtable,nrow=n1),
          X2table=matrix(X2table,nrow=n1),
          Dmaxtable=matrix(Dmaxtable,nrow=n1),
          Dijptable=matrix(Dijptable,nrow=n1),
          VarDijptable=matrix(VarDijptable,nrow=n1),
          x2=as.double(x2), seX2=as.double(seX2),
          rho=as.double(rho), seR=as.double(seR), optrho=as.integer(optrho),
          klinfo=as.double(klinfo),verbose=as.integer(verbose),PACKAGE="gap")

   n1 <- z$nalleles1
   n2 <- z$nalleles2
   h <- t(matrix(z$h,nrow=n2))
   n <- z$haplotypes
   Dp <- z$Dp
   VarDp <- z$VarDp
   Dijtable <- t(matrix(z$Dijtable,nrow=n2))
   VarDijtable <- t(matrix(z$VarDijtable,nrow=n2))
   X2table <- t(matrix(z$X2table,nrow=n2))
   Dmaxtable <- t(matrix(z$Dmaxtable,nrow=n2))
   Dijptable <- t(matrix(z$Dijptable,nrow=n2))
   VarDijptable <- t(matrix(z$VarDijptable,nrow=n2))
   ptable <- 1-pchisq(X2table,1)
   x2 <- z$x2
   seX2 <- z$seX2
   rho <- z$rho
   seR <- z$seR
   optrho <- z$optrho
   klinfo <- z$klinfo
   df <- (n1-1)*(n2-1)
   if(verbose)
   {
      cat("\nEstimated haplotype frequencies\n\n")
      print(h)
      cat("\nTable of D\n\n")
      print(Dijtable)
      cat("\nTable of Dmax\n\n")
      print(Dmaxtable)
      cat("\nTable of D'\n\n")
      print(Dijptable)
      cat("\nTable of SE(D')\n\n")
      print(sqrt(VarDijptable))
      cat("\nTable of Chi-squares (based on D)\n\n")
      print(X2table)
      cat("\nTable of p values\n\n")
      print(ptable)
      cat("\nChi-squared statistic=",x2,"df=",df,"p=",1-pchisq(x2,df),"\n\n")
      cat("\nGlobal disequilibrium statistics and their standard errors\n\n")
      cat("D' coefficient=",Dp,'SE=',sqrt(VarDp),"\n\n")
      cat("Kullback-Leibler information",klinfo,"\n")
   }
   invisible(list(n1=n1, n2=n2, h=h, n=n,
   Dp=Dp,VarDp=VarDp,Dijtable=Dijtable, VarDijtable=VarDijtable, 
   Dmaxtable=Dmaxtable,
   Dijptable=Dijptable, VarDijptable=VarDijptable,
   X2table=X2table, ptable=ptable,
   x2=x2, seX2=seX2, rho=rho, seR=seR, optrho=optrho, klinfo=klinfo))
}

#' Haplotype frequency estimation based on a genotype table of two multiallelic markers
#'
#' Haplotype frequency estimation using expectation-maximization algorithm based on a table of genotypes of two multiallelic markers.
#'
#' @param obs a table of genotype counts.
#' @param k number of alleles at marker 1.
#' @param l number of alleles at marker 2.
#'
#' The dimension of the genotype table should be k*(k+1)/2 x l*(l+1)/2.
#'
#' Modified from 2ld.c.
#'
#' @export
#' @return
#' The returned value is a list containing:
#' - h haplotype Frequencies.
#' - l0 log-likelihood under linkage equilibrium.
#' - l1 log-likelihood under linkage disequilibrium.
#'
#' @seealso [`genecounting`]
#'
#' @examples
#' \dontrun{
#' # an example with known genotype counts 
#' z <- klem(obs=1:9)
#' # an example with imputed genotypes at SH2B1
#' source(file.path(find.package("gap"),"scripts","SH2B1.R"),echo=TRUE)
#' }
#'
#' @author Jing Hua Zhao
#' @keywords htest

klem <- function(obs,k=2,l=2)
{
  if(length(obs)!=k*l*(k+1)*(l+1)/4) stop("incorrect length of genotype table")
  h <- rep(0,k*l)
  l0 <- l1 <- 0
  z <- .C("kbylem",obs=as.double(obs),nalleles1=as.integer(k),nalleles2=as.integer(l),
        Rh=as.double(h),l0=as.double(l0),l1=as.double(l1))
  invisible(list(h=z$Rh,l0=z$l0,l1=z$l1))
}
