# 10/4/2005, 13/4/2005
gc.control <- function(xdata=FALSE, convll=1,handle.miss=0,eps=0.000001,
                       tol=0.00000001, maxit=50,pl=0.001,assignment="assign.dat",verbose=T)
{
   list(xdata=xdata,convll=convll,handle.miss=handle.miss,eps=eps,tol=tol,
        maxit=maxit,pl=pl,assignment=assignment,verbose=verbose)
}

#' Gene counting for haplotype analysis
#'
#' Gene counting for haplotype analysis with missing data
#'
#' @param data genotype table.
#' @param weight a column of frequency weights.
#' @param loci an array containing number of alleles at each locus.
#' @param control is a function with the following arguments: 
#'      \enumerate{
#'      \item xdata. a flag indicating if the data involves X chromosome, if so, the
#'first column of data indicates sex of each subject: 1=male, 2=female. The marker 
#'data are no different from the autosomal version for females, but for males, two copies
#'of the single allele present at a given locus.
#'      \item convll. set convergence criteria according to log-likelihood, if its value set to 1
#'      \item handle.miss. to handle missing data, if its value set to 1
#'      \item eps. the actual convergence criteria, with default value 1e-5
#'      \item tol. tolerance for genotype probabilities with default value 1e-8
#'      \item maxit. maximum number of iterations, with default value 50
#'      \item pl. criteria for trimming haplotypes according to posterior probabilities
#'      \item assignment. filename containing haplotype assignment
#'      \item verbose. If TRUE, yields print out from the C routine
#'      }
#'
#' @export
#' @return The returned value is a list containing:
#' \describe{
#' \item{h}{haplotype frequency estimates under linkage disequilibrium (LD)}
#' \item{h0}{haplotype frequency estimates under linkage equilibrium (no LD)}
#' \item{prob}{genotype probability estimates}
#' \item{l0}{log-likelihood under linkage equilibrium}
#' \item{l1}{log-likelihood under linkage disequilibrium}
#' \item{hapid}{unique haplotype identifier (defunct, see gc.em)}
#' \item{npusr}{number of parameters according user-given alleles}
#' \item{npdat}{number of parameters according to observed}
#' \item{htrtable}{design matrix for haplotype trend regression (defunct, see gc.em)}
#' \item{iter}{number of iterations used in gene counting}
#' \item{converge}{a flag indicating convergence status of gene counting}
#' \item{di0}{haplotype diversity under no LD, defined as \eqn{1-\sum (h_0^2)}{1-sum (h0^2)}}
#' \item{di1}{haplotype diversity under LD, defined as \eqn{1-\sum (h^2))}{1-sum (h^2)}}
#' \item{resid}{residuals in terms of frequency weights = o - e}
#' }
#'
#' @references
#'
#' Zhao, J. H., Lissarrague, S., Essioux, L. and P. C. Sham (2002).
#' GENECOUNTING: haplotype analysis with missing genotypes.
#' Bioinformatics 18(12):1694-1695
#'
#' Zhao, J. H. and P. C. Sham (2003). Generic number systems and haplotype
#' analysis. Comp Meth Prog Biomed 70: 1-9
#'
#' Zhao, J. H. (2004). 2LD, GENECOUNTING and HAP: Computer programs for linkage
#' disequilibrium analysis. Bioinformatics, 20, 1325-1326 
#'
#' @seealso \code{\link[gap]{gc.em}}, \code{\link[gap]{LDkl}}
#'
#' @examples
#' \dontrun{
#' require(gap.datasets)
#' # HLA data
#' data(hla)
#' hla.gc <- genecounting(hla[,3:8])
#' summary(hla.gc)
#' hla.gc$l0
#' hla.gc$l1
#'
#' # ALDH2 data
#' data(aldh2)
#' control <- gc.control(handle.miss=1,assignment="ALDH2.out")
#' aldh2.gc <- genecounting(aldh2[,3:6],control=control)
#' summary(aldh2.gc)
#' aldh2.gc$l0
#' aldh2.gc$l1
#'
#' # Chromosome X data
#' # assuming allelic data have been extracted in columns 3-13
#' # and column 3 is sex
#' filespec <- system.file("tests/genecounting/mao.dat")
#' mao2 <- read.table(filespec)
#' dat <- mao2[,3:13]
#' loci <- c(12,9,6,5,3)
#' contr <- gc.control(xdata=TRUE,handle.miss=1)
#' mao.gc <- genecounting(dat,loci=loci,control=contr)
#' mao.gc$npusr
#' mao.gc$npdat
#' }
#' @author Jing Hua Zhao
#' @note adapted from GENECOUNTING.
#' @keywords models

genecounting <- function(data,weight=NULL,loci=NULL,control=gc.control())
{
  if (control$xdata)
  {
     sex <- data[,1]
     data <- data[,-1]
  }
  else
  sex <- rep(2,dim(data)[1])
  if(is.null(weight)) weight<-rep(1,dim(data)[1])
# precis<-1
# to call dpmach
# tol<-1.2
#  while(tol>1.0)
#  {
#     precis<-precis/2.0;
#     tol<-1.0+precis;
#  }
  precis<-.Machine$double.eps
  gid<-1:(dim(data)[1])
  nloci=dim(data)[2]/2
  if(is.null(loci))
  {
    loci<-rep(0,nloci)
    for (i in 1:nloci)
    {
        loci[i]=max(data[,c(2*i-1,2*i)],na.rm=TRUE)
    }
  }
  data<-as.matrix(data)
  data<-t(data)
  hapall<-1
  for(i in 1:nloci)
  {
    hapall<-hapall*loci[i];
  }
  h0<-h1<-hapid<-rep(0,hapall)
  obscom<-length(weight)
  prob<-rep(0,obscom)
  lnl0<-lnl1<-0
  npusr<-npdat<-rep(0,2)
# 13/11/2003
# change to reduce memory request
# htrtable<-matrix(rep(0,obscom*hapall),nrow=obscom)
  verbose<-0
  iter<-0
  converge<-0
  z <- .C("gcx",verbose=as.integer(control$verbose),
           Rhandlemissing=as.integer(control$handle.miss),
           Rconvll=as.integer(control$convll),
           Reps=as.double(control$eps),
           Rtol=as.double(control$tol),
           Rmaxit=as.integer(control$maxit),
           Rpl=as.double(control$pl),
           precis=as.double(precis),
           gid=as.integer(gid),
           Rnloci=as.integer(nloci),
           Rloci=as.integer(loci),
           Robscom=as.integer(obscom),
           Rhapall=as.integer(hapall),
           genotype=as.integer(data),
           count=as.integer(weight),
           Rxdata=as.integer(control$xdata),
           sex=as.integer(sex),
           hapid=as.integer(hapid),
           prob=as.double(prob),
           Rh0=as.double(h0),
           Rh1=as.double(h1),
           lnl0=as.double(lnl0),
           lnl1=as.double(lnl1),
           npusr=as.integer(npusr),
           npdat=as.integer(npdat),
           iter=as.integer(iter),
           converge=as.integer(converge),assignment=as.character(control$assignment),PACKAGE="gap"
           )
  x<-0
  hapid<-0
# Dprime<-sum(z$Rh0*abs(z$Rh1-z$Rh0))
# x<-t(matrix(z$htrtable/2,nrow=hapall))
# hapid<-apply(x,2,sum)>0
# x<-x[,hapid]
# hapid<-(1:hapall)[hapid]
  di0<-1-sum((z$Rh0)^2)
  di1<-1-sum((z$Rh1)^2)
  prob<-z$prob/sum(z$prob)
  resid<-weight*(1-prob)
  list(h=z$Rh1, h0=z$Rh0, prob=prob, l0=z$lnl0, l1=z$lnl1,
       hapid=hapid, npusr=z$npusr, npdat=z$npdat, htrtable=x,
       iter=z$iter,converge=z$converge,di0=di0,di1=di1,resid=resid)
}
