hap.control <- function(mb=0,pr=0,po=0.001,to=0.001,th=1,maxit=100,n=0,
      ss=0,rs=0,rp=0,ro=0,rv=0,sd=0,mm=0,mi=0,mc=50,ds=0.1,de=0,q=0,
      hapfile="hap.out",assignfile="assign.out")
{
   list(mb=mb,pr=pr,po=po,to=to,th=th,maxit=maxit,n=n,ss=ss,rs=rs,
       rp=rp,ro=ro,rv=rv,sd=sd,mm=mm,mi=mi,mc=mc,ds=ds,de=de,q=q,
       hapfile=hapfile,assignfile=assignfile)
}

#' Haplotype reconstruction
#'
#' Haplotype reconstruction using sorting and trimming algorithms.
#'
#' @param id a column of subject id.
#' @param data genotype table.
#' @param nloci number of loci.
#' @param loci number of alleles at all loci.
#' @param names locus names.
#' @param control is a function with the following arguments,
#'     \enumerate{
#'     \item mb Maximum dynamic storage to be allocated, in Mb
#'     \item pr Prior (ie population) probability threshold
#'     \item po Posterior probability threshold
#'     \item to Log-likelihood convergence tolerance
#'     \item th Posterior probability threshold for output
#'     \item maxit Maximum EM iteration
#'     \item n Force numeric allele coding (1/2) on output (off)
#'     \item ss Tab-delimited speadsheet file output (off)
#'     \item rs Random starting points for each EM iteration (off)
#'     \item rp Restart from random prior probabilities
#'     \item ro Loci added in random order (off)
#'     \item rv Loci added in reverse order (off)
#'     \item sd Set seed for random number generator (use date+time)
#'     \item mm Repeat final maximization multiple times
#'     \item mi Create multiple imputed datasets. If set >0
#'     \item mc  Number of MCMC steps between samples
#'     \item ds  Starting value of Dirichlet prior parameter
#'     \item de  Finishing value of Dirichlet prior parameter
#'     \item q Quiet operation (off)
#'     \item hapfile a file for haplotype frequencies
#'     \item assignfile a file for haplotype assignment
#'     }
#'
#' @details
#' The package can hanlde much larger number of multiallelic loci. 
#' For large sample size with relatively small number of multiallelic
#' loci, genecounting should be used.
#'
#' @export
#' @return The returned value is a list containing:
#' \describe{
#' \item{l1}{log-likelihood assuming linkage disequilibrium}
#' \item{converge}{convergence status, 0=failed, 1=succeeded}
#' \item{niter}{number of iterations}
#' }
#'
#' @references
#'
#' Clayton DG (2001) SNPHAP. https://github.com/chr1swallace/snphap.
#'
#' Zhao JH and W Qian (2003) Association analysis of unrelated individuals
#' using polymorphic genetic markers. RSS 2003, Hassalt, Belgium
#'
#' Zhao JH (2004). 2LD, GENECOUNTING and HAP: Computer programs for linkage
#' disequilibrium analysis. Bioinformatics 20: 1325-1326
#'
#' @seealso \code{\link[gap]{genecounting}}
#'
#' @examples
#' \dontrun{
#' require(gap.datasets)
#' # 4 SNP example, to generate hap.out and assign.out alone
#' data(fsnps)
#' hap(id=fsnps[,1],data=fsnps[,3:10],nloci=4)
#' dir()
#'
#' # to generate results of imputations
#' control <- hap.control(ss=1,mi=5,hapfile="h",assignfile="a")
#' hap(id=fsnps[,1],data=fsnps[,3:10],nloci=4,control=control)
#' dir()
#' }
#'
#' @note adapted from hap.
#' @keywords models

hap<-function(id,data,nloci,loci=rep(2,nloci),names=paste("loci",1:nloci,sep=""),
              control=hap.control())
{
  if (control$rv & control$ro) stop("rv and ro flags cannot both be set in hap.control\n");
# if (mi==0 & (mc | ds | de)) stop("mc, ds, de parameters are only legal if mi is set\n");
  if (control$rp & control$mm==0) stop("rp option only relevant with mm # option\n");
  nobs<-dim(data)[1]
  data<-as.matrix(data)
  if(length(id)!=dim(data)[1]) stop("id and data should have the same length")
  l1<-niter<-converge<-0

  z<-.C("hap_c",nobs=as.integer(nobs),idstr=as.character(id),data=as.character(t(data)),
        nloci=as.integer(nloci),loci=as.integer(loci),names=as.character(names),
        mb=as.double(control$mb),
        pr=as.double(control$pr),
        po=as.double(control$po),
        to=as.double(control$to),
        th=as.double(control$th),
        maxitt=as.double(control$maxit),
        n=as.integer(control$n),
        sst=as.integer(control$ss),
        rst=as.integer(control$rs),
        rp=as.integer(control$rp),
        ro=as.integer(control$ro),
        rv=as.integer(control$rv),
        sd=as.double(control$sd),
        mm=as.integer(control$mm),
        mi=as.integer(control$mi),
        mc=as.integer(control$mc),
        ds=as.double(control$ds),
        de=as.double(control$de),
        q=as.integer(control$q),
        l1=as.double(l1),
        niter=as.integer(niter),
        converged=as.integer(converge),
        hapfile=as.character(control$hapfile),
        assignfile=as.character(control$assignfile),
        PACKAGE="gap")

  list(l1=z$l1,converge=z$converged,niter=z$niter)
}
