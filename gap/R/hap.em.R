#' Gene counting for haplotype analysis
#'
#' Gene counting for haplotype analysis with missing data, adapted for hap.score.
#'
#' @param id a vector of individual IDs.
#' @param data Matrix of alleles, such that each locus has a pair of
#' adjacent columns of alleles, and the order of columns
#' corresponds to the order of loci on a  chromosome. If
#' there are K loci, then ncol(data) = 2*K. Rows represent
#' alleles for each subject.
#' @param locus.label Vector of  labels  for  loci,  of  length  K  (see definition of data matrix).
#' @param converge.eps Convergence criterion, based on absolute  change in log likelihood (lnlike).
#' @param maxiter Maximum number of iterations of EM.
#' @param miss.val missing value.
#'
#' @export
#' @return
#' List with components:
#' \describe{
#'  \item{converge}{Indicator of convergence of the EM algorithm
#'  (1=converged, 0 = failed).}
#'  \item{niter}{Number of iterations completed in the EM alogrithm.}
#'  \item{locus.info}{A list with  a  component for each locus.  Each
#'   component is also a list, and  the  items of a locus-
#'   specific list are the locus name and a vector for the
#'   unique alleles for the locus.}
#'  \item{locus.label}{Vector of  labels  for  loci,  of  length  K  (see
#'    definition of input values).}
#'  \item{haplotype}{Matrix of unique haplotypes. Each row represents a
#'   unique  haplotype, and the number of columns is the number of loci.}
#'  \item{hap.prob}{Vector of mle's of haplotype probabilities.  The ith
#'   element of hap.prob corresponds to the ith row of haplotype.}
#'  \item{lnlike}{Value of lnlike at last EM iteration (maximum lnlike if converged).}
#'  \item{indx.subj}{Vector for index of subjects, after  expanding  to
#'   all possible  pairs  of  haplotypes  for  each person. If
#'   indx=i, then i is the ith row of input matrix data. If the
#'   ith subject has  n possible  pairs  of haplotypes that
#'   correspond to their marker phenotype, then i is repeated n times.}
#'  \item{nreps}{Vector for the count of haplotype pairs that map to
#'   each subject's marker genotypes.}
#'  \item{hap1code}{Vector of codes for each subject's first haplotype.
#'   The values in hap1code are the row numbers of the unique
#'   haplotypes in the returned matrix haplotype.}
#'  \item{hap2code}{Similar to hap1code, but for  each  subject's  second haplotype.}
#'  \item{post}{Vector of posterior probabilities of pairs of
#'   haplotypes for a person, given thier marker phenotypes.}
#' }
#'
#' @seealso \code{\link[gap]{hap}}, \code{\link[gap]{LDkl}}
#'
#' @examples
#' \dontrun{
#' data(hla)
#' hap.em(id=1:length(hla[,1]),data=hla[,3:8],locus.label=c("DQR","DQA","DQB"))
#' }
#'
#' @author Jing Hua Zhao
#' @note Adapted from HAP.
#' @keywords models

hap.em <- function(id,data,locus.label=NA,converge.eps=0.000001,maxiter=500,miss.val=0)
{
  for(p in c("haplo.stats")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
        warning(paste("hap.em needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  tmp0 <- geno.recode(data,miss.val=miss.val)
  geno <- as.matrix(tmp0$grec)
  geno[is.na(geno)] <- 0
  data<-as.matrix(geno)
  nloci<-dim(data)[2]/2
  loci<-rep(0,nloci)
  for (i in 1:nloci)
  {
     loci[i]<-length(tmp0$alist[[i]]$allele) #  max(data[,c(2*i-1,2*i)],na.rm=TRUE)
  }
  if(all(is.na(locus.label))) {
     locus.label<- paste("loc-",1:nloci,sep="")
  }
##
  zz <- hap.control(ss=1)
  z<- hap(id=id,data=data,nloci,loci=loci,control=zz)
#
  tmp1<-read.table(zz$hapfile,header=T)
# unlink("hap.out")
  haplotype<-as.matrix(tmp1[,1:nloci])
  dimnames(haplotype)<-list(1:length(haplotype[,1]),locus.label)
  uhap<-hapid<-1:(dim(tmp1)[1])
  tmp1<-data.frame(tmp1,hapid)
  hap.prob<-tmp1[,nloci+1]
#
  tmp2<-read.table(zz$assignfile,header=T)
# unlink("assign.out")
  nrow<-dim(tmp2)[1]/2
  indx1<-2*1:nrow-1
  indx2<-2*1:nrow
  indx.subj<-tmp2[indx1,1]
  post<-tmp2[indx1,nloci+3]

  tmp<-merge(tmp1[,-(nloci+1)],tmp2[,-c(1,2)],sort=F)
  hap1<-tmp[indx1,(nloci+1)]
  hap2<-tmp[indx2,(nloci+1)]
  nreps<-tapply(indx.subj,indx.subj,length)

  list (lnlike=z$l1,hap.prob=hap.prob,indx.subj=indx.subj,post=post,
        hap1code=hap1,hap2code=hap2,haplotype=grec2g(haplotype,nloci,tmp0),nreps=nreps,
        converge=z$converge,niter=z$niter,uhap=uhap)
}

# 15-10-03 WAH
# 16-10-03 UCL office to add sort=F in merge
# 25-09-04 first attempt to robust label-handling
# 29-07-07 add library(haplo.stats)
