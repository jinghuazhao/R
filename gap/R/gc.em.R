#' Gene counting for haplotype analysis
#'
#' @param locus.label Vector of  labels  for  loci,  of  length  K  (see definition of data matrix).
#' @param converge.eps Convergence criterion, based on absolute  change in log likelihood (lnlike).
#' @param maxiter Maximum number of iterations of EM.
#' @param handle.miss a flag for handling missing genotype data, 0=no, 1=yes.
#' @param miss.val missing value.
#' @param control a function, see [genecounting].
#'
#' @details
#' Gene counting for haplotype analysis with missing data, adapted for hap.score
#'
#' @param data Matrix of alleles, such that each locus has a pair of
#' adjacent columns of alleles, and the order of columns
#' corresponds to the order of loci on a chromosome. If
#' there are K loci, then ncol(data) = 2*K. Rows represent
#' alleles for each subject.
#'
#' @export
#' @return
#' List with components:
#' - converge Indicator of convergence of the EM algorithm (1=converged, 0 = failed).
#' - niter Number of iterations completed in the EM alogrithm.
#' - locus.info A list with  a  component for each locus.  Each
#'   component is also a list, and  the  items of a locus-
#'   specific list are the locus name and a vector for the
#'   unique alleles for the locus.
#' - locus.label Vector of labels for loci, of length K (see definition of input values).
#' - haplotype Matrix of unique haplotypes. Each row represents a
#'   unique haplotype, and the number of columns is the number of loci.
#' - hap.prob Vector of mle's of haplotype probabilities.  The ith
#'   element of hap.prob corresponds to the ith row of haplotype.
#' - hap.prob.noLD Similar to hap.prob, but assuming no linkage disequilibrium.
#' - lnlike Value of lnlike at last EM iteration (maximum lnlike if converged).
#' - lr Likelihood ratio statistic to test no linkage disequilibrium among all loci.
#' - indx.subj Vector for index of subjects, after  expanding  to
#'   all possible  pairs  of  haplotypes  for  each person. If
#'   indx=i, then i is the ith row of input matrix data. If the
#'   ith subject has  n possible  pairs  of haplotypes that
#'   correspond to their marker phenotype, then i is repeated n times.
#' - nreps Vector for the count of haplotype pairs that map to
#'   each subject's marker genotypes.
#' - hap1code Vector of codes for each subject's first haplotype.
#'   The values in hap1code are the row numbers of the unique
#'   haplotypes in the returned matrix haplotype.
#' - hap2code Similar to hap1code, but for  each  subject's  second haplotype.
#' - post Vector of posterior probabilities of pairs of
#'   haplotypes for a person, given thier marker phenotypes.
#' - htrtable A table which can be used in haplotype trend regression.
#'
#' @references
#' \insertRef{zhao02}{gap}
#'
#' \insertRef{zhao03}{gap}
#'
#' @seealso [`genecounting`], [`LDkl`]
#'
#' @examples
#' \dontrun{
#' data(hla)
#' gc.em(hla[,3:8],locus.label=c("DQR","DQA","DQB"),control=gc.control(assignment="t"))
#' }
#'
#' @author Jing Hua Zhao
#' @note Adapted from GENECOUNTING.
#' @keywords models

gc.em <- function(data, locus.label=NA, converge.eps=0.000001, maxiter=500, 
         handle.miss=0, miss.val=0, control=gc.control())
{
  if (control$xdata) {
     sex <- data[,1]
     data <- data[,-1]
  }
  for(p in c("haplo.stats")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
        warning(paste("gc.em needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  tmp0<-geno.recode(data,miss.val=miss.val)
  geno<-tmp0$grec
  geno[is.na(geno)]<-0
  data<-as.matrix(geno)
  weight<-rep(1,dim(data)[1])
  nloci<-dim(data)[2]/2
  loci<-rep(0,nloci)
  for (i in 1:nloci)
  {
      loci[i]<-length(tmp0$alist[[i]]$allele) # max(data[,c(2*i-1,2*i)],na.rm=TRUE)
  }
  if(all(is.na(locus.label))) {
     locus.label<- paste("loc-",1:nloci,sep="")
  }
# to run genecounting
  if (control$xdata) data <- cbind(sex,data)
  data.gc<-genecounting(data,weight=weight,loci=loci,
           control=gc.control(xdata=control$xdata,
           eps=converge.eps,pl=0.001,maxit=maxiter,
           handle.miss=handle.miss,assignment=control$assignment,verbose=F))
  hap.prob<-data.gc$h
  hap.prob.noLD<-data.gc$h0
  lnlike<-data.gc$l1
  lr<-2*(data.gc$l1-data.gc$l0)
  df<-data.gc$npdat-sum(loci)-length(loci)
  niter<-data.gc$iter
  converge<-data.gc$converge
# to further extract information and obtain unique haplotypes
  hapas<-read.table(control$assignment)
  unlink(control$assignment)
  newnames<-c("subj","chr",locus.label,"post","hapid")
  names(hapas)<-newnames
  ncol<-nloci+4
  nrow<-dim(hapas)[1]/2
  indx1<-2*1:nrow-1
  indx2<-2*1:nrow
  indx.subj<-hapas$subj[indx1]
  hapdat<-hapas[,-c(1,2,ncol-1)]
  post<-hapas$post[indx1]
  hapid<-hapas$hapid
  one<-rep(1,nrow*2)
  hapdat<-cbind(hapdat,one)
  with(hapdat,{
  tmp<-by(hapdat,one,unique)
  haplotype<-as.matrix(tmp[[1]])
  tmp<-order(haplotype[,nloci+1])
  haplotype<-haplotype[tmp,1:(dim(haplotype)[2]-2)]
  dimnames(haplotype)<-list(1:length(haplotype[,1]),locus.label)
  hap1<-hapid[indx1]
  hap2<-hapid[indx2]
  uhap <- sort(unique(hapid))
  hap.prob<-hap.prob[uhap]
  nreps<-tapply(indx.subj,indx.subj,length)
# 13/11/2003
# haplotype trend regression
# assign.dat already has sequential number to avoid duplicate IDs
  idx.subj<-sort(unique(indx.subj))
  N<-length(idx.subj)
  P<-length(uhap)
  idx.subj<-cbind(1:N,idx.subj)
  idx.uhap<-cbind(1:P,uhap)
  htrtable<-matrix(rep(0,N*P),nrow=N)
  for(l in 1:nrow)
  {
    i<-idx.subj[,1][idx.subj[,2]==indx.subj[l]]
    j1<-idx.uhap[,1][idx.uhap[,2]==hap1[l]]
    htrtable[i,j1]<-htrtable[i,j1]+post[l]
    j2<-idx.uhap[,1][idx.uhap[,2]==hap2[l]]
    htrtable[i,j2]<-htrtable[i,j2]+post[l]
  }
  htrtable<-htrtable/2
  dimnames(htrtable)<-list(NULL,as.character(uhap))
  list(lnlike=lnlike,lr=lr,
       hap.prob=hap.prob,hap.prob.noLD=hap.prob.noLD,indx.subj=indx.subj,
       post=post,hap1code=hap1,hap2code=hap2,haplotype=grec2g(haplotype,nloci,tmp0),
       nreps=nreps,converge=converge,niter=niter,uhap=uhap,htrtable=htrtable)
  })
}
