#' Control for haplotype reconstruction
#'
#' @param mb Maximum dynamic storage to be allocated, in Mb.
#' @param pr Prior (ie population) probability threshold.
#' @param po Posterior probability threshold.
#' @param to Log-likelihood convergence tolerance.
#' @param th Posterior probability threshold for output.
#' @param maxit Maximum EM iteration.
#' @param n Force numeric allele coding (1/2) on output (off).
#' @param ss Tab-delimited speadsheet file output (off).
#' @param rs Random starting points for each EM iteration (off).
#' @param rp Restart from random prior probabilities.
#' @param ro Loci added in random order (off).
#' @param rv Loci added in reverse order (off).
#' @param sd Set seed for random number generator (use date+time).
#' @param mm Repeat final maximization multiple times.
#' @param mi Create multiple imputed datasets. If set >0.
#' @param mc Number of MCMC steps between samples.
#' @param ds Starting value of Dirichlet prior parameter.
#' @param de Finishing value of Dirichlet prior parameter.
#' @param q Quiet operation (off).
#' @param hapfile a file for haplotype frequencies.
#' @param assignfile a file for haplotype assignment.
#'
#' @export
#' @return
#' A list containing the parameter specifications to the function.

hap.control <- function(mb=0,pr=0,po=0.001,to=0.001,th=1,maxit=100,n=0,
      ss=0,rs=0,rp=0,ro=0,rv=0,sd=0,mm=0,mi=0,mc=50,ds=0.1,de=0,q=0,
      hapfile="hap.out",assignfile="assign.out")
{
   list(mb=mb,pr=pr,po=po,to=to,th=th,maxit=maxit,n=n,ss=ss,rs=rs,
       rp=rp,ro=ro,rv=rv,sd=sd,mm=mm,mi=mi,mc=mc,ds=ds,de=de,q=q,
       hapfile=hapfile,assignfile=assignfile)
}

