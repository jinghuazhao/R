#' Means and variances under 1- and 2- locus (biallelic) QTL model
#'
#' Function muvar() gives means and variances under 1-locus and 2-locus QTL model (simple); 
#' in the latter case it gives results from different avenues. This function is included for
#' experimental purpose and yet to be generalized.
#'
#' @param n.loci number of loci, 1=single locus, 2=two loci.
#' @param y1 the genotypic means of aa, Aa and AA.
#' @param p1 the frequency of the lower allele, or the that for the first locus under a 2-locus model.
#' @param y12 the genotypic means of aa, Aa and AA at the first locus and bb, Bb and BB at the second locus.
#' @param p2 the frequency of the lower allele at the second locus.
#'
#' @export
#' @return Currently it does not return any value except screen output; the results can be kept via R's sink()
#' command or via modifying the C/R codes.
#'
#' @references
#' Sham P (1998). Statistics in Human Genetics. Arnold
#'
#' @examples
#' \dontrun{
#' # the default 1-locus model
#' muvar(n.loci=1,y1=c(0,1,1),p1=0.5)
#'
#' # the default 2-locus model
#' muvar(n.loci=2,y12=c(1,1,1,1,1,0,0,0,0),p1=0.99,p2=0.9)
#' }
#'
#' @author Jing Hua Zhao
#' @note Adapted from an earlier C program written for the above book.
#' @keywords models utilties

muvar <- function(n.loci=1,y1=c(0,1,1),y12=c(1,1,1,1,1,0,0,0,0),p1=0.99,p2=0.9)
{
  if(n.loci==1)
    .C("onelocus",y1=as.single(y1),p1=as.single(p1),PACKAGE="gap")
  else
    .C("twolocus",y12=as.single(y12),p1=as.single(p1),p2=as.single(p2),PACKAGE="gap")
}
