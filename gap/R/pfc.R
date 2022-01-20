#' Probability of familial clustering of disease
#'
#' To calculate exact probability of familial clustering of disease
#'
#' @param famdata collective information of sib size, number of affected sibs and their frequencies.
#' @param enum a switch taking value 1 if all possible tables are to be enumerated.
#'
#' @export
#' @return The returned value is a list containing (tailp,sump,nenum are only available if enum=1:
#' \describe{
#'   \item{p}{the probabitly of familial clustering}
#'   \item{stat}{the deviances, chi-squares based on binomial and hypergeometric distributions, 
#'  the degrees of freedom should take into account the number of marginals used}
#'   \item{tailp}{the exact statistical significance}
#'   \item{sump}{sum of the probabilities used for error checking}
#'   \item{nenum}{the total number of tables enumerated}
#' }
#'
#' @references
#' Yu C, Zelterman D (2001) Exact inference for family disease clusters. Commun Stat -- Theory
#' Meth 30:2293-2305
#'
#' Yu C, Zelterman D (2002) Statistical inference for familial disease clusters. Biometrics
#' 58:481-491
#'
#' @seealso \code{\link[gap]{kin.morgan}}
#'
#' @examples
#' \dontrun{
#' # IPF among 203 siblings of 100 COPD patients from Liang KY, SL Zeger,
#' # Qaquish B. Multivariate regression analyses for categorical data
#' # (with discussion). J Roy Stat Soc B 1992, 54:3-40
#'
#' # the degrees of freedom is 15
#' famtest<-c(
#' 1, 0, 36,
#' 1, 1, 12,
#' 2, 0, 15,
#' 2, 1,  7,
#' 2, 2,  1,
#' 3, 0,  5,
#' 3, 1,  7,
#' 3, 2,  3,
#' 3, 3,  2,
#' 4, 0,  3,
#' 4, 1,  3,
#' 4, 2,  1,
#' 6, 0,  1,
#' 6, 2,  1,
#' 6, 3,  1,
#' 6, 4,  1,
#' 6, 6,  1)
#' test<-t(matrix(famtest,nrow=3))
#' famp<-pfc(test)
#' }
#'
#' @author Dani Zelterman, Jing Hua Zhao
#' @note Adapted from family.for by Dani Zelterman, 25/7/03
#' @keywords models

pfc <- function(famdata,enum=0)
{
  famsize<-dim(famdata)[1]
  pobs<-p<-tailp<-sump<-1.0
  stat<-rep(0,20)
  nenum<-0
  cat("Family frequency data read in:\n")
  cat("Sibs\tAffected\tFrequency\n")
  for(i in 1:famsize)
  {
     cat(famdata[i,1], "\t", famdata[i,2], "\t", famdata[i,3], "\n")
  }
  z<-.Fortran("family_",famdata=as.integer(matrix(famdata,ncol=3)),famsize=as.integer(famsize),
               pobs=as.double(pobs),p=as.double(p),stat=as.double(stat),toenum=as.integer(enum),
               tailp=as.double(tailp),sump=as.double(sump),nenum=as.double(nenum),PACKAGE="gap")
  cat("Probability of this table: ",z$pobs,"\n")
  if(enum==0) list(p=z$p,stat=z$stat[1:5])
  else list(p=z$p,stat=z$stat[1:5],tailp=z$tailp,sump=z$sump,nenum=z$nenum)
}
