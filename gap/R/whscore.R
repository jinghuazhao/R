#' Whittemore-Halpern scores for allele-sharing
#'
#' Allele sharing score statistics.
#'
#' @param allele a matrix of alleles of affected pedigree members.
#' @param type 0 = pairs, 1 = all.
#'
#' @export
#' @return The returned value is the value of score statistic.
#'
#' @references
#' Kruglyak L, Daly MJ, Reeve-Daly MP, Lander ES (1996) Parametric and Nonparametric
#' linkage analysis: a unified multipoint approach. Am. J. Hum. Genet. 58:1347-1363
#'
#' Whittemore AS, Halpern J (1994) A class of tests for linkage using affected 
#' pedigree members. Biometrics 50:118-127
#'
#' Whittemore AS, Halpern J (1994) Probability of gene identity by descent: 
#' computation and applications. Biometrics 50:109-117
#'
#'
#' @examples
#' \dontrun{
#' c<-matrix(c(1,1,1,2,2,2),ncol=2)
#' whscore(c,type=1)
#' whscore(c,type=2)
#' }
#'
#' @author Leonid Kruglyak, Jing Hua Zhao
#' @note adapted from GENEHUNTER.
#' @keywords utilities

whscore <- function(allele,type)
{
  n<-dim(allele)[1]
  s<-0
  if(type==1)
    z<-.C("score_pairs",data=as.integer(t(allele)),n=as.integer(n),arscore=as.double(s),PACKAGE="gap")
  else
    z<-.C("score_all",data=as.integer(t(allele)),n=as.integer(n),arscore=as.double(s),PACKAGE="gap")

  z$arscore
}
