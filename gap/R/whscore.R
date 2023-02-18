#' Whittemore-Halpern scores for allele-sharing
#'
#' @param allele a matrix of alleles of affected pedigree members.
#' @param type 0 = pairs, 1 = all.
#'
#' @details
#' Allele sharing score statistics.
#'
#' @export
#' @return
#' The returned value is the value of score statistic.
#'
#' @references
#' \insertRef{kruglyak96}{gap}
#'
#' \insertRef{wh94a}{gap}
#'
#' \insertRef{wh94b}{gap}
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
