#' Probability of familial clustering of disease
#'
#' @param famdata collective information of sib size, number of affected sibs and their frequencies.
#' @param n.sim number of simulations in a single Monte Carlo run.
#' @param n.loop total number of Monte Carlo runs.
#'
#' @details
#' To calculate probability of familial clustering of disease using Monte Carlo simulation.
#'
#' @export
#' @return
#' The returned value is a list containing:
#' - n.sim a copy of the number of simulations in a single Monte Carlo run.
#' - n.loop the total number of Monte Carlo runs.
#' - p the observed p value.
#' - tailpl accumulated probabilities at the lower tails.
#' - tailpu simulated p values.
#'
#' @references
#' \insertRef{yu01}{gap}
#'
#' @seealso [`pfc`]
#'
#' @examples
#' \dontrun{
#' # Li FP, Fraumeni JF Jr, Mulvihill JJ, Blattner WA, Dreyfus MG, Tucker MA,
#' # Miller RW. A cancer family syndrome in twenty-four kindreds.
#' # Cancer Res 1988, 48(18):5358-62. 
#'
#' # family_size  #_of_affected frequency
#'
#' famtest<-c(
#' 1, 0, 2,
#' 1, 1, 0,
#' 2, 0, 1,
#' 2, 1, 4,
#' 2, 2, 3,
#' 3, 0, 0,
#' 3, 1, 2,
#' 3, 2, 1,
#' 3, 3, 1,
#' 4, 0, 0,
#' 4, 1, 2,
#' 5, 0, 0,
#' 5, 1, 1,
#' 6, 0, 0,
#' 6, 1, 1,
#' 7, 0, 0,
#' 7, 1, 1,
#' 8, 0, 0,
#' 8, 1, 1,
#' 8, 2, 1,
#' 8, 3, 1,
#' 9, 3, 1)
#'
#' test<-matrix(famtest,byrow=T,ncol=3)
#'
#' famp<-pfc.sim(test)
#' }
#'
#' @author Chang Yu, Dani Zelterman
#' @note Adapted from runi.for from Change Yu, 5/6/4
#' @keywords models

pfc.sim <- function(famdata,n.sim=1000000,n.loop=1)
{
  famsize<-dim(famdata)[1]
  obsp<-0
  tailpl<-tailpu<-rep(0,n.loop)
  z<-.Fortran("runifamily",famdata=as.integer(matrix(famdata,ncol=3)),famsize=as.integer(famsize),
               nsim=as.integer(n.sim),ncycle=as.integer(n.loop),
               obsp=as.double(obsp),tailpl=as.double(tailpl),tailpu=as.integer(tailpu),PACKAGE="gap")

  list(n.sim=n.sim,n.loop=n.loop,p=z$obsp,tailpl=z$tailpl,tailpu=z$tailpu/n.sim)
}
