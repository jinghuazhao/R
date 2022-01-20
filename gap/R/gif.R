#' Kinship coefficient and genetic index of familiality
#'
#' The genetic index of familality is defined as the mean kinship between
#' all pairs of individuals in a set multiplied by 100,000. Formally, it 
#' is defined as 
#' \deqn{100,000 \times \frac{2}{n(n-1)}\sum_{i=1}^{n-1}\sum_{j=i+1}^n k_{ij}}{100,000 x 2/[n(n-1)]\sum_(i=1)^(n-1)\sum_(j=i+1)^n k_(ij)}
#' where \eqn{n} is the number of individuals in the set and \eqn{k_{ij}} is the
#' kinship coefficient between individuals \eqn{i} and \eqn{j}.
#'
#' The scaling is purely for convenience of presentation.
#'
#' @param data the trio data of a pedigree.
#' @param gifset a subgroup of pedigree members.
#'
#' @export
#' @return The returned value is a list containing:
#' \describe{
#' \item{gifval}{the genetic index of familiarity.}
#' }
#'
#' @references
#' Gholami K, Thomas A (1994) A linear time algorithm for calculation of
#' multiple pairwise kinship coefficients and genetic index of familiality.
#' Comp Biomed Res 27:342-350
#'
#' @seealso \code{\link[gap]{pfc}}
#'
#' @examples
#' \dontrun{
#' test<-c(
#'  5,      0,      0,
#'  1,      0,      0,
#'  9,      5,      1,
#'  6,      0,      0,
#' 10,      9,      6,
#' 15,      9,      6,
#' 21,     10,     15,
#'  3,      0,      0,
#' 18,      3,     15,
#' 23,     21,     18,
#'  2,      0,      0,
#'  4,      0,      0,
#'  7,      0,      0,
#'  8,      4,      7,
#' 11,      5,      8,
#' 12,      9,      6,
#' 13,      9,      6,
#' 14,      5,      8,
#' 16,     14,      6,
#' 17,     10,      2,
#' 19,      9,     11,
#' 20,     10,     13,
#' 22,     21,     20)
#' test<-matrix(test,ncol=3,byrow=TRUE)
#' gif(test,gifset=c(20,21,22))
#'
#' # all individuals
#' gif(test,gifset=1:23)
#' }
#'
#' @author Alun Thomas, Jing Hua Zhao
#' @note Adapted from gif.c, testable with -Dexecutable as standalone program, 
#' which can be use for any pair of indidivuals
#' @keywords datagen

gif <- function(data,gifset)
{
  famsize<-dim(data)[1]
  giflen<-length(gifset)
  gifval<-0
  z<-.C("gif_c",data=as.integer(t(data)),famsize=as.integer(famsize),
        gifset=as.integer(array(gifset)),giflen=as.integer(giflen),
        gifval=as.double(gifval),PACKAGE="gap")

  list(gifval=z$gifval)
}
