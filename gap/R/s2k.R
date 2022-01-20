#' Statistics for 2 by K table
#'
#' This function calculates one-to-others and maximum accumulated chi-squared
#' statistics for a 2 by K contingency table.
#'
#' @param y1 a vector containing the first row of a 2 by K contingency table.
#' @param y2 a vector containing the second row of a 2 by K contingency table.
#'
#' @export
#' @return
#' The returned value is a list containing:
#' \describe{
#' \item{x2a}{the one-to-other chisquare.}
#' \item{x2b}{the maximum accumulated chisquare.}
#' \item{col1}{the column index for x2a.}
#' \item{col2}{the column index for x2b.}
#' \item{p}{the corresponding p value.}
#' }
#'
#' @references
#' Hirotsu C, Aoki S, Inada T, Kitao Y (2001) An exact test for the association 
#' between the disease and alleles at highly polymorphic loci with particular interest 
#' in the haplotype analysis. Biometrics 57:769-778
#'
#' @examples
#' \dontrun{
#' # an example from Mike Neale
#' # termed 'ugly' contingency table by Patrick Sullivan
#' y1 <- c(2,15,16,35,132,30,25,7,12,24,10,10,0)
#' y2 <- c(0, 6,31,49,120,27,15,8,14,25, 3, 9,3)
#'
#' result <- s2k(y1,y2)
#' }
#' @author Chihiro Hirotsu, Jing Hua Zhao
#' @note The lengths of y1 and y2 should be the same.
#' @keywords models

s2k <- function(y1,y2)
{
  if (length(y1)!=length(y2)) stop ("wrong number of elements")
  tablen<-length(y1)
  data<-c(y1,y2)
  x2a<-x2b<-0
  col1<-col2<-1
  p<-0
  z<-.C("x22k",data=as.integer(array(data)), tablen=as.integer(tablen),
        x2a=as.double(x2a), x2b=as.double(x2b), col1=as.integer(col1),
        col2=as.integer(col2), p=as.double(p),PACKAGE="gap")
  c1<-z$col1
  c2<-z$col2
  cat("\nthe maximum accumulated table below and above",c1,"\n")
  a<-sum(y1[1:c1])
  b<-sum(y1[-(1:c1)])
  c<-sum(y2[1:c1])
  d<-sum(y2[-(1:c1)])
  cat(a,b,a+b,"\n",c,d,c+d,"\n",a+c,b+d,a+b+c+d,"\n")
  cat("\nthe 1-to-other table with and without column",c2,"\n")
  a<-y1[c2]
  b<-sum(y1[-c2])
  c<-y2[c2]
  d<-sum(y2[-c2])
  cat(a,b,a+b,"\n",c,d,c+d,"\n",a+c,b+d,a+b+c+d,"\n")
  cat("\n")  

  list(x2a=z$x2a,x2b=z$x2b,col1=z$col1,col2=z$col2,p=z$p)
}
