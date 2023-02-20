#' A function used to experiment with ETDT
#'
#' @details
#' It fits a Bradley-Terry Model to a squared contingency table
#' @noRd

# JH Zhao 5/1/2004, 12/3/2005

# This generates format as required by ETDT
# adapted from bt.sas (16/05/1999)
toETDT <- function(a)
{
   n <- dim(a)[1]
   nn1 <- n*(n-1)/2;
   C<-tr<-i1<-i2<-rep(0,nn1)
   ij <- matrix(rep(0,n*nn1),nrow=nn1)
   vi <- rep(0,n)
   l <- 1
   for(i in 1:(n-1))
   {
      for(j in (i+1):n)
      {
        C[l] <- a[i,j]+a[j,i]
        tr[l] <- a[i,j]
        i1[l] <- i
        i2[l] <- j
        ij[l,i] <- 1
        ij[l,j] <- -1
        l <- l + 1
      }
      vi[i] <- i
   }
   cbind(i1,i2,C,tr,ij)
}

#' Bradley-Terry model for contingency table
#'
#' @param x the data table.
#'
#' @details
#' This function calculates statistics under Bradley-Terry model. 
#'
#' @export
#' @return
#' The returned value is a list containing:
#' - y A column of 1.
#' - count the frequency count/weight.
#' - allele the design matrix.
#' - bt.glm a glm.fit object.
#' - etdt.dat a data table that can be used by ETDT.
#'
#' @references
#' \insertRef{bradley52}{gap}
#'
#' \insertRef{sham95}{gap}
#'
#' \insertRef{copeman95}{gap}
#' @seealso [`mtdt`]
#'
#' @examples
#' \dontrun{
#' x <- matrix(c(0,0, 0, 2, 0,0, 0, 0, 0, 0, 0, 0,
#'               0,0, 1, 3, 0,0, 0, 2, 3, 0, 0, 0,
#'               2,3,26,35, 7,0, 2,10,11, 3, 4, 1,
#'               2,3,22,26, 6,2, 4, 4,10, 2, 2, 0,
#'               0,1, 7,10, 2,0, 0, 2, 2, 1, 1, 0,
#'               0,0, 1, 4, 0,1, 0, 1, 0, 0, 0, 0,
#'               0,2, 5, 4, 1,1, 0, 0, 0, 2, 0, 0,
#'               0,0, 2, 6, 1,0, 2, 0, 2, 0, 0, 0,
#'               0,3, 6,19, 6,0, 0, 2, 5, 3, 0, 0,
#'               0,0, 3, 1, 1,0, 0, 0, 1, 0, 0, 0,
#'               0,0, 0, 2, 0,0, 0, 0, 0, 0, 0, 0,
#'               0,0, 1, 0, 0,0, 0, 0, 0, 0, 0, 0),nrow=12)
#'
#' # Bradley-Terry model, only deviance is available in glm
#' # (SAS gives score and Wald statistics as well)
#' bt.ex<-bt(x)
#' anova(bt.ex$bt.glm)
#' summary(bt.ex$bt.glm)
#' }
#'
#' @author Jing Hua Zhao
#' @note Adapted from a SAS macro for data in the example section.
#' @keywords models

bt <- function(x)
{
  n <- dim(x)[1]
  if (n==1) stop("the dimension of the square table should be at least 2")
  else if (n==2) allele <- matrix(c(1,-1,-1,1),nrow=2,byrow=T)
  else if (n>2) {
     one <- rep(1,n-1)
     d <- - diag(n-1)
     allele <- cbind(one,d)
     for (i in 2:(n-1))
     {
       indx <- i:(n-1)
       left <- d[,-indx]
       right <- d[,indx]
       a <- cbind(left,one,right)
       allele <- rbind(allele,a)
     }
     allele <- rbind(allele,cbind(d,one))
  }
  y <- rep(1,n*(n-1))
  count <- x[lower.tri(x)|upper.tri(x)]
  bt.glm <- glm(y~-1+allele,weights=count,family="binomial")
  list(y=y,count=count,allele=allele,bt.glm=bt.glm,etdt.dat=toETDT(x))
}
