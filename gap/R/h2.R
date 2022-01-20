k <- function(r,N,adjust=TRUE)
{
  r2 <- r^2
  n <- N-1
  k1 <- ifelse(adjust,r-r*(1-r2)/2/n,r)
  k2 <- (1-r2)^2/n*(1+11*r2/2/n)
  invisible(c(k1,k2))
}

#' Heritability estimation according to twin correlations
#'
#' Heritability and variance estimation according to twin pair correlations.
#'
#' @param mzDat a data frame for monzygotic twins (MZ).
#' @param dzDat a data frame for dizygotic twins (DZ).
#' @param rmz correlation for MZ twins.
#' @param rdz correlation for DZ twins.
#' @param nmz sample size for MZ twins.
#' @param ndz sample size for DZ twins.
#' @param selV names of variables for twin and cotwin.
#'
#' @details The example section shows how to obtain bootstrap 95\%CI.
#'
#' @export
#' @return
#' The returned value is a matrix containing heritability and their variance estimations
#' for "h2","c2","e2","vh","vc","ve".
#'
#' @references
#' Keeping ES. Introduction to Statistical Inference, Dover Pulications, Inc. 1995
#'
#' @examples
#' \dontrun{
#'
#' ACE_CI <- function(mzData,dzData,n.sim=5,selV=NULL,verbose=TRUE)
#' {
#' ACEr_twinData <- h2(mzDat=mzData,dzDat=dzData,selV=selV)
#' print(ACEr_twinData)
#'
#' nmz <- dim(mzData)[1]
#' ndz <- dim(dzData)[1]
#' a <- ar <- vector()
#' set.seed(12345)
#' for(i in 1:n.sim)
#' {
#'   cat("\rRunning # ",i,"/", n.sim,"\r",sep="")
#'   sampled_mz <- sample(1:nmz, replace=TRUE)
#'   sampled_dz <- sample(1:ndz, replace=TRUE)
#'   mzDat <- mzData[sampled_mz,]
#'   dzDat <- dzData[sampled_dz,]
#'   ACEr_i <- h2(mzDat=mzDat,dzDat=dzDat,selV=selV)
#'   if(verbose) print(ACEr_i)
#'   ar <- rbind(ar,ACEr_i)
#' }
#' cat("\n\nheritability according to correlations\n\n")
#' ar <- as.data.frame(ar)
#' m <- mean(ar,na.rm=TRUE)
#' s <- sd(ar,na.rm=TRUE)
#' allr <- data.frame(mean=m,sd=s,lcl=m-1.96*s,ucl=m+1.96*s)
#' print(allr)
#' }
#'
#' selVars <- c('bmi1','bmi2')
#'
#' library(mvtnorm)
#' n.sim <- 500
#' cat ("\nThe first study\n\n")
#' mzm <- as.data.frame(rmvnorm(195, c(22.75,22.75),
#'                      matrix(2.66^2*c(1, 0.67, 0.67, 1), 2)))
#' dzm <- as.data.frame(rmvnorm(130, c(23.44,23.44),
#'                      matrix(2.75^2*c(1, 0.32, 0.32, 1), 2)))
#' mzw <- as.data.frame(rmvnorm(384, c(21.44,21.44),
#'                      matrix(3.08^2*c(1, 0.72, 0.72, 1), 2)))
#' dzw <- as.data.frame(rmvnorm(243, c(21.72,21.72),
#'                      matrix(3.12^2*c(1, 0.33, 0.33, 1), 2)))
#' names(mzm) <- names(dzm) <- names(mzw) <- names(dzw) <- c("bmi1","bmi2")
#' ACE_CI(mzm,dzm,n.sim,selV=selVars,verbose=FALSE)
#' ACE_CI(mzw,dzw,n.sim,selV=selVars,verbose=FALSE)
#'
#' }
#' @author Jing Hua Zhao
#' @keywords htest

h2_mzdz <- function(mzDat=NULL,dzDat=NULL,rmz=NULL,rdz=NULL,nmz=NULL,ndz=NULL,selV=NULL)
{
  if(!is.null(mzDat))
  {
    r1 <- cor(mzDat[selV[1]],mzDat[selV[2]], use="complete")
    n1 <- length(!is.na(c(mzDat[selV[1]],mzDat[selV[2]])))
  } else {
    if(is.null(rmz)|is.null(nmz)) stop("Either raw data or correlation/sample size is neeeded")
    r1 <- rmz
    n1 <- nmz
  }
  if(!is.null(dzDat))
  {
    r2 <- cor(dzDat[selV[1]],dzDat[selV[2]], use="complete")
    n2 <- length(!is.na(c(dzDat[selV[1]],dzDat[selV[2]])))
  } else {
    if(is.null(rdz)|is.null(ndz)) stop("Either raw data or correlation/sample size is neeeded")
    r2 <- rdz
    n2 <- ndz
  }
  kmz <- k(r1,n1)
  k1mz <- kmz[1]
  k2mz <- kmz[2]
  kdz <- k(r2,n2)
  k1dz <- kdz[1]
  k2dz <- kdz[2]
  h2 <- 2 * (k1mz - k1dz)
  vh <- 4 * (k2mz + k2dz)
  c2 <- 2 * k1dz - k1mz
  vc <- 4 * k2dz + k2mz
  e2 <- 1 - k1mz
  ve <- k2mz
  ACEr_est <- as.matrix(c(h2,c2,e2,vh,vc,ve))
  rownames(ACEr_est) <- c("h2","c2","e2","vh","vc","ve")
  invisible(t(ACEr_est))
}
