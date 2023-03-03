#' Heritability estimation according to twin correlations
#'
#' @param mzDat a data frame for monzygotic twins (MZ).
#' @param dzDat a data frame for dizygotic twins (DZ).
#' @param rmz correlation for MZ twins.
#' @param rdz correlation for DZ twins.
#' @param nmz sample size for MZ twins.
#' @param ndz sample size for DZ twins.
#' @param selV names of variables for twin and cotwin.
#'
#' @details
#' Given MZ/DZ data or their correlations and sample sizes, it obtains
#' heritability and variance estimates under an ACE model as in
#' \doi{10.1038/s41562-023-01530-y} and \insertCite{keeping95;textual}{gap}.
#'
#' @export
#' @return
#' A data.frame with variables h2, c2, e2, vh2, vc2, ve2.
#'
#' @examples
#' \dontrun{
#  # simulated BMI data for men and women with bootstap
#' library(mvtnorm)
#' set.seed(12345)
#' mzm <- as.data.frame(rmvnorm(195, c(22.75,22.75),
#'                      matrix(2.66^2*c(1, 0.67, 0.67, 1), 2)))
#' dzm <- as.data.frame(rmvnorm(130, c(23.44,23.44),
#'                      matrix(2.75^2*c(1, 0.32, 0.32, 1), 2)))
#' mzw <- as.data.frame(rmvnorm(384, c(21.44,21.44),
#'                      matrix(3.08^2*c(1, 0.72, 0.72, 1), 2)))
#' dzw <- as.data.frame(rmvnorm(243, c(21.72,21.72),
#'                      matrix(3.12^2*c(1, 0.33, 0.33, 1), 2)))
#' selVars <- c('bmi1','bmi2')
#' names(mzm) <- names(dzm) <- names(mzw) <- names(dzw) <- selVars
#' ACE_CI <- function(mzData,dzData,n.sim=5,selV=NULL,verbose=TRUE)
#' {
#'   ACE_obs <- h2_mzdz(mzDat=mzData,dzDat=dzData,selV=selV)
#'   cat("\n\nheritability according to correlations\n\n")
#'   print(format(ACE_obs,digits=3),row.names=FALSE)
#'   nmz <- nrow(mzData)
#'   ndz <- nrow(dzData)
#'   r <- data.frame()
#'   for(i in 1:n.sim)
#'   {
#'     cat("\rRunning # ",i,"/", n.sim,"\r",sep="")
#'     sampled_mz <- sample(1:nmz, replace=TRUE)
#'     sampled_dz <- sample(1:ndz, replace=TRUE)
#'     mzDat <- mzData[sampled_mz,]
#'     dzDat <- dzData[sampled_dz,]
#'     ACE_i <- h2_mzdz(mzDat=mzDat,dzDat=dzDat,selV=selV)
#'     if (verbose) print(ACE_i)
#'     r <- rbind(r,ACE_i)
#'   }
#'   m <- apply(r,2,mean,na.rm=TRUE)
#'   s <- apply(r,2,sd,na.rm=TRUE)
#'   allr <- data.frame(mean=m,sd=s,lcl=m-1.96*s,ucl=m+1.96*s)
#'   print(format(allr,digits=3))
#' }
#' ACE_CI(mzm,dzm,n.sim=500,selV=selVars,verbose=FALSE)
#' ACE_CI(mzw,dzw,n.sim=500,selV=selVars,verbose=FALSE)
#' }
#'
#' @references
#' \insertRef{elks12}{gap}
#'
#' \insertAllCited{}
#' @keywords htest

h2_mzdz <- function(mzDat=NULL,dzDat=NULL,rmz=NULL,rdz=NULL,nmz=NULL,ndz=NULL,selV=NULL)
{
  if(!is.null(mzDat))
  {
    r1 <- cor(mzDat[selV[[1]]],mzDat[[selV[2]]], use="complete")
    n1 <- nrow(mzDat[stats::complete.cases(mzDat),])
  } else {
    if(is.null(rmz)|is.null(nmz)) stop("Either raw data or correlation/sample size is neeeded")
    r1 <- rmz
    n1 <- nmz
  }
  if(!is.null(dzDat))
  {
    r2 <- cor(dzDat[[selV[1]]],dzDat[[selV[2]]], use="complete")
    n2 <- nrow(dzDat[selV][stats::complete.cases(dzDat),])
  } else {
    if(is.null(rdz)|is.null(ndz)) stop("Either raw data or correlation/sample size is neeeded")
    r2 <- rdz
    n2 <- ndz
  }
  h2 <- 2 * (r1 - r2)
  c2 <- 2 * r2 - r1
  e2 <- 1 - r1
  vmz <- 1 / (n1 - 1)
  vdz <- 1 / (n2 - 1)
  vh2 <- 4 * (vmz + vdz)
  vc2 <- 4 * vdz + vmz
  ve2 <- vmz
  ACEr_est <- data.frame(h2,c2,e2,vh2,vc2,ve2)
  invisible(ACEr_est)
}
