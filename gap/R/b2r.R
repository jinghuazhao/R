# 18-8-2009 MRC-Epid JHZ

covfun <- function(r,n)
{
   rst <- rts <- r[1]
   rsu <- rus <- r[2]
   rtu <- rut <- r[3]
   rsv <- rvs <- r[4]
   rtv <- rvt <- r[5]
   ruv <- rvu <- r[6]
   cov.r <- (0.5*rst*ruv*(rsu^2+rsv^2+rtu^2+rtv^2)+rsu*rtv+rsv*rtu
            -(rst*rsu*rsv+rts*rtu*rtv+rus*rut*ruv+rus*rvt*rvu))/n
}

#' Obtain correlation coefficients and their variance-covariances
#'
#' This function converts linear regression coefficients of phenotype on
#' single nucleotide polymorphisms (SNPs) into Pearson correlation coefficients
#' with their variance-covariance matrix. It is useful as a preliminary step 
#' for meta-analyze SNP-trait associations at a given region. Between-SNP
#' correlations (e.g., from HapMap) are required as auxiliary information.
#'
#' @param b the vector of linear regression coefficients.
#' @param s the corresponding vector of standard errors.
#' @param rho triangular array of between-SNP correlation.
#' @param n the sample size.
#'
#' @export
#' @return The returned value is a list containing:
#' \describe{
#' \item{r}{the vector of correlation coefficients}
#' \item{V}{the variance-covariance matrix of correlations}
#' }
#'
#' @references
#' Becker BJ (2004). Multivariate meta-analysis. in Tinsley HEA,
#' Brown SD (Ed.) Handbook of Applied Multivariate Statistics and
#' Mathematical Modeling (Chapter 17, pp499-525). Academic Press.
#'
#' Casella G, Berger RL (2002). Statistical Inference, 2nd Edition, Duxbury.
#'
#' Elston RC (1975). On the correlation between correlations. Biometrika 62:133-40
#'
#' @seealso \code{\link[gap]{mvmeta}}, \code{\link[gap]{LD22}}
#'
#' @examples
#' \dontrun{
#' n <- 10
#' r <- c(1,0.2,1,0.4,0.5,1)
#' b <- c(0.1,0.2,0.3)
#' s <- c(0.4,0.3,0.2)
#' bs <- b2r(b,s,r,n)
#' }
#' 
#' @author Jing Hua Zhao
#' @keywords datagen

b2r <- function(b,s,rho,n)
{
   covfun2 <- function(r1,r2,r12) invisible(r1*r2*r12^2+(r1^2+r2^2-1)*(r1*r2-2*r12))
   m <- length(b)
   t <- b/s
   t2 <- t^2
   r2 <- t2/(n-2+t2)
   r <- sqrt(r2)
   V <- matrix(NA,m,m)
   for(i in 1:m)
   {
     for(j in 1:i)
     {
        l <- min(i,j)
        u <- max(i,j)
        l <- l+(u-1)*u/2
        r1 <- r[i]
        r2 <- r[j]
        r12 <- rho[l]
        V[i,j] <- V[j,i] <- covfun2(r1,r2,r12)
     }
   }
   invisible(list(r=r,V=V/2/n))
}
