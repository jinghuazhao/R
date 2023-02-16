#' Multivariate meta-analysis based on generalized least squares
#'
#' @param b the parameter estimates.
#' @param V the triangular variance-covariance matrix.
#'
#' @details
#' This function accepts a data matrix of parameter estimates and their variance-covariance matrix
#' from individual studies and obtain a generalized least squares (GLS) estimate and heterogeneity statistic.
#'
#' For instance, this would be appropriate for combining linear correlation coefficients of single
#' nucleotide polymorphisms (SNPs) for a given region.
#'
#' @export
#' @return The returned value is a list containing:
#' - d the compact parameter estimates.
#' - Psi the compact covariance-covariance matrix.
#' - X the design matrix.
#' - beta the pooled parameter estimates.
#' - cov.beta the pooled variance-covariance matrix.
#' - X2 the Chi-squared statistic for heterogeneity.
#' - df the degrees(s) of freedom.
#' - p the p value.
#'
#' @references
#' \insertRef{hartung08}{gap}
#'
#' @seealso [`metareg`]
#'
#' @examples
#' \dontrun{
#' # example 11.3 from Hartung et al.
#' #
#' b <- matrix(c(
#' 0.808, 1.308, 1.379, NA, NA,
#' NA, 1.266, 1.828, 1.962, NA,
#' NA, 1.835, NA, 2.568, NA,
#' NA, 1.272, NA, NA, 2.038,
#' 1.171, 2.024, 2.423, 3.159, NA,
#' 0.681, NA, NA, NA, NA),ncol=5, byrow=TRUE)
#'
#' psi1 <- psi2 <- psi3 <- psi4 <- psi5 <- psi6 <- matrix(0,5,5)
#'
#' psi1[1,1] <- 0.0985
#' psi1[1,2] <- 0.0611
#' psi1[1,3] <- 0.0623
#' psi1[2,2] <- 0.1142
#' psi1[2,3] <- 0.0761
#' psi1[3,3] <- 0.1215
#'
#' psi2[2,2] <- 0.0713
#' psi2[2,3] <- 0.0539
#' psi2[2,4] <- 0.0561
#' psi2[3,3] <- 0.0938
#' psi2[3,4] <- 0.0698
#' psi2[4,4] <- 0.0981
#'
#' psi3[2,2] <- 0.1228
#' psi3[2,4] <- 0.1119
#' psi3[4,4] <- 0.1790
#'
#' psi4[2,2] <- 0.0562
#' psi4[2,5] <- 0.0459
#' psi4[5,5] <- 0.0815
#'
#' psi5[1,1] <- 0.0895
#' psi5[1,2] <- 0.0729
#' psi5[1,3] <- 0.0806
#' psi5[1,4] <- 0.0950
#' psi5[2,2] <- 0.1350
#' psi5[2,3] <- 0.1151
#' psi5[2,4] <- 0.1394
#' psi5[3,3] <- 0.1669
#' psi5[3,4] <- 0.1609
#' psi5[4,4] <- 0.2381
#'
#' psi6[1,1] <- 0.0223
#'
#' V <- rbind(psi1[upper.tri(psi1,diag=TRUE)],psi2[upper.tri(psi2,diag=TRUE)],
#' psi3[upper.tri(psi3,diag=TRUE)],psi4[upper.tri(psi4,diag=TRUE)],
#' psi5[upper.tri(psi5,diag=TRUE)],psi6[upper.tri(psi6,diag=TRUE)])
#'
#' mvmeta(b,V)
#' }
#' @author Jing Hua Zhao
#' @keywords datagen

mvmeta <- function(b,V)
{
   for(p in c("magic","MASS")) {
      if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
         if (!requireNamespace(p, quietly = TRUE))
         warning(paste("mvmeta needs package `", p, "' to be fully functional; please install", sep=""))
      }
   }
   n1 <- dim(b)[1]
   n2 <- dim(b)[2]
   d <- as.vector(t(b))
   d <- d[!is.na(d)]
   X <- vector()
   Vi <- matrix(NA,n2,n2)
   for (i in 1:n1)
   {
      dz <- !is.na(b[i,])
      dx <- diag(ifelse(dz,1,NA))
      X <- rbind(X,dx[dz,])
      Vi[upper.tri(Vi,diag=TRUE)] <- V[i,]
      if (i==1) Psi <- magic::adiag(Vi[dz,dz])
      else Psi <- magic::adiag(Psi,Vi[dz,dz])
   }
   dl <- length(d)
   for (i in 1:dl) Psi[i:dl,i] <- Psi[i,i:dl]
#  Psi <- replace(Psi,is.na(Psi),0)
   cpd <- t(X) %*% MASS::ginv(Psi) %*% X
   beta <- MASS::ginv(cpd) %*% t(X) %*% MASS::ginv(Psi) %*% d
   cov.beta <- MASS::ginv(cpd)
   X2 <- t(d) %*% MASS::ginv(Psi) %*% d - t(beta) %*% cpd %*% beta
   df <- length(d)-n2
   p <- 1-pchisq(X2,df)

   invisible (list(d=d,Psi=Psi,X=X,beta=beta,cov.beta=cov.beta,X2=X2,df=df,p=p))
}
