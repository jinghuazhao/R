#' Heritability and its variance
#'
#' @md
#' @param V Variance estimates.
#' @param VCOV Variance-covariance matrix.
#' @param verbose Detailed output.
#' @export
#' @return A list of phenotypic variance/heritability estimates and their variances.

h2G <- function(V,VCOV,verbose=TRUE)
{
  VG <- V[1]
  Ve <- V[2]
  Vp <- VG + Ve
  VVG <- VCOV[1,1]
  VVe <- VCOV[2,2]
  cVGVe <- VCOV[2,1]
  h2G <- VG / Vp
  VVp <- VVG + VVe + 2 * cVGVe
  cVpVG <- VVG + cVGVe
  Varh2G <- VR(VG,VVG,Vp,VVp,cVpVG)
  if (verbose) {
     cat("Vp =",Vp,"SE =",sqrt(VVp),"\n")
     cat("h2G =",h2G,"SE =",sqrt(Varh2G),"\n")
  }
  invisible(list(Vp=Vp,VVp=VVp,h2G=h2G,Varh2G=Varh2G))
}

