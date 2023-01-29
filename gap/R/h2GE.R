#' Heritability and its variance when there is an environment component
#'
#' @param V Variance estimates.
#' @param VCOV Variance-covariance matrix.
#' @param verbose Detailed output.
#' @export
#' @return A list of phenotypic variance/heritability/GxE interaction esimates and their variances.

h2GE <- function(V, VCOV, verbose=TRUE)
 {
   n <- length(V)
   if (n==3)
   {
     VG <- V[1]
     VGE <- V[2]
     Ve <- V[3]
     Vp <- VG + VGE + Ve
     VVG <- VCOV[1, 1]
     VVGE <- VCOV[2, 2]
     VVe <- VCOV[3, 3]
     cVGVGE <- VCOV[2, 1]
     cVGVe <- VCOV[3, 1]
     cVGEVe <- VCOV[3, 2]
     s <- VVG + VVGE + VVe
     s12 <- 2 * cVGVGE
     s13 <- 2 * cVGVe
     s23 <- 2 * cVGEVe
     VVp <- s + s12 + s13 + s23
     cVpVG <- VVG + cVGVGE + cVGVe
     cVpVGE <- cVGVGE + VVGE + cVGEVe
     h2G <- VG/Vp
     Varh2G <- VR(VG, VVG, Vp, VVp, cVpVG)
     h2GE <- VGE/Vp
     Varh2GE <- VR(VGE, VVGE, Vp, VVp, cVpVGE)
     if (verbose) {
         cat("The old function results: \n",
             "Vp =", Vp, "SE =", sqrt(VVp), "\n",
             "h2G =", h2G, "SE =", sqrt(Varh2G), "\n",
             "h2GE =", h2GE, "SE =", sqrt(Varh2GE), "\n")
     }
   }
   VarV <- diag(VCOV)
   # phenotypic variance
   P <- sum(V)
   VarP <- sum(VCOV)
   # heritabilities for individual components
   GGE <- V[-n]
   s <- vector()
   for(j in 1:(n-1)) s[j] <- VR(GGE[j],VarV[j],P,VarP,sum(VCOV[j,]))
   G <- GGE[1]
   h2G <- G/P
   Varh2G <- s[1]
   GE <- GGE[-1]
   h2GE <- GE/P
   Varh2GE <- s[-1]
   # total interaction variances/heritability
   m <- -c(1,n)
   sGE <- sum(GE)
   VarsGE <- sum(VCOV[m,m])
   CovsGEP <- sum(VCOV[1,m])+VarsGE+sum(VCOV[n,m])
   h2sGE <- sGE/P
   Varh2sGE <- VR(sGE,VarsGE,P,VarP,CovsGEP)
   # genotypic variance/heritability including interactions
   sGGE <- sum(GGE)
   VarsGGE <- sum(VCOV[-n,-n])
   CovsGGE <- sum(VCOV[n,-n])
   h2sGGE <- sGGE/P
   Varh2sGGE <- VR(sGGE,VarsGGE,P,VarP,VarsGGE+CovsGGE)
   if (verbose) {
      cat("\nVariance, SE\n")
      source <- c("G",paste("GxE",1:(n-2),sep=""),"E")
      print(data.frame(V,SE=sqrt(diag(VCOV)),row.names=source))
      cat("\nV(G) (SE) =", sGGE, sqrt(VarsGGE),
          "including V(GE) (SE) =", sGE, sqrt(VarsGE),
          "\nV(P) (SE) =", P, sqrt(VarP), "\n")
      cat("\nHeritability, SE\n")
      h2 <- data.frame(h2=c(h2G,h2GE),SE=sqrt(s),row.names=source[-n])
      print(h2)
      cat("\nsum of h2GE (SE) =", h2sGE, sqrt(Varh2sGE),
          "\nsum of h2G and h2GE (SE) =", h2sGGE, sqrt(Varh2sGGE), "\n")
   }
   invisible(list(G=G, VarG=VarV[1],
                  GE=GE, VarGE=VarV[m],
                  P=P, VarP=VarP,
                  sGE=sGE, VarsGE=VarsGE,
                  sGGE=sGGE, VarsGGE=VarsGGE,
                  h2G=h2G, Varh2G=Varh2G,
                  h2GE=h2GE, Varh2GE=Varh2GE,
                  h2sGE=h2sGE, Varh2sGE=Varh2sGE,
                  h2sGGE=h2sGGE, Varh2sGGE=Varh2sGGE))
}
