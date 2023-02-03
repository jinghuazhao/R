#' Disease prevalences in cases and controls
#'
#' @param model disease model (one of "multiplicative","additive","recessive","dominant","overdominant").
#' @param GRR genotype relative risk.
#' @param p1 disease allele frequency.
#' @param K disease prevalence in the whole population.
#'
#' @details
#' KCC calculates disease prevalences in cases and controls for a given genotype relative risk,
#' allele frequency and prevalencen of the disease in the whole population. It is used by tscc
#' and pbsize2.
#' 
#' @export
#' @return
#' A list of two elements:
#' - pprime prevlence in cases.
#' - p prevalence in controls.
#'

KCC <- function(model,GRR,p1,K)
# 6-6-2018 MRC-Epid JHZ
{
   model.idx <- charmatch(model,c("multiplicative","additive","recessive","dominant","overdominant"))
   if(is.na(model.idx)) stop("Invalid model type")
   if(model.idx == 0) stop("Ambiguous model type")
   multiplicative <- c(1,GRR,GRR*GRR)
   additive <- c(1,GRR,2*GRR-1)
   recessive <- c(1,1,GRR)
   dominant <- c(1,GRR,GRR)
   overdominant <- c(GRR,1,GRR)
   f <- switch(model.idx,multiplicative,additive,recessive,dominant,overdominant)
   scale <- K/(f[1]*(1-p1)^2+f[2]*2*p1*(1-p1)+f[3]*p1^2)
   f <- f*scale
#  if(f[3]>1) stop("misspecified model")
   pprime <- (f[3]*p1^2+f[2]*p1*(1-p1))/K
   p <- ((1-f[3])*p1^2+(1-f[2])*p1*(1-p1))/(1-K)
   invisible(list(pprime=pprime,p=p))
}
