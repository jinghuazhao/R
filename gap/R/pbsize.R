#' Power for population-based association design
#'
#' This function implements Long et al. (1997) statistics for population-based association
#' design. This is based on a contingency table test and accurate level of significance can be
#' obtained by Fisher's exact test.
#'
#' @param kp population disease prevalence.
#' @param gamma genotype relative risk assuming multiplicative model.
#' @param p frequency of disease allele.
#' @param alpha type I error rate.
#' @param beta type II error rate.
#'
#' @export
#' @return The returned value is scaler containing the required sample size.
#'
#' @references
#' Scott WK, Pericak-Vance MA, et al. (1997). Genetic analysis of complex diseases. Science 275: 1327.
#'		
#' Long AD, Grote MN, Langley CH (1997). Genetic analysis of complex traits. Science 275: 1328.
#'
#' Rosner B (2000). Fundamentals of Biostatistics, 5th Edition, Duxbury.
#'
#' Armitage P, Colton T (2005). Encyclopedia of Biostatistics, 2nd Edition, Wiley.
#'
#' @seealso \code{\link[gap]{fbsize}}
#'
#' @examples
#' kp <- c(0.01,0.05,0.10,0.2)
#' models <- matrix(c(
#'     4.0, 0.01,
#'     4.0, 0.10,
#'     4.0, 0.50, 
#'     4.0, 0.80,
#'     2.0, 0.01,
#'     2.0, 0.10,
#'     2.0, 0.50,
#'     2.0, 0.80,
#'     1.5, 0.01,    
#'     1.5, 0.10,
#'     1.5, 0.50,
#'     1.5, 0.80), ncol=2, byrow=TRUE)
#' outfile <- "pbsize.txt"
#' cat("gamma","p","p1","p5","p10","p20\n",sep="\t",file=outfile)
#' for(i in 1:dim(models)[1])
#' {
#'   g <- models[i,1]
#'   p <- models[i,2]
#'   n <- vector()
#'   for(k in kp) n <- c(n,ceiling(pbsize(k,g,p)))
#'   cat(models[i,1:2],n,sep="\t",file=outfile,append=TRUE)
#'   cat("\n",file=outfile,append=TRUE)
#' } 
#' table5 <- read.table(outfile,header=TRUE,sep="\t")
#' unlink(outfile)
#'
#' # Alzheimer's disease
#' g <- 4.5
#' p <- 0.15
#' alpha <- 5e-8
#' beta <- 0.2
#' z1alpha <- qnorm(1-alpha/2)   # 5.45
#' z1beta <- qnorm(1-beta)
#' q <- 1-p
#' pi <- 0.065                   # 0.07 and zbeta generate 163
#' k <- pi*(g*p+q)^2
#' s <- (1-pi*g^2)*p^2+(1-pi*g)*2*p*q+(1-pi)*q^2
#' # LGL formula
#' lambda <- pi*(g^2*p+q-(g*p+q)^2)/(1-pi*(g*p+q)^2)
#' # mine
#' lambda <- pi*p*q*(g-1)^2/(1-k)
#' n <- (z1alpha+z1beta)^2/lambda
#' cat("\nPopulation-based result: Kp =",k, "Kq =",s, "n =",ceiling(n),"\n")
#'
#' @author Jing Hua Zhao
#' extracted from rm.c.
#' @keywords misc

pbsize <- function (kp, gamma=4.5, p=0.15, alpha=5e-8, beta=0.2)
# population-based sample size
# alpha=5e-8, beta=0.8, alpha would give 5% genome-wide significance level
# x2alpha = 29.72 (Q=29.7168)
#
# lambda is the NCP from the marginal table
# pi is the pr(Affected|aa)
{
  z1alpha <- qnorm(1-alpha/2)
  zbeta <- qnorm(beta)
  q <- 1-p
  pi <- kp/(gamma*p+q)^2
  lambda <- pi*p*q*(gamma-1)^2/(1-kp)
  n <- (z1alpha-zbeta)^2/lambda
  n
}
