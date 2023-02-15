#' Hardy-Weinberg equlibrium test for a multiallelic marker
#'
#' @param data A rectangular data containing the genotype, or an array of genotype counts.
#' @param data.type An option taking values "allele", "genotype", "count"  if data is alleles, genotype or genotype count.
#' @param yates.correct A flag indicating if Yates' correction is used for Pearson \eqn{\chi^2}{chi-squared} statistic.
#' @param miss.val A list of missing values.
#'
#' @details
#' Hardy-Weinberg equilibrium test.
#'
#' This function obtains Hardy-Weinberg equilibrium test statistics. It can
#' handle data coded as allele numbers (default), genotype identifiers (by
#' setting data.type="genotype") and counts corresponding to individual genotypes
#' (by setting data.type="count") which requires that genotype counts for all n(n+1) possible genotypes,
#' with n being the number of alleles.
#'
#' For highly polymorphic markers when asymptotic results do not hold, please resort to [`hwe.hardy`].
#'
#' @export
#' @return
#' The returned value is a list containing:
#' - allele.freq Frequencies of alleles.
#' - x2 Pearson \eqn{\chi^2}{chi-square}.
#' - p.x2 p value for \eqn{\chi^2}{chi-square}.
#' - lrt Log-likelihood ratio test statistic.
#' - p.lrt p value for lrt.
#' - df Degree(s) of freedom.
#' - rho \eqn{\sqrt{\chi^2/N}}{sqrt{chi-square/N}} the contingency table coefficient.
#'
#' @seealso [hwe.hardy]
#'
#' @examples
#' \dontrun{
#' a <- c(3,2,2)
#' a.out <- hwe(a,data.type="genotype")
#' a.out
#' a.out <- hwe(a,data.type="count")
#' a.out
#' require(haplo.stats)
#' data(hla)
#' hla.DQR <- hwe(hla[,3:4])
#' summary(hla.DQR)
#' # multiple markers
#' s <- vector()
#' for(i in seq(3,8,2))
#' {
#'   hwe_i <- hwe(hla[,i:(i+1)])
#'   s <- rbind(s,hwe_i)
#' }
#' s
#' }
#'
#' @author Jing Hua Zhao
#' @keywords htest

hwe <- function(data, data.type="allele", yates.correct=FALSE, miss.val=0)
{
  data.int <- charmatch(data.type,c("allele","genotype","count"))
  if (is.na(data.int)) stop("Invalid data type")
  if (data.int==0) stop("Ambiguous data type")
  x2 <- lrt <- 0
  if (data.int<3)
  {
     # remove individuals with missing data
     # obtain maximum number of alleles and table of genotype counts
     if (data.int==1)
     {
         a1 <- data[,1]
         a2 <- data[,2]
     }
     else
     {
         tmp <- g2a(data)
         a1 <- tmp[,1]
         a2 <- tmp[,2]
         # need recoding since some may be overcoded.
     }
     tmp <- allele.recode(a1,a2,miss.val)
     a1 <- tmp$a1
     a2 <- tmp$a2
     n.allele <- length(tmp$allele.label)
     geno.table <- table(a2g(a1,a2))
     n <- sum(geno.table, na.rm=T)
     n2 <- 2.0 * n
     allele.freq <- table(c(a1,a2)) / n2
     dimnames(allele.freq)[[1]] <- tmp$allele.label
     geno.obs <- as.numeric(names(geno.table))
     geno.len <- length(geno.table)
     n.genotype <- n.allele * (n.allele + 1) /2
     # nonobserved contributes data
     geno.all <- array(0,n.genotype)
     for(i in 1:geno.len)
     {
        j <- geno.obs[i]
        geno.all[j] <- geno.table[i]
     }
     for (i in 1:n.genotype)
     {
         o <- geno.all[i]
         z <- g2a(i)
         l <- z[1]
         u <- z[2]
         e <- ifelse(l==u,1,2) * allele.freq[l] * allele.freq[u] * n
         x2 <- x2 + ifelse(yates.correct, (abs (o - e) - 0.5)^2 / e, (o - e)^2 / e)
         if (o>0) lrt <- lrt + o * log(o / e)
     }
  }
  else
  {
      n.genotype <- length(data)
      n <- sum(data)
      e <- array(0,n.genotype)
      n.allele <- (sqrt(1 + 8 * n.genotype) - 1) / 2
      allele.freq <- array(0,n.allele)
      k <- 0
      for (i in 1:n.allele)
      {
          for (j in 1:i)
          {
              k <- k + 1
              allele.freq[i] <- allele.freq[i] + data[k]
              allele.freq[j] <- allele.freq[j] + data[k]
          }
      }
      allele.freq <- allele.freq / 2.0 / n
      k <- 0
      for (i in 1:n.allele)
          for (j in 1:i)
          {
              k <- k + 1
              e[k] <- ifelse(i==j,1,2) * allele.freq[i] * allele.freq[j] * n
          }
      for (i in 1:n.genotype)
      {
          o <- data[i]
          if (e[i]>0)
          {
            x2 <- x2 + ifelse(yates.correct, (abs(o - e[i]) - 0.5)^2 / e[i], (o - e[i])^2 / e[i])
            if (o > 0) lrt <- lrt + o * log(o / e[i])
          }
      }
  }
  df <- n.genotype - n.allele
  lrt <- 2.0 * lrt
  rho <- sqrt(x2 / n)
  cat("Pearson x2=",round(x2,3),", df=",df,", p=",1-pchisq(x2,df),sep="\t")
  cat("\n")
  list (allele.freq=allele.freq,x2=x2, p.x2=1-pchisq(x2,df),
        lrt=lrt, p.lrt=1-pchisq(lrt,df), df=df, rho=rho)
}

# 08-02-2004 Working but miss.value needs to be added later
# 09-02-2004 Add using genotype counts
# 11-02-2004 Works ok in all three modes
# 17-02-2004 Add is.miss function
# 04-10-2004 tidy is.count/is.genotype with data.type and document allele.freq
# 06-11-2004 simply code
# 08-11-2004 debug, refer to code of August 2004 and done
# 09-11-2004 correct allele lables for allele frequencies
# 10-11-2004 fix lrt statistic and p.lrt
# 16-01-2005 fix n.genotype
