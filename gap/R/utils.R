#' Allele recoding
#'
#' @param a1 first allele.
#' @param a2 second allele.
#' @param miss.val missing value.
#'
#' @export

allele.recode <- function (a1, a2, miss.val = NA)
{
    n <- length(a1)
    if (is.factor(a1))
        a1 <- as.character(a1)
    if (is.factor(a2))
        a2 <- as.character(a2)
    is.ch <- is.character(a1) | is.character(a2)
    if (is.ch) {
        t <- factor(c(a1, a2), exclude = miss.val)
    }
    if (!is.ch) {
        lev <- sort(unique(c(a1, a2)))
        t <- factor(c(a1, a2), levels = lev, exclude = miss.val)
    }
    allele.label <- levels(t)
    t <- as.numeric(t)
    a1 <- t[1:n]
    a2 <- t[(n + 1):(2 * n)]
    return(list(a1 = a1, a2 = a2, allele.label = allele.label))
}

#' Genotype recoding
#'
#' @param geno genotype.
#' @param miss.val missing value.
#'
#' @export

geno.recode <- function (geno, miss.val = 0)
{
    n.loci <- ncol(geno)/2
    alist <- vector("list", n.loci)
    grec <- NULL
    for (i in 1:n.loci) {
        t <- (i - 1) * 2 + 1
        tmp <- allele.recode(geno[, t], geno[, t + 1], miss.val = miss.val)
        grec <- cbind(grec, tmp$a1, tmp$a2)
        alist[[i]] <- list(allele = tmp$allele.label)
    }
    return(list(grec = grec, alist = alist))
}

#' Allele-to-genotype conversion
#'
#' @param a1 first allele.
#' @param a2 second allele.
#'
#' @export

a2g <- function(a1,a2)
{
  i <- ifelse(a1 < a2,a2,a1)
  j <- ifelse(a1 < a2,a1,a2)
  genocode <- ifelse (j==0, 0, i*(i-1)/2 + j)
  invisible (genocode)

}

#' Conversion of a genotype identifier to alleles
#'
#' @param g a genotype identifier.
#'
#' @export

g2a <- function (g)
{
    i <- 1 + floor((sqrt(8 * g + 1) - 1)/2)
    j <- g - i * (i - 1)/2
    i <- ifelse(j == 0, i - 1, i)
    j <- ifelse(j == 0, i, j)
    return(cbind(j,i))
}

#' A function to read GRM binary files
#'
#' @param prefix file root.
#' @param AllN a logical variable.
#' @param size size.
#'
#' @details Modified from GCTA documentation
#'
#' @export

ReadGRMBin <- function(prefix, AllN=FALSE, size=4)
{
  BinFileName <- paste(prefix,".grm.bin",sep="")
  NFileName <- paste(prefix,".grm.N.bin",sep="")
  IDFileName <- paste(prefix,".grm.id",sep="")
  id <- read.table(IDFileName)
  n <- dim(id)[1]
  BinFile <- file(BinFileName, "rb");
  grm <- readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  close(BinFile)
  NFile <- file(NFileName, "rb");
  if(AllN) N <- readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  else N <- readBin(NFile, n=1, what=numeric(0), size=size)
  close(NFile)
  i <- sapply(1:n, function(i) i*(i+1)/2)
  GRM <- matrix(NA,n,n)
  GRM[upper.tri(GRM,diag=TRUE)] <- grm
  GRM[lower.tri(GRM)] <- t(GRM)[lower.tri(GRM)]
  invisible(list(grm=grm, id=id, N=N, GRM=GRM))
}

#' A function to write GRM binary file
#'
#' @param prefix file root.
#' @param grm a GRM.
#' @param N Sample size.
#' @param id id.
#' @param size size.
#'
#' @export

WriteGRMBin <- function(prefix, grm, N, id, size=4)
{
  BinFileName <- paste(prefix,".grm.bin",sep="")
  NFileName <- paste(prefix,".grm.N.bin",sep="")
  IDFileName <- paste(prefix,".grm.id",sep="")
  grm.bin <- file(BinFileName, "wb")
  writeBin(grm,grm.bin,size=size)
  close(grm.bin)
  grm.N.bin <- file(NFileName, "wb")
  writeBin(N,grm.N.bin,size=size)
  close(grm.N.bin)
  write.table(id,IDFileName,col.names=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
}

#' A function to read GRM file
#' @param prefix file root.
#'
#' @export

ReadGRM <- function(prefix=51)
{
  idfile <- paste(prefix,".grm.id",sep="")
  id <- read.delim(idfile,header=FALSE,sep="\t")
  N <- dim(id)[1]
  L <- N * (N + 1) /2
  grmfile <- paste(prefix,".grm.gz",sep="")
  gz <- gzfile(grmfile)
  grm_lines <- readLines(gz)
  close(gz)
  GRM_list <- sapply(grm_lines,strsplit,"\t")
  M <- as.integer(lapply(GRM_list,"[[",3))
  grm <- as.numeric(lapply(GRM_list,"[[",4))
  GRM <- matrix(NA,N,N)
  GRM[upper.tri(GRM,diag=TRUE)] <- grm
  GRM[lower.tri(GRM)] <- t(GRM)[lower.tri(GRM)]
  invisible(list(GRM=GRM,N=M,id=id,grm=grm))
}

#' A function to write GRM file
#' @param prefix file root.
#' @param id id.
#' @param N sample size.
#' @param GRM a GRM.
#'
#' @export

WriteGRM <- function(prefix=51,id,N,GRM)
{
  M <- N
  N <- dim(GRM)[1]
  k2l <- matrix(NA,N*(N+1)/2,4)
  L <- 1
  for(i in 1:N)
  {
    for(j in 1:i)
    {
      k2l[L,] <- c(i,j,M[L],GRM[i,j])
      L <- L + 1
    }
  }
  idfile <- paste(prefix,".grm.id",sep="")
  write(t(id),file=idfile,2,sep="\t")
  grmfile <- paste(prefix,".grm.gz",sep="")
  gz <- gzfile(grmfile,"w")
  write.table(k2l,file=gz,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  close(gz)
}

#' @title Chromosomal lengths for build 36
#' @description Data are used in other functions.
#' @docType data
#' @keywords datasets
#' @format A vector containing lengths of chromosomes.
#' @details generated from GRCh.R.
#' @rdname hg18

"hg18"

#' @title Chromosomal lengths for build 37
#' @description Data are used in other functions.
#' @format A vector containing lengths of chromosomes.
#' @rdname hg19

"hg19"

#' @title Chromosomal lengths for build 38
#' @description Data are used in other functions.
#' @format A vector containing lengths of chromosomes.
#' @rdname hg38

"hg38"

#' Effect size and standard error from confidence interval
#'
#' @param ci confidence interval (CI). The delimiter between lower and upper limit is either a hyphen (-) or en dash (\enc{–}{-}).
#' @param logscale a flag indicating the confidence interval is based on a log-scale.
#' @param alpha Type 1 error.
#'
#' @details
#' Effect size is a measure of strength of the relationship between two variables in a population or parameter estimate of that population.
#' Without loss of generality, denote `m` and `s` to be the mean and standard deviation of a sample from \eqn{N(\mu,\sigma^2)}{N(mu,sigma^2}).
#' Let \eqn{z \sim N(0,1)}{z ~ N(0,1)} with cutoff point \eqn{z_\alpha}{z_alpha}, confidence limits `L`, `U` in a CI are defined as follows,
#' \deqn{
#' \begin{aligned}
#' L & = m - z_\alpha s \cr
#' U & = m + z_\alpha s
#' \end{aligned}
#' }{L = m - z_alpha s, U = m + z_alpha s}
#' \eqn{\Rightarrow}{==>} \eqn{U + L = 2 m}, \eqn{U - L=2 z_\alpha s}{U - L = 2 z_alpha s}. Consequently,
#' \deqn{
#' \begin{aligned}
#' m & = \frac{U + L}{2} \cr
#' s & = \frac{U - L}{2 z_\alpha}
#' \end{aligned}
#' }{m = (U + L)/2, s = (U - L)/(2 z_alpha)}
#' Effect size in epidemiological studies on a binary outcome is typically reported as odds ratio from a logistic regression
#' or hazard ratio from a Cox regression, \eqn{L\equiv\log(L)}{L ==> log(L)}, \eqn{U\equiv\log(U)}{U ==> log(U)}.
#'
#' @export
#' @return
#' Based on CI, the function provides a list containing estimates
#' - m effect size (log(OR))
#' - s standard error
#' - direction a decrease/increase (-/+) sign such that `sign(m)`=-1, 0, 1, is labelled "-", "0", "+", respectively as in PhenoScanner.
#'
#' @examples
#' # rs3784099 and breast cancer recurrence/mortality
#' ms <- ci2ms("1.28-1.72")
#' print(ms)
#' # Vector input
#' ci2 <- c("1.28-1.72","1.25-1.64")
#' ms2 <- ci2ms(ci2)
#' print(ms2)

ci2ms <- function(ci,logscale=TRUE,alpha=0.05)
{
  lci <- strsplit(gsub("\u2013","-",ci),"-")
  l <- as.numeric(lapply(lci,"[",1))
  u <- as.numeric(lapply(lci,"[",2))
  d <- abs(qnorm(alpha/2))
  if (!logscale) {
     s <- (u-l)/2/d
     m <- (l+u)/2
  } else {
     s <- (log(u)-log(l))/2/d
     m <- (log(l)+log(u))/2
  }
  direction <- sapply(m,function(x) {if(sign(x)==-1) "-" else if(sign(x)==0) "0" else "+"})
  invisible(list(m=m,s=s,direction=direction))
}
