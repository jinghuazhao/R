#' Haplotype trend regression
#'
#' @param y a vector of phenotype.
#' @param x a haplotype table.
#' @param n.sim the number of permutations.
#'
#' @details
#' Haplotype trend regression (with permutation)
#'
#' @return
#' The returned value is a list containing:
#' - f the F statistic for overall association.
#' - p the p value for overall association.
#' - fv the F statistics for individual haplotypes.
#' - pi the p values for individual haplotypes.
#' 
#' @references
#' \insertRef{zaykin02}{gap}
#'
#' \insertRef{xie05}{gap}
#'
#' @seealso [`hap.score`]
#'
#' @examples
#' \dontrun{
#' # 26-10-03
#' # this is now part of demo
#' test2<-read.table("test2.dat")
#' y<-test2[,1]
#' x<-test2[,-1]
#' y<-as.matrix(y)
#' x<-as.matrix(x)
#' htr.test2<-htr(y,x)
#' htr.test2
#' htr.test2<-htr(y,x,n.sim=10)
#' htr.test2
#'
#' # 13-11-2003
#' require(gap.datasets)
#' data(apoeapoc)
#' apoeapoc.gc<-gc.em(apoeapoc[,5:8])
#' y<-apoeapoc$y
#' for(i in 1:length(y)) if(y[i]==2) y[i]<-1
#' htr(y,apoeapoc.gc$htrtable)
#'
#' # 20-8-2008
#' # part of the example from useR!2008 tutorial by Andrea Foulkes
#' # It may be used beyond the generalized linear model (GLM) framework
#' HaploEM <- haplo.em(Geno,locus.label=SNPnames)
#' HapMat <- HapDesign(HaploEM)
#' m1 <- lm(Trait~HapMat)
#' m2 <- lm(Trait~1)
#' anova(m2,m1)
#' }
#'
#' @author Dimitri Zaykin, Jing Hua Zhao
#' @note adapted from emgi.cpp, a pseudorandom number seed will be added on.
#' @keywords regression

htr <- function(y, x, n.sim=0)
{
   mlr <- function(y,x)
   {
      N <- length(y)
      L <- dim(x)[2]
      l <- L - 1
      x1 <- cbind(rep(1,N),x[,1:l])
      b1 <- MASS::ginv(t(x1) %*% x1) %*% t(x1) %*% y
      y2 <- sum(y)^2
      redss <- y2 / N
      dtssc <- t(y) %*% y - redss
      dssm <- t(b1) %*% t(x1) %*% y
      reg.ss <- dssm - redss
      reg.df <- L - 1
      err.df <- (N - 1) - reg.df
      err.ss <- dtssc - reg.ss
      reg.ms <- reg.ss / reg.df
      err.ms <- err.ss / err.df
      fstat <- reg.ms / err.ms
      pv.fstat <- 1- pf(fstat, reg.df, err.df)
      fv <- rep(1,L)
      pi <- rep(1,L)
      for (i in 1:L)
      {
          x2 <- cbind(rep(1,N),x[,i])
          b2 <- MASS::ginv(t(x2) %*% x2) %*% t(x2) %*% y
          dssm <- t(b2) %*% t(x2) %*% y
          reg.ss <- dssm - redss
          reg.df <- dim(x2)[2] - 1
          err.df <- (N - 1) - reg.df
          err.ss <- dtssc - reg.ss
          reg.ms <- reg.ss / reg.df
          err.ms <- err.ss / err.df
          b.fstat <- reg.ms / err.ms
          b.pv <- 1 - pf(b.fstat, reg.df, err.df)
          fv[i] <- b.fstat
          pi[i] <- b.pv
      }
      list(f=fstat,p=pv.fstat,fv=fv,pi=pi)
   }

   N <- length(y)
   L <- dim(x)[2]
   z0 <- mlr(y,x)
   if (n.sim==0) list(f=z0$f,p=z0$p,fv=z0$fv,pi=z0$pi)
   else
   {
      p <- 0
      pi <- rep(0,L)
      for (i in 1:n.sim)
      {
          rand.ord <- order(runif(N))
          y <- y[rand.ord]
          z <- htr(y,x)
          if (z$f >= z0$f)  p <- p + 1
          for (j in 1:L) if (z$fv[j] >= z0$fv[j]) pi[j] <- pi[j] + 1
      }
      p <- p / n.sim
      for (j in 1:L) pi[j] <- pi[j] / n.sim
      list(f=z0$f,p=p,fv=z0$fv,pi=pi)
   }
}
