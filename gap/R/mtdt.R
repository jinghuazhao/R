#' Transmission/disequilibrium test of a multiallelic marker
#'
#' This function calculates transmission-disequilibrium statistics involving multiallelic marker.
#'
#' Inside the function are tril and triu used to obtain lower and upper triangular matrices.
#'
#' @param x the data table.
#' @param n.sim the number of simulations.
#'
#' @export
#' @return
#' It returned list contains the following components:
#' \describe{
#'   \item{SE}{Spielman-Ewens Chi-square from the observed data}
#'   \item{ST}{Stuart or score Statistic from the observed data}
#'   \item{pSE}{the simulated p value}
#'   \item{sSE}{standard error of the simulated p value}
#'   \item{pST}{the simulated p value}
#'   \item{sST}{standard error of the simulated p value}
#' }
#'
#' @references
#' Miller MB (1997) Genomic scanning and the transmission/disequilibrium test: 
#' analysis of error rates. Genet. Epidemiol. 14:851-856
#'
#' Sham PC (1997) Transmission/disequilibrium tests for multiallelic loci. 
#' Am. J. Hum. Genet. 61:774-778
#'
#' Spielman RS, Ewens WJ (1996) The TDT and other family-based tests for
#' linkage disequilibrium and association. Am. J. Hum. Genet. 59:983-989
#'
#' Zhao JH, Sham PC, Curtis D (1999) A program for the Monte Carlo evaluation 
#' of significance of the extended transmission/disequilibrium test. 
#' Am. J. Hum. Genet. 64:1484-1485
#'
#' @seealso{\code{\link[gap]{bt}}}
#'
#' @examples
#' \dontrun{
#' # Copeman et al (1995) Nat Genet 9: 80-5
#'
#' x <- matrix(c(0,0, 0, 2, 0,0, 0, 0, 0, 0, 0, 0,
#'               0,0, 1, 3, 0,0, 0, 2, 3, 0, 0, 0,
#'               2,3,26,35, 7,0, 2,10,11, 3, 4, 1,
#'               2,3,22,26, 6,2, 4, 4,10, 2, 2, 0,
#'               0,1, 7,10, 2,0, 0, 2, 2, 1, 1, 0,
#'               0,0, 1, 4, 0,1, 0, 1, 0, 0, 0, 0,
#'               0,2, 5, 4, 1,1, 0, 0, 0, 2, 0, 0,
#'               0,0, 2, 6, 1,0, 2, 0, 2, 0, 0, 0,
#'               0,3, 6,19, 6,0, 0, 2, 5, 3, 0, 0,
#'               0,0, 3, 1, 1,0, 0, 0, 1, 0, 0, 0,
#'               0,0, 0, 2, 0,0, 0, 0, 0, 0, 0, 0,
#'               0,0, 1, 0, 0,0, 0, 0, 0, 0, 0, 0),nrow=12)
#'
#' # See note to bt for the score test obtained by SAS
#'
#' mtdt(x)
#' }
#'
#' @author Mike Miller, Jing Hua Zhao
#' @keywords models

mtdt <- function(x,n.sim=0)
{
#
# lower triangular matrix 5/6/2004
#
  tril <- function(t)
  {
    a <- t
    a[upper.tri(t)] <- 0
    a
  }

#
# upper triangular matrix 5/6/2004
#
  triu <- function(t)
  {
    a <- t
    a[lower.tri(t)] <- 0
    a
  }

#
# Spielman & Ewens' statistic (diagnoal elements kept)
# SAS 20/07/1999
#
  T <- x
  T <- T-diag(diag(T))
  m <- dim(T)[1]
#
# Find non-zero off-diagonal elements of table.
#
  b <- matrix(rep(0,m*m),nrow=m)
  i <- j <- matrix(rep(0,m*(m-1)/2),ncol=1)
  for (ii in 1:m)
  {
     for (jj in 1:(ii-1)) b[ii,jj]=T[ii,jj]+T[jj,ii]
  }
  l <- 0
  for (ii in 1:(m-1))
  {
     for (jj in (ii + 1):m)
        if (b[jj,ii]>0)
        {
           l <- l + 1
           i[l] <- jj
           j[l] <- ii
        }
  }
#
# Count number of different heterozygotes.
#
  NH <- l
#
# Count number of each heterozygote.
#
  t0 <- T
  tc <- apply(t0,2,sum)
  tr <- apply(t0,1,sum)
  c0 <- (m-1)/m
  se0 <- c0*sum((tc-tr)^2/(tc+tr))
  if (n.sim>0)
  {
    C <- rep(0,NH)
    for (k in 1:NH) C[k] <- T[i[k],j[k]] + T[j[k],i[k]]
    X <- rep(0,n.sim)
    for (k in 1:n.sim)
    {
       for (L in 1:NH)
       {
          T[i[L],j[L]] <- rbinom(1,C[L],0.5) # sum(runif(C(L))<0.5)
          T[j[L],i[L]] <- C[L] - T[i[L],j[L]]
       }
       tc <- apply(T,2,sum)
       tr <- apply(T,1,sum)
       X[k] <- c0*sum((tc-tr)^2/(tc+tr))
    }
    MCp <- 0
    for (k in 1:n.sim) if (X[k]>=se0) MCp <- MCp + 1
    pSE <- MCp/n.sim
    sSE <- sqrt(pSE*(1-pSE)/n.sim)
    cat('Spielman-Ewens Chi-square and empirical p (se): ', se0, pSE, sSE, "\n")
  } else cat('Spielman-Ewens Chi-square: ', se0, "\n")
#
# Simulate tables and compute TDT chi-square statistics.
#
# Should diag(T) is kept, the statistic will be similar to Spielman-Ewens'
# se.check

  T <- x

# This is according to Mike Miller's Matlab program
# Produce inverse (IV) of variance-covariance matrix (V).  This is
# constant across repeated samples in the Monte Carlo simulation.

  V <- diag(apply(T,1,sum)+apply(T,2,sum))-tril(T)-t(triu(T))-t(tril(T))-triu(T)
  IV <- solve(V[1:(m-1),1:(m-1)])

  T0 <- T
  d0=apply(T0[1:(m-1),1:(m-1)],1,sum)-apply(T0[1:(m-1),1:(m-1)],2,sum)
  x0 <- t(d0)%*%IV%*%d0
  se.check <- x0

  T <- T - diag(diag(T))
  V <- diag(apply(T,1,sum)+apply(T,2,sum))-tril(T)-t(triu(T))-t(tril(T))-triu(T)
  IV <- solve(V[1:(m-1),1:(m-1)])

  d0=apply(T0[1:(m-1),1:(m-1)],1,sum)-apply(T0[1:(m-1),1:(m-1)],2,sum)
  st0 <- t(d0)%*%IV%*%d0
  if (n.sim>0)
  {
    C <- rep(0,NH)
    for (k in 1:NH) C[k] <- T[i[k],j[k]] + T[j[k],i[k]]
    X <- rep(0,n.sim)
    for (k in 1:n.sim)
    {
       for (L in 1:NH)
       {
          T[i[L],j[L]] <- rbinom(1,C[L],0.5) # sum(runif(C(L))<0.5)
          T[j[L],i[L]] <- C[L]-T[i[L],j[L]]
       }
       d <- apply(T[1:(m-1),1:(m-1)],1,sum)-apply(T[1:(m-1),1:(m-1)],2,sum)
       X[k] <- t(d)%*%IV%*%d
    }
    MCp <- 0
    for (k in 1:n.sim) if (X[k]>st0) MCp <- MCp + 1
    pST <- MCp/n.sim
    sST <- sqrt(pST*(1-pST)/n.sim)
    cat('Stuart Chi-square and p (se): ', st0, pST, sST,"\n")
  } else {
    cat('Stuart Chi-square ',st0, "\n")
    cat('Value of Chi-square if diagonal elements are kept: ',se.check,"\n")
  }
  if (n.sim>0) list(SE=se0,pSE=pSE,sSE=sSE,ST=st0,pST=pST,sST=sST)
  else list(SE=se0,ST=st0)
}
