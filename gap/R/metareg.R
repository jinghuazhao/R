#' Fixed and random effects model for meta-analysis
#'
#' @param data Data frame to be used.
#' @param N Number of studies.
#' @param verbose A control for screen output.
#' @param prefixb Prefix of estimate; default value is "b".
#' @param prefixse Prefix of standard error; default value is "se".
#' The function accepts a wide format data with estimates as \eqn{b1,...,bN}
#' and standard errors as \eqn{se1,...,seN}. More generally, they can be specified
#' by prefixes in the function argument.
#'
#' @details
#' Given \eqn{k=n} studies with \eqn{b_1, ..., b_N} being \eqn{\beta}'s and
#' \eqn{se_1, ..., se_N} standard errors from regression, the fixed effects
#' model uses inverse variance weighting such that \eqn{w_1=1/se_1^2}, ...,
#' \eqn{w_N=1/se_N^2} and the combined \eqn{\beta} as the weighted average,
#' \eqn{\beta_f=(b_1*w_1+...+b_N*w_N)/w}, with \eqn{w=w_1+...+w_N} 
#' being the total weight, the se for this estimate is \eqn{se_f=\sqrt{1/w}}.
#' A normal z-statistic is obtained as \eqn{z_f=\beta_f/se_f}, and the
#' corresponding p value \eqn{p_f=2*pnorm(-abs(z_f))}. For the random effects
#' model, denote \eqn{q_w=w_1*(b_1-\beta_f)^2+...+w_N*(b_N-\beta_f)^2}
#' and \eqn{dl=max(0,(q_w-(k-1))/(w-(w_1^2+...+w_N^2)/w))}, corrected
#' weights are obtained such that \eqn{{w_1}_c=1/(1/w_1+dl)}, ..., 
#' \eqn{{w_N}_c=1/(1/w_N+dl)}, totaling \eqn{w_c={w_1}_c+...+{w_N}_c}.
#' The combined \eqn{\beta} and se are then \eqn{\beta_r=(b_1*{w_1}_c+...+b_N*{w_N}_c)/w_c}
#' and \eqn{se_r=\sqrt(1/w_c)}{se_r=sqrt(1/wc)}, leading to a z-statistic
#' \eqn{z_r=\beta_r/se_r} and a p-value \eqn{p_r=2*pnorm(-abs(z_r))}. Moreover, a
#' p-value testing for heterogeneity is \eqn{p_{heter}=pchisq(q_w,k-1,lower.tail=FALSE)}.
#'
#' @export
#' @return
#' The returned value is a data frame with the following variables:
#' - p_f P value (fixed effects model).
#' - p_r P value (random effects model).
#' - beta_f regression coefficient.
#' - beta_r regression coefficient.
#' - se_f standard error.
#' - se_r standard error.
#' - z_f z value.
#' - z_r z value.
#' - p_heter heterogeneity test p value.
#' - i2 \eqn{I^2}{I^2} statistic.
#' - k No of tests used.
#' - eps smallest double-precision number.
#'
#' @references
#' \insertRef{higgins03}{gap}
#'
#' @examples
#' \dontrun{
#' abc <- data.frame(chromosome=1,rsn='abcd',startpos=1234,
#'                   b1=1,se1=2,p1=0.1,b2=2,se2=6,p2=0,b3=3,se3=8,p3=0.5)
#' metareg(abc,3)
#' abc2 <- data.frame(b1=c(1,2),se1=c(2,4),b2=c(2,3),se2=c(4,6),b3=c(3,4),se3=c(6,8))
#' print(metareg(abc2,3))
#' }
#'
#' @author Shengxu Li, Jing Hua Zhao
#' @note Adapted from a SAS macro, 23-7-2009 MRC-Epid JHZ
#' @keywords models

metareg <- function(data,N,verbose="Y",prefixb="b",prefixse="se")
{
   M <- dim(data)[1]
   eps <- .Machine$double.eps
   B <- SE <- W <- OK <- WC <- matrix(0,M,N)
   for (j in 1:M)
   {
       for (i in 1:N)
       {
           B[j,i] <- data[paste(prefixb,i,sep="")][j,]
           SE[j,i] <- data[paste(prefixse,i,sep="")][j,]
       }
   }
   K <- BW <- SW <- SWW <- QW <- DL <- P_HETER <- I2 <- rep(0,M)
   BETA_F <- SE_F <- Z_F <- P_F <- BETA_R <- SE_R <- Z_R <- P_R <- rep(0,M)
   for (j in 1:M)
   {
       for (i in 1:N)
       {
           OK[j,i] <- 0
           S1 <- as.numeric(B[j,i])
           S2 <- as.numeric(SE[j,i])
           if (!is.na(S1+S2))
           {
              OK[j,i] <- 1
              T <- max(S2,eps)
              K[j] <- K[j] + 1
              W[j,i] <- 1/T/T
              BW[j] <- BW[j] + S1*W[j,i]
              SW[j] <- SW[j] + W[j,i]
              SWW[j] <- SWW[j]+W[j,i]^2
           }
       }
       if (K[j]>1)
       {
          BETA_F[j] <- BW[j]/SW[j]
          SE_F[j] <- sqrt(1/SW[j])
          Z_F[j] <- BETA_F[j]/SE_F[j]
          P_F[j] <- 2*pnorm(-abs(Z_F[j]))
      }
      for (i in 1:N) if(OK[j,i]==1) QW[j] <- QW[j]+W[j,i]*(B[j,i]-BETA_F[j])^2
      if (K[j]>1)
      {
         P_HETER[j] <- pchisq(QW[j],K[j]-1,lower.tail=FALSE)
         DL[j] <- max(0,(QW[j]-(K[j]-1))/(SW[j]-SWW[j]/SW[j]))
         bwr <- swc <- 0
         for (i in 1:N)
         {
             if (OK[j,i]==1)
             {
                WC[j,i] <- 1/(1/W[j,i]+DL[j])
                bwr <- bwr + B[j,i]*WC[j,i]
                swc <- swc + WC[j,i]
             }
         }
         BETA_R[j] <- bwr/swc
         SE_R[j] <- sqrt(1/swc)
         Z_R[j] <- BETA_R[j]/SE_R[j]
         P_R[j] <- 2*pnorm(-abs(Z_R[j]))
         I2[j] <- (QW[j]-K[j]+1)/QW[j]
         if (I2[j]<0) I2[j] <- 0
      }
      if (toupper(verbose)=="Y")
         cat ("\nMeta-analysis of", N, "studies:\n\n",
           "p_f=", P_F[j], "\n",
           "p_r=", P_R[j], "\n",
           "beta_f=", BETA_F[j], "\n",
           "beta_r=", BETA_R[j], "\n",
           "se_f=",SE_F[j], "\n",
           "se_r=", SE_R[j], "\n",
           "z_f=", Z_F[j], "\n",
           "z_r=", Z_R[j], "\n",
           "p_heter=", P_HETER[j], "\n",
           "i2=", I2[j], "\n",
           "k=", K[j], "\n",
           "eps=",eps, "\n")
   }
   if (toupper(verbose)=="Y")
      cat ("\nwhere\n\n",
           "p_f=P value (fixed effects model)  \n",
           "p_r=P value (random effects model) \n",
           "beta_f=regression coefficient      \n",
           "beta_r=regression coefficient      \n",
           "se_f=standard error                \n",
           "se_r=standard error                \n",
           "z_f=z value                        \n",
           "z_r=z value                        \n",
           "p_heter=heterogeneity test p value \n",
           "i2=I^2                             \n",
           "k=No of tests used                 \n",
           "eps=smallest double-precision number\n")
   invisible(data.frame(beta_f=BETA_F,se_f=SE_F,z_f=Z_F,beta_r=BETA_R,se_r=SE_R,p_heter=P_HETER,i2=I2,k=K,eps=eps))
}
