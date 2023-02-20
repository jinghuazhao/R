#' Q-Q plot for uniformly distributed random variable
#'
#' @param u a vector of uniformly distributed random variables.
#' @param type string option to specify distribution: "unif"=uniform, "exp"=exponential.
#' @param logscale to use logscale.
#' @param base the base of the log function.
#' @param col color for points.
#' @param lcol color for the diagonal line.
#' @param ci logical option to show confidence interval.
#' @param alpha 1-confidence level, e.g., 0.05.
#' @param ... other options as appropriae for the qqplot function.
#'
#' @details
#' This function produces Q-Q plot for a random variable following uniform distribution with or
#' without using log-scale. Note that the log-scale is by default for type "exp", which is a plot based on
#' exponential order statistics. This appears to be more appropriate than the commonly used procedure whereby
#' the expected value of uniform order statistics is directly log-transformed.
#'
#' @export
#' @return
#' The returned value is a list with components of a qqplot:
#' - x expected value for uniform order statistics or its -log(,base) counterpart.
#' - y observed value or its -log(,base) counterpart.
#'
#' @references
#' \insertRef{balakrishnan03}{gap}
#'
#' \insertRef{casella02}{gap}
#'
#' \insertRef{davison03}{gap}
#'
#' @seealso [`qqfun`]
#'
#' @examples
#' \dontrun{
#' # Q-Q Plot for 1000 U(0,1) r.v., marking those <= 1e-5
#' u_obs <- runif(1000)
#' r <- qqunif(u_obs,pch=21,bg="blue",bty="n")
#' u_exp <- r$y
#' hits <- u_exp >= 2.30103
#' points(r$x[hits],u_exp[hits],pch=21,bg="green")
#' legend("topleft",sprintf("GC.lambda=\%.4f",gc.lambda(u_obs)))
#' }
#'
#' @author Jing Hua Zhao
#' @keywords hplot distribution univar

qqunif <- function(u,type="unif",logscale=TRUE,base=10,col=palette()[4],lcol=palette()[2],ci=FALSE,alpha=0.05,...)
{
  u <- u[!is.na(u)]
  n <- length(u)
  xlabel <- ifelse(logscale,paste("-log",base,"(Expected)",sep=""), "Expected")
  ylabel <- ifelse(logscale,paste("-log",base,"(Observed)",sep=""), "Observed")
  if (ci) c <- abs(qnorm(alpha/2))
  if (type=="exp") {
     n1 <- 1/(n:1)
     n2 <- n1^2
     lambda <- 1/log(base)
     m <- cumsum(n1)*lambda
     z <- qqplot(m,-log(u,base),xlab=xlabel,ylab=ylabel,col=col,...)
     if (ci)
     {
        v <- cumsum(n2)*lambda^2
        s <- sqrt(v)
        lcl <- m - c*s
        ucl <- m + c*s
        lines(m,lcl)
        lines(m,ucl)
     }
  } else if (type=="unif") {
     m <- (1:n)/(n+1)
     if (logscale) z <- qqplot(-log(m,base),-log(u,base),xlab=xlabel,ylab=ylabel,col=col,...)
     else z <- qqplot(m,u,xlab=xlabel,ylab=ylabel,col=col,...)
     if (ci)
     {
        v <- (1:n)*(n-(1:n)+1)/(n+1)^2/(n+2)
        s <- sqrt(v)
        lcl <- m - c*s
        ucl <- m + c*s
        lid <- (lcl>0)
        uid <- (ucl<=1)
        if (logscale)
        {
           a <- -log(m[lid],base)
           b <- -log(lcl[lid],base)
           c <- -log(m[uid],base)
           d <- -log(ucl[uid],base)
        } else {
           a <- m[lid]
           b <- lcl[lid]
           c <- m[uid]
           d <- ucl[uid]
        }
        lines(a,b)
        lines(c,d)
     }
  } else stop ("invalid type")
  abline(0,1,col=lcol)
# polygon(c(m,rev(m)),c(ucl,rev(lcl)),col="gray")
  invisible(z)
}

#09-11-2009 refine U(0,1) and exp(lambda)
#27-08-2008 first attempt of CI
