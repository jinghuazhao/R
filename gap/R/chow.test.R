#' Chow's test for heterogeneity in two regressions
#'
#' @param y1 a vector of dependent variable.
#' @param x1 a matrix of independent variables.
#' @param y2 a vector of dependent variable.
#' @param x2 a matrix of independent variables.
#' @param x a known matrix of independent variables.
#'
#' @details
#' Chow's test is for differences between two or more regressions.  Assuming that
#' errors in regressions 1 and 2 are normally distributed with zero mean and
#' homoscedastic variance, and they are independent of each other, the test of
#' regressions from sample sizes \eqn{n_1} and \eqn{n_2} is then carried out using
#' the following steps.  1.  Run a regression on the combined sample with size
#' \eqn{n=n_1+n_2} and obtain within group sum of squares called \eqn{S_1}.  The
#' number of degrees of freedom is \eqn{n_1+n_2-k}, with \eqn{k} being the number
#' of parameters estimated, including the intercept.  2.  Run two regressions on
#' the two individual samples with sizes \eqn{n_1} and \eqn{n_2}, and obtain their
#' within group sums of square \eqn{S_2+S_3}, with \eqn{n_1+n_2-2k} degrees of
#' freedom.  3.  Conduct an \eqn{F_{(k,n_1+n_2-2k)}} test defined by \deqn{F =
#' \frac{[S_1-(S_2+S_3)]/k}{[(S_2+S_3)/(n_1+n_2-2k)]}} If the \eqn{F} statistic
#' exceeds the critical \eqn{F}, we reject the null hypothesis that the two
#' regressions are equal.
#'
#' In the case of haplotype trend regression, haplotype frequencies from combined
#' data are known, so can be directly used.
#'
#' @export
#' @return The returned value is a vector containing (please use subscript to access them):
#' - F the F statistic.
#' - df1 the numerator degree(s) of freedom.
#' - df2 the denominator degree(s) of freedom.
#' - p the p value for the F test.
#'
#' @references
#' \insertRef{chow60}{gap}
#'
#' @seealso [`htr`]
#'
#' @examples
#' \dontrun{
#' dat1 <- matrix(c(
#'      1.2, 1.9, 0.9,
#'      1.6, 2.7, 1.3,
#'      3.5, 3.7, 2.0,
#'      4.0, 3.1, 1.8,
#'      5.6, 3.5, 2.2,
#'      5.7, 7.5, 3.5,
#'      6.7, 1.2, 1.9,
#'      7.5, 3.7, 2.7,
#'      8.5, 0.6, 2.1,
#'      9.7, 5.1, 3.6), byrow=TRUE, ncol=3)
#'
#' dat2 <- matrix(c(
#'      1.4, 1.3, 0.5,
#'      1.5, 2.3, 1.3,
#'      3.1, 3.2, 2.5,
#'      4.4, 3.6, 1.1,
#'      5.1, 3.1, 2.8,
#'      5.2, 7.3, 3.3,
#'      6.5, 1.5, 1.3,
#'      7.8, 3.2, 2.2,
#'      8.1, 0.1, 2.8,
#'      9.5, 5.6, 3.9), byrow=TRUE, ncol=3)
#'
#' y1<-dat1[,3]
#' y2<-dat2[,3]
#' x1<-dat1[,1:2]
#' x2<-dat2[,1:2]
#' chow.test.r<-chow.test(y1,x1,y2,x2)
#' # from http://aoki2.si.gunma-u.ac.jp/R/
#' }
#'
#' @author Shigenobu Aoki, Jing Hua Zhao
#' @note adapted from chow.R.
#' @keywords htest

chow.test <- function(y1,x1,y2,x2,x=NULL)
{
	mlr <- function(xy)
	{
		N <- nrow(xy)
		P <- ncol(xy)
		R <- cor(xy)
		b <- solve(R[-P,-P], R[,P][-P])
		variance <- var(xy)
		vr <- diag(variance)
		vr <- (vr/vr[P])[-P]
		b <- b/sqrt(vr)
		sse <- (variance[P, P]-(var(xy)[,P][-P])%*%b)*(N-1)
		sse
	}
        xy1<-cbind(x1,y1)
        xy2<-cbind(x2,y2)
	sse12 <- mlr(xy1)+mlr(xy2)
# in case of pooled x is known
        if(!is.null(x))
        {
          xy <- cbind(x,c(y1,y2))
          sse <- mlr(xy)
        }
        else 
	sse <- mlr(rbind(xy1, xy2))
        df1 <- ncol(xy1)
	df2 <- nrow(xy1)+nrow(xy2)-2*(df1)
	f <- (sse-sse12)*df2/(df1*sse12)
	p <- pf(f, df1, df2, lower.tail=FALSE)
	z <- c(f, df1, df2, p)
	names(z) <- c("F value", "d.f.1", "d.f.2", "P value")
	z
}
