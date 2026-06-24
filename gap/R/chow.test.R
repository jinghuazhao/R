#' Chow Test for Equality of Regression Models
#'
#' @param y1 Numeric vector of responses for the first sample.
#' @param x1 Numeric matrix of predictors for the first sample.
#' @param y2 Numeric vector of responses for the second sample.
#' @param x2 Numeric matrix of predictors for the second sample.
#' @param x Optional pooled design matrix for the combined sample.
#' When supplied, the pooled regression is fitted using `x`
#' instead of `rbind(x1, x2)`.
#'
#' @details
#' Chow's test assesses whether two linear regression models have
#' identical regression coefficients.
#'
#' Let `n1` and `n2` denote the sample sizes of the two groups and
#' let `k` denote the number of estimated parameters, including the
#' intercept.
#'
#' The test proceeds as follows:
#'
#' 1. Fit a regression model to the combined data and obtain the
#'    residual sum of squares `S1`.
#'
#' 2. Fit separate regression models to the two groups and obtain
#'    residual sums of squares `S2` and `S3`.
#'
#' 3. Compute the Chow statistic
#'
#'    \deqn{
#'    F =
#'    \frac{\left[S_1-(S_2+S_3)\right]/k}
#'         {\left[(S_2+S_3)/(n_1+n_2-2k)\right]}
#'    }
#'
#' Under the null hypothesis that the regression coefficients are
#' identical in the two groups, the statistic follows an
#' `F(k, n1 + n2 - 2k)` distribution.
#'
#' In haplotype trend regression applications, pooled haplotype
#' frequencies may be known in advance. In that case, the pooled
#' design matrix can be supplied directly through `x`.
#'
#' @export
#' @return
#' An object of class `"htest"` containing:
#'
#' - `statistic`: Chow F statistic.
#' - `parameter`: numerator and denominator degrees of freedom.
#' - `p.value`: p-value for the test.
#' - `method`: description of the test.
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
#'      9.7, 5.1, 3.6
#' ), byrow = TRUE, ncol = 3)
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
#'      9.5, 5.6, 3.9
#' ), byrow = TRUE, ncol = 3)
#'
#' y1 <- dat1[, 3]
#' y2 <- dat2[, 3]
#' x1 <- dat1[, 1:2]
#' x2 <- dat2[, 1:2]
#'
#' chow.test(y1, x1, y2, x2)
#' }
#'
#' @author
#' Shigenobu Aoki, Jing Hua Zhao
#'
#' @keywords htest regression
#'
chow.test <- function(y1, x1, y2, x2, x = NULL)
{
    rss <- function(y, X)
    {
        fit <- lm.fit(cbind(1, as.matrix(X)), y)
        sum(fit$residuals^2)
    }
    k <- ncol(x1) + 1L
    rss.sep <-
        rss(y1, x1) +
        rss(y2, x2)
    rss.pool <-
        if (is.null(x))
            rss(c(y1, y2), rbind(x1, x2))
        else
            rss(c(y1, y2), x)
    df1 <- k
    df2 <- length(y1) + length(y2) - 2 * k
    F <- ((rss.pool - rss.sep) / df1) /
         (rss.sep / df2)
    p <- pf(F, df1, df2, lower.tail = FALSE)
    structure(
        list(
            statistic = c(F = F),
            parameter = c(df1 = df1, df2 = df2),
            p.value = p,
            method = "Chow test for equality of regression coefficients"
        ),
        class = "htest"
    )
}
