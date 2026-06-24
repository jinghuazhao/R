#' Haplotype Trend Regression
#'
#' @param y Numeric phenotype vector.
#' @param x Haplotype dosage matrix.
#' @param n.sim Number of permutations (default 0).
#'
#' @details
#' Performs haplotype trend regression (HTR) with optional permutation testing.
#' The method follows \insertCite{zaykin02;textual}{gap}.
#'
#' The model fits a linear regression of phenotype on haplotype dosage
#' indicators and computes:
#'
#' - Overall F statistic for joint haplotype effects
#' - Individual haplotype F statistics
#' - Optional permutation-based p-values for both
#'
#' Permutation p-values are computed using
#' \eqn{(b + 1) / (B + 1)} to avoid zero-valued estimates.
#'
#' @export
#' @return
#' A list with the following components:
#'
#' **f**  
#' Overall F statistic for global haplotype association.
#'
#' **p**  
#' P-value for the overall test (asymptotic or permutation-based).
#'
#' **fv**  
#' Vector of F statistics for individual haplotypes.
#'
#' **pi**  
#' Vector of p-values for individual haplotypes.
#'
#' @references
#' \insertRef{zaykin02}{gap}
#'
#' \insertRef{xie05}{gap}
#'
#' @examples
#' \dontrun{
#' ## 26-10-2003
#' ## This example is now part of the demo.
#' \dontrun{
#' if (!requireNamespace("gap.examples", quietly = TRUE))
#'     remotes::install_github("jinghuazhao/R/gap.examples")
#'
#' filespec <- file.path(
#'     find.package("gap.examples"),
#'     "tests", "htr", "test2.dat"
#' )
#' test2 <- read.table(filespec)
#' y <- test2[, 1]
#' x <- test2[, -1]
#' y <- as.matrix(y)
#' x <- as.matrix(x)
#' htr.test2 <- htr(y, x)
#' htr.test2
#' htr.test2 <- htr(y, x, n.sim = 10)
#' htr.test2
#'
#' ## 13-11-2003
#' library(gap.datasets)
#' data(apoeapoc)
#'
#' apoeapoc.gc <- gc.em(apoeapoc[, 5:8])
#' y <- apoeapoc$y
#' y[y == 2] <- 1
#'
#' htr(y, apoeapoc.gc$htrtable)
#'
#' ## 20-08-2008
#' ## Part of the useR! 2008 tutorial by Andrea Foulkes.
#' ## The approach may be used beyond the generalized linear
#' ## model (GLM) framework.
#' HaploEM <- haplo.em(Geno, locus.label = SNPnames)
#' HapMat <- HapDesign(HaploEM)
#'
#' m1 <- lm(Trait ~ HapMat)
#' m2 <- lm(Trait ~ 1)
#'
#' anova(m2, m1)
#' }
#'
htr <- function(y, x, n.sim = 0L)
{
    # ---------- checks ----------
    y <- as.numeric(y)
    x <- as.matrix(x)
    if (length(y) != nrow(x))
        stop("length(y) must equal nrow(x)")
    if (!is.numeric(x))
        stop("x must be numeric")
    if (n.sim < 0L)
        stop("n.sim must be non-negative")
    n.sim <- as.integer(n.sim)
    # ---------- core analysis ----------
    mlr <- function(y, x)
    {
        N <- length(y)
        L <- ncol(x)
        if (L < 2)
            stop("x must contain at least two haplotype columns")
        # HTR parameterization:
        # last haplotype is the reference category
        X <- cbind("(Intercept)" = 1,
                   x[, seq_len(L - 1), drop = FALSE])
        fit <- lm.fit(X, y)
        rss1 <- sum(fit$residuals^2)
        rss0 <- sum((y - mean(y))^2)
        reg.df <- L - 1
        err.df <- N - ncol(X)
        reg.ss <- rss0 - rss1
        reg.ms <- reg.ss / reg.df
        err.ms <- rss1 / err.df
        fstat <- reg.ms / err.ms
        pv.fstat <- pf(fstat,
                       df1 = reg.df,
                       df2 = err.df,
                       lower.tail = FALSE)
        # ---------- single-haplotype tests ----------
        fv <- numeric(L)
        pi <- numeric(L)
        for (i in seq_len(L))
        {
            Xi <- cbind("(Intercept)" = 1,
                        x[, i, drop = FALSE])
            fit.i <- lm.fit(Xi, y)
            rss.i <- sum(fit.i$residuals^2)
            reg.ss.i <- rss0 - rss.i
            reg.df.i <- 1L
            err.df.i <- N - ncol(Xi)
            f.i <- (reg.ss.i / reg.df.i) /
                   (rss.i / err.df.i)
            p.i <- pf(f.i,
                      df1 = reg.df.i,
                      df2 = err.df.i,
                      lower.tail = FALSE)
            fv[i] <- f.i
            pi[i] <- p.i
        }
        list(
            f  = unname(as.numeric(fstat)),
            p  = unname(as.numeric(pv.fstat)),
            fv = unname(as.numeric(fv)),
            pi = unname(as.numeric(pi))
        )
    }
    # ---------- observed statistics ----------
    z0 <- mlr(y, x)
    if (n.sim == 0L)
        return(z0)
    # ---------- permutations ----------
    overall.count <- 0L
    hap.count <- integer(length(z0$fv))
    N <- length(y)
    for (b in seq_len(n.sim))
    {
        yperm <- y[sample.int(N)]
        z <- mlr(yperm, x)
        if (z$f >= z0$f)
            overall.count <- overall.count + 1L
        hap.count <- hap.count + (z$fv >= z0$fv)
    }
    ## Phipson & Smyth (2010):
    ## permutation p-values should never be zero
    p.emp <- (overall.count + 1) / (n.sim + 1)
    pi.emp <- (hap.count + 1) / (n.sim + 1)
    list(
        f  = z0$f,
        p  = p.emp,
        fv = z0$fv,
        pi = pi.emp
    )
}
