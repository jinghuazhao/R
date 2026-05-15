#' Effect-size / Odds-ratio forest plot
#'
#' Create a publication-ready forest plot for model effect estimates.
#' The function supports both linear effect sizes (e.g. regression betas)
#' and exponentiated effects (e.g. odds ratios or hazard ratios).
#'
#' @description
#' The function accepts parameter estimates and their standard errors
#' from one or more models and produces a horizontal forest plot with
#' confidence intervals.
#'
#' Two plotting modes are supported:
#'
#' * `transform="none"` (default): plots estimates on the linear scale  
#'   (typical for GWAS or Mendelian randomisation beta coefficients).
#'
#' * `transform="exp"`: plots exponentiated estimates on a log10 axis  
#'   (typical for odds ratios or hazard ratios). Confidence intervals are
#'   computed on the log scale and back-transformed.
#'
#' @param ESdat Data frame with **three columns**:
#'   \describe{
#'     \item{id}{Model or trait label}
#'     \item{b}{Effect estimate (beta or log(OR)/log(HR))}
#'     \item{se}{Standard error of the estimate}
#'   }
#' @param alpha Type-I error rate for the confidence interval
#'   (default 0.05 for 95% CI).
#' @param fontsize Base font size used in the plot.
#' @param transform Either `"none"` (linear scale) or `"exp"`
#'   (exponentiated scale).
#' @param xlab Optional x-axis label. If `NULL`, a sensible default is used.
#'
#' @details
#' Confidence intervals are computed as
#'
#' \deqn{ estimate \pm z_{\alpha/2} \times SE }
#'
#' When `transform="exp"`, estimates are interpreted as log(OR) or log(HR)
#' and are exponentiated before plotting. The x-axis is displayed on a
#' log10 scale and the reference line is placed at 1.
#'
#' This function replaces an earlier base-R implementation and provides a
#' consistent interface for GWAS, Mendelian randomisation, and
#' epidemiological regression analyses.
#'
#' @return A `ggplot2` plot object.
#'
#' @examples
#' ## Example 1: Linear effect sizes (GWAS / MR)
#' rs12075 <- data.frame(
#'   id=c("CCL2","CCL7","CCL8","CCL11","CCL13","CXCL6","Monocytes"),
#'   b=c(0.1694,-0.0899,-0.0973,0.0749,0.189,0.0816,0.0338387),
#'   se=c(0.0113,0.013,0.0116,0.0114,0.0114,0.0115,0.00713386)
#' )
#' ESplot(rs12075)
#'
#' ## Example 2: Odds ratios
#' dat <- data.frame(
#'   id=c("Basic","Adjusted","Moderate","Heavy","Other"),
#'   b=log(c(4.5,3.5,2.5,1.5,1)),
#'   se=c(0.2,0.1,0.2,0.3,0.2)
#' )
#' ESplot(dat, transform="exp")
#'
#' @author Jing Hua Zhao
#' @keywords hplot
#' @export
#'
ESplot <- function(
  ESdat,
  alpha = 0.05,
  fontsize = 12,
  transform = c("none","exp"),
  xlab = NULL
)
{
  transform <- match.arg(transform)
  id <- ESdat[,1]
  ES <- ESdat[,2]
  SE <- ESdat[,3]
  z <- abs(qnorm(alpha/2))
  lcl <- ES - z*SE
  ucl <- ES + z*SE

  if (transform == "exp") {
    ES  <- exp(ES); lcl <- exp(lcl); ucl <- exp(ucl)
    ref <- 1
    if (is.null(xlab)) xlab <- "Odds ratio"
  } else {
    ref <- 0
    if (is.null(xlab)) xlab <- "Effect size"
  }
  y <- 1:nrow(ESdat)
  ESdata <- data.frame(ES,y,lcl,ucl,id)
  requireNamespace("ggplot2", quietly=TRUE)
  p <- ggplot2::ggplot(ESdata, ggplot2::aes(x=ES, y=y)) +
    ggplot2::geom_point(size=2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin=lcl, xmax=ucl),
      width=0.15, linewidth=0.6, orientation="y"
    ) +
    ggplot2::scale_y_continuous(
      breaks=y, labels=id, name="", trans="reverse"
    ) +
    ggplot2::geom_vline(xintercept=ref, linetype="dashed", alpha=.5) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text=ggplot2::element_text(size=fontsize, colour="black"),
      panel.grid=ggplot2::element_blank(),
      axis.line.x  = ggplot2::element_line(colour="black"),
      axis.ticks.x = ggplot2::element_line(colour="black"),
      panel.spacing=ggplot2::unit(1,"lines")
    )
  if (transform == "exp") {
    ticks <- c(0.25,0.5,1,2,4,8)
    p <- p + ggplot2::scale_x_log10(name=xlab, breaks=ticks, labels=ticks)
  } else {
    p <- p + ggplot2::scale_x_continuous(
      name=xlab, breaks=scales::pretty_breaks(n=6)
    )
  }
  p
}
