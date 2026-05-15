#' Mendelian Randomization wrapper (IVW, Egger, Weighted Median, Penalised WM)
#'
#' @param data Data frame containing SNP summary statistics in wide format.
#' @param X Character string. Exposure trait name.
#' @param Y Character vector. Outcome trait names.
#' @param alpha Significance level used for confidence intervals (default 0.05).
#' @param other_plots Logical. If TRUE, produces metafor funnel and forest plots.
#'
#' Performs two–sample Mendelian Randomization using summary statistics
#' stored in a wide data frame with columns named:
#'
#' \itemize{
#'   \item SNP
#'   \item b.trait   : effect estimate
#'   \item SE.trait  : standard error
#' }
#'
#' For each outcome in `Y`, the function:
#' \enumerate{
#'   \item Extracts instruments for exposure `X`
#'   \item Computes IVW, MR-Egger, Weighted Median and Penalised Weighted Median estimates
#'   \item Computes Cochran’s Q heterogeneity statistic
#'   \item Produces MR scatter plots
#' }
#'
#' @return Invisible list with components:
#' \describe{
#'   \item{r}{Matrix of MR estimates for each exposure–outcome pair.}
#'   \item{plots}{List of ggplot objects (MR scatter plots).}
#' }
#'
#' @details
#' Methods implemented:
#' \itemize{
#'   \item IVW (inverse-variance weighted regression)
#'   \item MR-Egger regression
#'   \item Weighted median estimator
#'   \item Penalised weighted median
#'   \item Cochran’s Q heterogeneity statistic (metafor)
#' }
#'
#' Required packages:
#' \code{metafor}, \code{ggplot2}, \code{cowplot}.
#' Functions \code{weighted.median()} and \code{mr.boot()} must be available
#' in the environment (e.g. from Mendelian randomization toolkits).
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("metafor", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("cowplot", quietly = TRUE)) {
#' txt <- "
#' rs188743906  0.6804 0.1104  0.00177 0.01660        NA        NA
#' rs2289779   -0.0788 0.0134  0.00104 0.00261 -0.007543 0.0092258
#' rs117804300 -0.2281 0.0390 -0.00392 0.00855  0.109372 0.0362219
#' rs7033492   -0.0968 0.0147 -0.00585 0.00269  0.022793 0.0119903
#' rs10793962   0.2098 0.0212  0.00378 0.00536 -0.014567 0.0138196
#' rs635634    -0.2885 0.0153 -0.02040 0.00334  0.077157 0.0117123
#' rs176690    -0.0973 0.0142  0.00293 0.00306 -0.000007 0.0107781
#' rs147278971 -0.2336 0.0378 -0.01240 0.00792  0.079873 0.0397491
#' rs11562629   0.1155 0.0181  0.00960 0.00378 -0.010040 0.0151460
#' "
#' v <- c("SNP","b.LIF.R","SE.LIF.R","b.FEV1","SE.FEV1","b.CAD","SE.CAD")
#' mrdat <- setNames(as.data.frame(scan(text = txt,
#'                   what = list("",0,0,0,0,0,0), quiet=TRUE)), v)
#' res <- mr(mrdat, "LIF.R", c("CAD","FEV1"), other_plots=TRUE)
#' r <- as.data.frame(res$r, stringsAsFactors=FALSE)
#' rownames(r) <- r$IV
#' r$IV <- NULL
#' r[] <- lapply(r, function(x) as.numeric(as.character(x)))
#' b_cols  <- grep("^b[A-Z]", names(r), value=TRUE)
#' se_cols <- paste0("se", substring(b_cols, 2))
#' keep <- se_cols %in% names(r)
#' b_cols <- b_cols[keep]; se_cols <- se_cols[keep]
#' if(length(b_cols) > 0){
#'   pvals <- sapply(seq_along(b_cols), function(i)
#'     2 * pnorm(-abs(r[[b_cols[i]]] / r[[se_cols[i]]])))
#'   pvals <- as.data.frame(pvals)
#'   colnames(pvals) <- substring(b_cols, 2)
#'   pvals[] <- lapply(pvals, format.pval, digits=3, eps=1e-4)
#'   results_table <- cbind(round(r,3), pvals)
#' } else {
#'   results_table <- round(r,3)
#' }
#' results_table
#' res$plots
#' }
#' }
#'
#' @export
mr <- function (data, X, Y, alpha = 0.05, other_plots = FALSE)
{
    cval <- qnorm(alpha/2, lower.tail = FALSE)
    nxy <- character(length(Y))
    results_list <- list()
    ES_plot <- vector("list", length(Y))
    k <- 1
    for (i in Y) {
        dat <- data.frame(
            SNP  = data[["SNP"]],
            bzx  = data[[paste0("b.", X)]],
            SEzx = data[[paste0("SE.", X)]],
            bzy  = data[[paste0("b.", i)]],
            SEzy = data[[paste0("SE.", i)]]
        )
        dat <- dat[complete.cases(dat) & dat$bzx != 0, ]
        SNP  <- dat$SNP
        bzx  <- dat$bzx
        SEzx <- dat$SEzx
        bzy  <- dat$bzy
        SEzy <- dat$SEzy
        nxy[k] <- paste0(X, ".", i)
        m1 <- lm(bzy ~ bzx - 1, weights = SEzy^-2)
        coef.summary.m1 <- coef(summary(m1))
        sebIVW_reg <- coef.summary.m1[1,2] / min(summary(m1)$sigma, 1)
        m <- lm(bzy ~ bzx, weights = SEzy^-2)
        coef.summary.m <- coef(summary(m))
        bEGGER   <- coef.summary.m[2,1]
        sebEGGER <- coef.summary.m[2,2] / min(summary(m)$sigma, 1)
        intEGGER <- coef.summary.m[1,1]
        seintEGGER <- coef.summary.m[1,2]
        bIV <- bzy / bzx
        weights <- (SEzy / bzx)^-2
        bIVW <- sum(bzy * bzx * SEzy^-2) / sum(bzx^2 * SEzy^-2)
        bWM  <- weighted.median(bIV, weights)
        sebWM <- mr.boot(bzx, SEzx, bzy, SEzy, weights)
        penalty <- pchisq(weights * (bIV - bIVW)^2, df = 1, lower.tail = FALSE)
        penalty.weights <- weights * pmin(1, penalty * 20)
        bPWM  <- weighted.median(bIV, penalty.weights)
        sebPWM <- mr.boot(bzx, SEzx, bzy, SEzy, penalty.weights)
        res <- metafor::rma(bIV, 1/sqrt(weights), method = "FE")
        CochQ  <- res$QE
        CochQp <- res$QEp
        graph_title <- paste("SNP effect on", i)
        x_title <- paste("Effect on", X)
        y_title <- paste("Effect on", i)
        if (other_plots) {
            metafor::funnel(res, main = graph_title, xlab = y_title)
            abline(v = 0, lty = "dashed", col = "red")
            metafor::forest(res, main = graph_title, addfit = FALSE,
                            slab = SNP, xlab = y_title)
        }
        bzxLCL <- bzx - cval * SEzx
        bzxUCL <- bzx + cval * SEzx
        bzyLCL <- bzy - cval * SEzy
        bzyUCL <- bzy + cval * SEzy
        group4 <- data.frame(
            intercept = c(0, intEGGER, 0, 0),
            slope = c(bIVW, bEGGER, bWM, bPWM),
            colour = c("red","blue","green","orange"),
            method = c("IVW","Egger","WM","PWM")
        )
        ES_plot[[k]] <- ggplot2::ggplot(
            data.frame(bzx,SEzx,bzy,SEzy,bzxLCL,bzxUCL,bzyLCL,bzyUCL),
            ggplot2::aes(x=bzx,y=bzy)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw() +
            cowplot::theme_cowplot(12) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin=bzyLCL,ymax=bzyUCL)) +
            ggplot2::geom_errorbarh(ggplot2::aes(xmin=bzxLCL,xmax=bzxUCL)) +
            ggplot2::geom_abline(data=group4,
                                 ggplot2::aes(intercept=intercept,slope=slope,colour=method),
                                 size=1) +
            ggplot2::scale_colour_manual(values=group4$colour) +
            ggplot2::theme(legend.position="bottom") +
            ggplot2::ggtitle(graph_title) +
            ggplot2::xlab(x_title) +
            ggplot2::ylab(y_title)
        plot(ES_plot[[k]])
        results_list[[k]] <- c(
            nxy[k], bIVW, sebIVW_reg, CochQ, CochQp,
            bEGGER, sebEGGER, intEGGER, seintEGGER,
            bWM, sebWM, bPWM, sebPWM
        )
        k <- k + 1
    }
    r <- do.call(rbind, results_list)
    colnames(r) <- c("IV","bIVW","sebIVW","CochQ","CochQp","bEGGER",
                     "sebEGGER","intEGGER","seintEGGER","bWM",
                     "sebWM","bPWM","sebPWM")
    invisible(list(r = r, plots = ES_plot))
}

#' A function for obtaining weighted median with interpolation
#' @author Adapted code by Felix Day 16/9/2015
#' @references
#' \insertRef{bowden15}{gap}
#' @noRd
#'
weighted.median <- function(x, w)
{
   x.order <- x[order(x)]
   w.order <- w[order(x)]
   ws <- cumsum(w.order)-0.5*w.order
   ws <- ws / sum(w.order)
   below <- max(which(ws < 0.5))
   x.order[below] + (x.order[below+1]-x.order[below]) * (0.5-ws[below]) / (ws[below+1]-ws[below])
}

mr.boot = function(bXG, sebXG, bYG, sebYG, w, n.boot=1000, method="median")
{
   if (interactive()) n.boot <- 1000 else n.boot <- 100
   boot <- vector('numeric')
   for(i in 1:n.boot)
   {
       bXG.boot <- rnorm(length(bXG), mean = bXG, sd = sebXG)
       bYG.boot <- rnorm(length(bYG), mean = bYG, sd = sebYG)
       if(method=="Egger")
       {
         bYG.boot = bYG.boot*sign(bXG.boot)
         bXG.boot = abs(bXG.boot)
         boot[i] = coef(summary(lm(bYG.boot~bXG.boot,weights=sebYG^-2)))[2,1]
       } else {
         bIV.boot <- bYG.boot / bXG.boot
         boot[i] <- weighted.median(bIV.boot, w)
      }
   }
   sd(boot)
}
