#' Mendelian Randomization forest plot
#'
#' @param dat A data.frame with outcome id, effect size and standard error.
#' @param sm  Summary measure such as OR, RR, MD.
#' @param title Title of the meta-analysis.
#' @param ... Additional arguments passed to meta::forest().
#'
#' @details
#' Wrapper around meta::forest() for multi-outcome Mendelian Randomization.
#' Works for binary and continuous outcomes, with or without summary statistics.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## Example data -----------------------------------------------------------
#' tnfb <- '
#'              "multiple sclerosis"  0.69058600 0.059270400
#'    "systemic lupus erythematosus"  0.76687500 0.079000500
#'          "sclerosing cholangitis"  0.62671500 0.075954700
#'   "juvenile idiopathic arthritis" -1.17577000 0.160293000
#'                       "psoriasis"  0.00582586 0.000800016
#'            "rheumatoid arthritis" -0.00378072 0.000625160
#'      "inflammatory bowel disease" -0.14334200 0.025272500
#'          "ankylosing spondylitis" -0.00316852 0.000626225
#'                  "hypothyroidism" -0.00432054 0.000987324
#'               "allergic rhinitis"  0.00393075 0.000926002
#'          "IgA glomerulonephritis" -0.32696600 0.105262000
#'                   "atopic eczema" -0.00204018 0.000678061
#' '
#'
#' tnfb <- as.data.frame(scan(file = textConnection(tnfb), what = list("",0,0)))
#' names(tnfb) <- c("outcome","Effect","StdErr")
#' tnfb$outcome <- gsub("\\b(^[a-z])","\\U\\1", tnfb$outcome, perl = TRUE)
#'
#'
#' ## 1) Default meta-style forest plot (b, SE, CI + weights) ----------------
#' mr_forestplot(
#'   tnfb,
#'   colgap.forest.left = "0.05cm",
#'   fontsize = 14,
#'   leftcols  = c("studlab","effect","seTE","ci"),
#'   leftlabs  = c("Outcome","b","SE","95% CI"),
#'   rightcols = c("w.common","w.random"),
#'   rightlabs = c("Weight (FE)","Weight (RE)"),
#'   common = FALSE, random = FALSE,
#'   print.I2 = FALSE, print.pval.Q = FALSE, print.tau2 = FALSE,
#'   spacing = 1.6, digits.TE = 2, digits.seTE = 2
#' )
#'
#'
#' ## 2) MR summary (OR + CI only) -------------------------------------------
#' mr_forestplot(
#'   tnfb,
#'   sm = "OR",
#'   backtransf = TRUE,
#'   colgap.forest.left = "0.05cm",
#'   fontsize = 14,
#'   leftcols  = "studlab",
#'   leftlabs  = "Outcome",
#'   rightcols = c("effect","ci"),
#'   rightlabs = c("OR","95% CI"),
#'   sortvar = tnfb$Effect,
#'   common = FALSE, random = FALSE,
#'   print.I2 = FALSE, print.pval.Q = FALSE, print.tau2 = FALSE,
#'   spacing = 1.6
#' )
#'
#'
#' ## 3) MR summary with P-values --------------------------------------------
#' mr_forestplot(
#'   tnfb,
#'   sm = "OR",
#'   backtransf = TRUE,
#'   colgap.forest.left = "0.05cm",
#'   fontsize = 14,
#'   leftcols  = "studlab",
#'   leftlabs  = "Outcome",
#'   rightcols = c("effect","ci","pval"),
#'   rightlabs = c("OR","95% CI","P"),
#'   digits = 3, digits.pval = 2, scientific.pval = TRUE,
#'   sortvar = tnfb$Effect,
#'   common = FALSE, random = FALSE,
#'   print.I2 = FALSE, print.pval.Q = FALSE, print.tau2 = FALSE,
#'   addrow = TRUE,
#'   spacing = 1.6
#' )
#' }

mr_forestplot <- function (dat, sm = "", title = "", ...)
{
    requireNamespace("meta")
    outcome <- dat[, 1]
    Effect <- dat[, 2]
    StdErr <- dat[, 3]
    mg <- meta::metagen(Effect, StdErr, sprintf("%s", outcome),
        sm = sm, title = title, method.tau.ci = "")
    meta::forest(mg, ...)
    with(mg, cat("Q =", Q, "df =", df.Q, "p =", pval.Q, "I2 =",
        I2, "lower.I2 =", lower.I2, "upper.I2 =", upper.I2, "\n"))
}
