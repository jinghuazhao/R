#' Mendelian Randomization forest plot
#'
#' @param dat A data.frame with outcome id, effect size and standard error.
#' @param sm Summary measure such as OR, RR, MD.
#' @param title Title of the meta-analysis.
#' @param ... Other options for meta::forest().
#'
#' @details
#' This is a wrapper of meta::forest() for multi-outcome Mendelian Randomization. It allows for the flexibility of both binary and continuous outcomes with and without summary level statistics.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  tnfb <- '
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
#'          "IGA glomerulonephritis" -0.32696600 0.105262000
#'                   "atopic eczema" -0.00204018 0.000678061
#'  '
#'  tnfb <- as.data.frame(scan(file=textConnection(tnfb),what=list("",0,0)))
#'  names(tnfb) <-c("outcome","Effect","StdErr")
#'
#'  # default output
#'  mr_forestplot(tnfb[c("outcome","Effect","StdErr")],colgap.forest.left = "0.5cm",fontsize=14,leftlabs=c("Outcome","b","SE"), spacing=1.6)
#'  # no summary level statistics
#'  mr_forestplot(tnfb[c("outcome","Effect","StdErr")], colgap.forest.left = "0.5cm", fontsize=14, leftlabs=c("Outcome","b","SE"), rightlabs="ci",
#'                plotwidth="5inch",sm="OR",
#'                common=FALSE, random=FALSE, print.I2=FALSE, print.pval.Q=FALSE, print.tau2=FALSE, addrow=TRUE, backtransf=TRUE, spacing=1.6)
#' }

mr_forestplot <- function(dat,sm="MD",title="",...)
{
   requireNamespace("meta")
   outcome <- dat[,1]
   Effect <- dat[,2]
   StdErr <- dat[,3]
   mg <- meta::metagen(Effect,StdErr,sprintf("%s",outcome),sm=sm,title=title)
   meta::forest(mg,...)
   with(mg,cat("Q =", Q, "df =", df.Q, "p =", pval.Q, "I2 =", I2, "lower.I2 =", lower.I2, "upper.I2 =", upper.I2, "\n"))
}
