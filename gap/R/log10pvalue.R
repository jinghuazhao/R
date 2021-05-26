#' log10(p) for a P value including its scientific format
#'
#' @md
#' @param p value.
#' @param base base part in scientific format.
#' @param exponent exponent part in scientific format.
#' @export
#' @return log10(P)
#' @examples
#' log10pvalue(1e-323)
#' log10pvalue(base=1,exponent=-323)

log10pvalue <- function(p=NULL,base=NULL,exponent=NULL)
{
  if(!is.null(p))
  {
    p <- format(p,scientific=TRUE)
    p2 <- strsplit(p,"e")
    base <- as.numeric(lapply(p2,"[",1))
    exponent <- as.numeric(lapply(p2,"[",2))
  } else if(is.null(base) | is.null(exponent)) stop("base and exponent should both be specified")
  log10(base)+exponent
}
