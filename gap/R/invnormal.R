#' Inverse normal transformation
#'
#' @param x Data with missing values.
#' @export
#' @return
#' Transformed value.
#' @examples
#' x <- 1:10
#' z <- invnormal(x)
#' plot(z,x,type="b")

invnormal <- function(x)
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))) 
