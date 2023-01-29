#' Effect-size plot
#'
#' @param ESdat A data frame consisting of model id, parameter estimates and standard errors.
#' @param alpha Type-I error rate used to construct 100(1-alpha) confidence interval.
#' @param fontsize size of font.
#' 
#' @details
#' The function accepts parameter estimates and their standard errors for a range of models.
#'
#' @export
#' @return A high resolution plot object.
#'
#' @examples
#' rs12075 <- data.frame(id=c("CCL2","CCL7","CCL8","CCL11","CCL13","CXCL6","Monocytes"),
#'                       b=c(0.1694,-0.0899,-0.0973,0.0749,0.189,0.0816,0.0338387),
#'                       se=c(0.0113,0.013,0.0116,0.0114,0.0114,0.0115,0.00713386))
#' ESplot(rs12075)
#'
#' # The function replaces an older implementation.
#' within(data.frame(
#'        id=c("Basic model","Adjusted","Moderately adjusted","Heavily adjusted","Other"),
#'        b=log(c(4.5,3.5,2.5,1.5,1)),
#'        se=c(0.2,0.1,0.2,0.3,0.2)
#' ), {
#'    lcl <- exp(b-1.96*se)
#'    ucl <- exp(b+1.96*se)
#'    x <- seq(-2,8,length=length(id))
#'    y <- 1:length(id)
#'    plot(x,y,type="n",xlab="",ylab="",axes=FALSE)
#'    points((lcl+ucl)/2,y,pch=22,bg="black",cex=3)
#'    segments(lcl,y,ucl,y,lwd=3,lty="solid")
#'    axis(1,cex.axis=1.5,lwd=0.5)
#'    par(las=1)
#'    abline(v=1)
#'    axis(2,labels=id,at=y,lty="blank",hadj=0.2,cex.axis=1.5)
#'    title("A fictitious plot")
#' })
#' @author Jing Hua Zhao
#' @keywords hplot

ESplot <- function(ESdat,alpha=0.05,fontsize=12)
{
   id <- ESdat[,1]
   ES <- ESdat[,2]
   SE <- ESdat[,3]
   z <- abs(qnorm(alpha / 2))
   lcl <- ES - z * SE
   ucl <- ES + z * SE
   N <- nrow(ESdat)
   y <- 1:N
   ESdata <- data.frame(ES,y,lcl,ucl,id)
   requireNamespace("ggplot2", quietly=TRUE)
   f <- ggplot2::ggplot(data=ESdata, ggplot2::aes(y, x=ES))+
        ggplot2::geom_point()+
        ggplot2::geom_errorbarh(ggplot2::aes(xmax = ucl, xmin = lcl, height=0.001))+
        ggplot2::scale_x_continuous(name="Effect size")+
        ggplot2::scale_y_continuous(breaks=y,label=with(ESdata,id),name="",trans="reverse")+
        ggplot2::geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)+
        ggplot2::theme_minimal()+
        ggplot2::theme(text=ggplot2::element_text(size=fontsize, color="black"))+
        ggplot2::theme(panel.grid=ggplot2::element_blank())+
        ggplot2::theme(panel.spacing = ggplot2::unit(1, "lines"))
   f
}
