ESplot <- function(ESdat,SE=TRUE,logscale=TRUE,alpha=0.05,xlim=c(-2,8),v=1,...)
{
   n <- dim(ESdat)[1]
   models <- ESdat[,1]
   ES <- ESdat[,2]
   if (SE)
   {
      z <- abs(qnorm(alpha / 2))
   # effect size
      SE <- ESdat[,3]
      LCL <- ES - z * SE
      UCL <- ES + z * SE
      if (logscale)
      {
      # log-scale, i.e., OR
         LCL <- exp(LCL)
         UCL <- exp(UCL)
      }
   } else {
   # draw directly
      LCL <- ESdat[,3]
      UCL <- ESdat[,4]
   }
   x <- seq(xlim[1],xlim[2],length=n)
   y <- 1:n
   plot(x,y,type="n",xlab="",ylab="",axes=FALSE,...)
   points((LCL+UCL)/2,y,pch=22,bg="black",cex=3,...)
   segments(LCL,y,UCL,y,lwd=3,lty="solid")
   axis(1,cex.axis=1.5,lwd=0)
   par(las=1)
   abline(v=v,...)
   axis(2,labels=models,at=y,lty="blank",hadj=0.2,cex.axis=1.5)
}
