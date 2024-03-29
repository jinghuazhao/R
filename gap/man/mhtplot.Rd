% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mhtplot.R
\name{mhtplot}
\alias{mhtplot}
\title{Manhattan plot}
\usage{
mhtplot(data, control = mht.control(), hcontrol = hmht.control(), ...)
}
\arguments{
\item{data}{a data frame with three columns representing chromosome, position and p values.}

\item{control}{A control function named mht.control() with the following arguments:
\itemize{
\item type a flag with value "p" or "l" indicating if points or lines are to be drawn.
\item usepos a flag to use real chromosomal positions as composed to ordinal positions with default value FALSE.
\item logscale a flag to indicate if p value is to be log-transformed with default value TRUE.
\item base the base of the logarithm with default value 10.
\item cutoffs the cut-offs where horizontal line(s) are drawn with default value NULL.
\item colors the color for different chromosome(s), and random if unspecified with default values NULL.
\item labels labels for the ticks on x-axis with default value NULL.
\item srt degree to which labels are rotated with default value of 45.
\item gap gap between chromosomes with default value NULL.
\item cex cex for the data points.
\item yline Margin line position.
\item xline Margin line position.
}}

\item{hcontrol}{A control function named hmht.control() with the following arguments:
\itemize{
\item data. chunk of data to be highlighted with default value NULL.
\item colors. colors for annotated genes.
\item yoffset. offset above the data point showing most significant p value with default value 0.5.
\item cex shrinkage factor for data points with default value 1.5.
\item boxed if the label for the highlited region with default value FALSE.
}}

\item{...}{other options in compatible with the R plot function.}
}
\value{
The plot is shown on or saved to the appropriate device.
}
\description{
Manhattan plot
}
\details{
To generate Manhattan plot, e.g., of genomewide significance (p values) and
a random variable that is uniformly distributed. By default, a log10-transformation is applied.
Note that with real chromosomal positions, it is also appropriate to plot and some but not all chromosomes.

It is possible to specify options such as xlab and ylim when the plot is requested for data in other context.
}
\examples{
\dontrun{
# foo example
test <- matrix(c(1,1,4,1,1,6,1,10,3,2,1,5,2,2,6,2,4,8),byrow=TRUE,6)
mhtplot(test)
mhtplot(test,mht.control(logscale=FALSE))

# fake example with Affy500k data
affy <-c(40220, 41400, 33801, 32334, 32056, 31470, 25835, 27457, 22864, 28501, 26273,
         24954, 19188, 15721, 14356, 15309, 11281, 14881, 6399, 12400, 7125, 6207)
CM <- cumsum(affy)
n.markers <- sum(affy)
n.chr <- length(affy)
test <- data.frame(chr=rep(1:n.chr,affy),pos=1:n.markers,p=runif(n.markers))

# to reduce size of the plot
# bitmap("mhtplot.bmp",res=72*5)
oldpar <- par()
par(cex=0.6)
colors <- rep(c("blue","green"),11)
# other colors, e.g.
# colors <- c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green",
#             "cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray",
#             "magenta","red")
mhtplot(test,control=mht.control(colors=colors),pch=19,srt=0)
title("A simulated example according to EPIC-Norfolk QCed SNPs")
axis(2)
axis(1,pos=0,labels=FALSE,tick=FALSE)
abline(0,0)
# dev.off()
par(oldpar)

mhtplot(test,control=mht.control(usepos=TRUE,colors=colors,gap=10000),pch=19,bg=colors)
title("Real positions with a gap of 10000 bp between chromosomes")
box()

png("manhattan.png",height=3600,width=6000,res=600)
opar <- par()
par(cex=0.4)
ops <- mht.control(colors=rep(c("lightgray","lightblue"),11),srt=0,yline=2.5,xline=2)
require(gap.datasets)
mhtplot(mhtdata[,c("chr","pos","p")],ops,xlab="",ylab="",srt=0)
axis(2,at=1:16)
title("An adaptable plot as .png")
par(opar)
dev.off()

data <- with(mhtdata,cbind(chr,pos,p))
glist <- c("IRS1","SPRY2","FTO","GRIK3","SNED1","HTR1A","MARCH3","WISP3","PPP1R3B",
         "RP1L1","FDFT1","SLC39A14","GFRA1","MC4R")
hdata <- subset(mhtdata,gene\\%in\\%glist)[c("chr","pos","p","gene")]
color <- rep(c("lightgray","gray"),11)
glen <- length(glist)
hcolor <- rep("red",glen)
par(las=2, xpd=TRUE, cex.axis=1.8, cex=0.4)
ops <- mht.control(colors=color,yline=1.5,xline=3,labels=paste("chr",1:22,sep=""),
                   srt=270)
hops <- hmht.control(data=hdata,colors=hcolor)
mhtplot(data,ops,hops,pch=19)
axis(2,pos=2,at=1:16)
title("Manhattan plot with genes highlighted",cex.main=1.8)

mhtplot(data,mht.control(cutoffs=c(4,6,8,16)),pch=19)
title("Another plain Manhattan plot")

# Miami plot

test <- within(test, {pr=1-p})
miamiplot(test,chr="chr",bp="pos",p="p",pr="pr")
}
}
\seealso{
\code{\link{qqunif}}
}
\author{
Jing Hua Zhao
}
\keyword{hplot}
