#' Manhattan plot
#'
#' @param data a data frame with three columns representing chromosome, position and p values.
#' @param control A control function named mht.control() with the following arguments:
#' - type a flag with value "p" or "l" indicating if points or lines are to be drawn.
#' - usepos a flag to use real chromosomal positions as composed to ordinal positions with default value FALSE.
#' - logscale a flag to indicate if p value is to be log-transformed with default value TRUE.
#' - base the base of the logarithm with default value 10.
#' - cutoffs the cut-offs where horizontal line(s) are drawn with default value NULL.
#' - colors the color for different chromosome(s), and random if unspecified with default values NULL.
#' - labels labels for the ticks on x-axis with default value NULL.
#' - srt degree to which labels are rotated with default value of 45.
#' - gap gap between chromosomes with default value NULL.
#' - cex cex for the data points.
#' - yline Margin line position.
#' - xline Margin line position.
#' @param hcontrol A control function named hmht.control() with the following arguments:
#' - data. chunk of data to be highlighted with default value NULL.
#' - colors. colors for annotated genes.
#' - yoffset. offset above the data point showing most significant p value with default value 0.5.
#' - cex shrinkage factor for data points with default value 1.5.
#' - boxed if the label for the highlited region with default value FALSE.
#' @param ... other options in compatible with the R plot function.
#'
#' @details
#' To generate Manhattan plot, e.g., of genomewide significance (p values) and
#' a random variable that is uniformly distributed. By default, a log10-transformation is applied.
#' Note that with real chromosomal positions, it is also appropriate to plot and some but not all chromosomes.
#'
#' It is possible to specify options such as xlab and ylim when the plot is requested for data in other context.
#'
#' @export
#' @return
#' The plot is shown on or saved to the appropriate device.
#' @seealso [`qqunif`]
#' @examples
#' \dontrun{
#' # foo example
#' test <- matrix(c(1,1,4,1,1,6,1,10,3,2,1,5,2,2,6,2,4,8),byrow=TRUE,6)
#' mhtplot(test)
#' mhtplot(test,mht.control(logscale=FALSE))
#'
#' # fake example with Affy500k data
#' affy <-c(40220, 41400, 33801, 32334, 32056, 31470, 25835, 27457, 22864, 28501, 26273,
#'          24954, 19188, 15721, 14356, 15309, 11281, 14881, 6399, 12400, 7125, 6207)
#' CM <- cumsum(affy)
#' n.markers <- sum(affy)
#' n.chr <- length(affy)
#' test <- data.frame(chr=rep(1:n.chr,affy),pos=1:n.markers,p=runif(n.markers))
#'
#' # to reduce size of the plot
#' # bitmap("mhtplot.bmp",res=72*5)
#' oldpar <- par()
#' par(cex=0.6)
#' colors <- rep(c("blue","green"),11)
#' # other colors, e.g.
#' # colors <- c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green",
#' #             "cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray",
#' #             "magenta","red")
#' mhtplot(test,control=mht.control(colors=colors),pch=19,srt=0)
#' title("A simulated example according to EPIC-Norfolk QCed SNPs")
#' axis(2)
#' axis(1,pos=0,labels=FALSE,tick=FALSE)
#' abline(0,0)
#' # dev.off()
#' par(oldpar)
#'
#' mhtplot(test,control=mht.control(usepos=TRUE,colors=colors,gap=10000),pch=19,bg=colors)
#' title("Real positions with a gap of 10000 bp between chromosomes")
#' box()
#'
#' png("manhattan.png",height=3600,width=6000,res=600)
#' opar <- par()
#' par(cex=0.4)
#' ops <- mht.control(colors=rep(c("lightgray","lightblue"),11),srt=0,yline=2.5,xline=2)
#' require(gap.datasets)
#' mhtplot(mhtdata[,c("chr","pos","p")],ops,xlab="",ylab="",srt=0)
#' axis(2,at=1:16)
#' title("An adaptable plot as .png")
#' par(opar)
#' dev.off()
#'
#' data <- with(mhtdata,cbind(chr,pos,p))
#' glist <- c("IRS1","SPRY2","FTO","GRIK3","SNED1","HTR1A","MARCH3","WISP3","PPP1R3B",
#'          "RP1L1","FDFT1","SLC39A14","GFRA1","MC4R")
#' hdata <- subset(mhtdata,gene\%in\%glist)[c("chr","pos","p","gene")]
#' color <- rep(c("lightgray","gray"),11)
#' glen <- length(glist)
#' hcolor <- rep("red",glen)
#' par(las=2, xpd=TRUE, cex.axis=1.8, cex=0.4)
#' ops <- mht.control(colors=color,yline=1.5,xline=3,labels=paste("chr",1:22,sep=""),
#'                    srt=270)
#' hops <- hmht.control(data=hdata,colors=hcolor)
#' mhtplot(data,ops,hops,pch=19)
#' axis(2,pos=2,at=1:16)
#' title("Manhattan plot with genes highlighted",cex.main=1.8)
#'
#' mhtplot(data,mht.control(cutoffs=c(4,6,8,16)),pch=19)
#' title("Another plain Manhattan plot")
#'
#' # Miami plot
#'
#' test <- within(test, {pr=1-p})
#' miamiplot(test,chr="chr",bp="pos",p="p",pr="pr")
#'}
#' @author Jing Hua Zhao
#' @keywords hplot

mhtplot <- function(data, control=mht.control(), hcontrol=hmht.control(), ...) {
  for(p in c("grid")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
        warning(paste("mhtplot needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  data2 <- data[!apply(is.na(data),1,any),]
  n2 <- dim(data2[1])
  chr <- data2[,1]
  pos <- newpos <- data2[,2]
  p <- data2[,3]
  tablechr <- table(chr)
  allchr <- as.vector(tablechr)
  n.chr <- length(allchr)
  type <- control$type
  usepos <- control$usepos
  logscale <- control$logscale
  base <- control$base
  cutoffs <- control$cutoffs
  colors <- control$colors
  labels <- control$labels
  srt <- control$srt
  gap <- control$gap
  pcex <- control$cex
  yline <- control$yline
  xline <- control$xline
  colorlist <- colors()
  if(is.null(colors)) colors <- sample(colorlist,n.chr)
  tablechr <- unique(chr)
  if(is.null(labels)) labels <- tablechr
  if(is.null(gap)) gap <- 0
  if (!is.null(hcontrol$data))
  {
     hdata <- hcontrol$data
     hdata2 <- hdata[!apply(is.na(hdata),1,any),]
     hchr <- hdata2[,1]
     hpos <- hnewpos <- hdata2[,2]
     hp <- hdata2[,3]
     hname <- hdata2[,4]
     hcolors <- hcontrol$colors
     hyoffs <- hcontrol$yoffs
     hboxed <- hcontrol$boxed
     hcex <- hcontrol$cex
     htablechr <- unique(hchr)
     hn.chr <- length(htablechr)
     hlabels <- unique(hname)
     htablechrname <- unique(data.frame(hchr, hname))
     if(is.null(hcolors)) hcolors <- rep("red",length(hlabels))
     else hcolors <- hcontrol$colors
  }
  CMindex <- cumsum(allchr)
  for(i in 1:n.chr)
  {
     u <- CMindex[i]
     l <- CMindex[i]-allchr[i]+1
     chr <- l:u
     if (usepos) d <- diff(pos[chr]) else d <- rep(1,allchr[i]-1)
     newpos[chr] <- c(gap,d)
  }
  CM <- cumsum(as.numeric(newpos))
  args <- list(...)
  if ("ylim"%in%names(args)) dp <- seq(args$ylim[1],args$ylim[2],length=sum(allchr))
  else dp <- seq(min(p),max(p),length=sum(allchr))
  if (logscale) y <- -log(dp,base) else y <- dp
  y1 <- min(y)
  par(xaxt="n",yaxt="n")
  xy <- xy.coords(CM,y)
  plot(xy$x,xy$y,type="n",ann=FALSE,axes=FALSE,...)
  axis(1)
  axis(2)
  par(xaxt="s",yaxt="s")
  for(i in 1:n.chr)
  {
     u <- CMindex[i]
     l <- CMindex[i]-allchr[i]+1
     chr <- l:u
     cat("Plotting points ",l,"-",u,"\n");
     if (logscale) y <- -log(p[chr],base) else y <- p[chr]
     col.chr <- colors[i]
     if(type=="l") lines(CM[chr],y,col=col.chr,cex=pcex,...)
     else points(CM[chr],y,col=col.chr,cex=pcex,...)
     text(ifelse(i==1,CM[1],CM[l]),y1,pos=1,offset=1.5,labels[i],srt=srt,cex=0.5,...)
  }
  j <- 1
  for(i in 1:n.chr)
  {
     u <- CMindex[i]
     l <- CMindex[i]-allchr[i]+1
     chr <- l:u
     if (logscale) y <- -log(p[chr],base) else y <- p[chr]
     col.chr <- colors[i]
     if (!is.null(hcontrol$data))
     {
        chrs <- htablechrname[tablechr[i]==htablechrname[,1],]
        if(dim(chrs)[1]>0) {
          hchrs <- as.character(chrs[,2])
          for(k in 1:length(hchrs))
          {
             hregion <- hpos[hchr==chrs[k,1]&hname==hchrs[k]]
             hl <- chr[pos[chr]==hregion[1]]
             hu <- chr[pos[chr]==hregion[length(hregion)]]
             cat("  ... highlighting",hl,"-",hu,hchrs[k],"\n")
             l1 <- hl-l+1
             l2 <- hu-l+1
             col.chr[l1:l2] <- hcolors[j]
             if (hboxed) textbox(hchrs[k], name="tbt", vp=grid::viewport(x = CM[chr][l1]/max(CM), y = (max(y[l1:l2]) +  hyoffs)/max(y)))
             else text(CM[chr][l1],max(y[l1:l2]+hyoffs),hchrs[k],cex=1,font=3)
             points(CM[l+(l1:l2)],y[l1:l2],col=col.chr[l1:l2],cex=pcex,...)
             j <- j+1
          }
        }
     }
  }
  if(!is.null(cutoffs)) segments(0,cutoffs,n2+gap*n.chr,cutoffs) # abline(h=cutoffs)
  if ("ylab"%in%names(args)) mtext(args$ylab,2,line=yline,las=0) else
     mtext(ifelse(logscale,paste("-log",base,"(Observed value)",sep=""),"Observed value"),2,line=yline,las=0,cex=0.5)
  if ("xlab"%in%names(args)) xlabel <- args$xlab else
      xlabel <- ifelse(is.null(names(chr)),"Chromosome",names(chr))
  mtext(xlabel,1,line=xline,las=0,cex=0.5) }

#3-9-2013 MRC-Epid JHZ

#' Controls for mhtplot
#'
#' Parameter specification through function
#'
#' @param type Type of plot.
#' @param usepos A flag.
#' @param logscale A flag for log-scale.
#' @param base Base of log.
#' @param cutoffs Cutoffs of P-value, etc.
#' @param colors Colours for chromosomes.
#' @param labels Labels for chromosomes.
#' @param srt Rotation degrees.
#' @param gap Gap between data points.
#' @param cex Scaling factor of data points.
#' @param yline Vertical adjustment.
#' @param xline Horiztonal adjustment.
#'
#' @export
#' @return A list as above.

mht.control <- function(type="p", usepos=FALSE, logscale=TRUE, base=10, cutoffs=NULL, colors=NULL,
                        labels=NULL, srt=45, gap=NULL, cex=0.4, yline=3, xline=3)
               list(type=type, usepos=usepos, logscale=logscale, base=base, cutoffs=cutoffs, colors=colors,
                    labels=labels, srt=srt, gap=gap, cex=cex, yline=yline, xline=xline)

#' Controls for highlights
#'
#' Specification of highlights
#'
#' @param data Data.
#' @param colors Colors.
#' @param yoffset Y offset.
#' @param cex Scaling factor.
#' @param boxed Label in box.
#'
#' @export
#' @return A list as above.

hmht.control <- function(data=NULL, colors=NULL, yoffset=0.25, cex=1.5, boxed=FALSE)
                list(data=data,colors=colors,yoffset=yoffset,cex=cex,boxed=boxed)
