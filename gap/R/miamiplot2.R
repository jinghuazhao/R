#' Miami Plot
#'
#' Creates a Miami plot to compare results from two genome-wide association analyses.
#'
#' @param gwas1 The first of two GWAS datasets to plot, in the upper region.
#' @param gwas2 The second of two GWAS datasets to plot, in the lower region.
#' @param name1 The name of the first dataset, plotted above the upper plot region. Defaults to \samp{"GWAS 1"}.
#' @param name2 The name of the second dataset, plotted below the lower plot region. Defaults to \samp{"GWAS 2"}.
#' @param chr1 The name of the column containing chromosome number in \env{gwas1}. Defaults to \samp{"chr"}.
#' @param chr2 The name of the column containing chromosome number in \env{gwas2}. Defaults to \samp{"chr"}.
#' @param pos1 The name of the column containing SNP position in \env{gwas1}. Defaults to \samp{"pos"}.
#' @param pos2 The name of the column containing SNP position in \env{gwas2}. Defaults to \samp{"pos"}.
#' @param p1 The name of the column containing p-values in \env{gwas1}. Defaults to \samp{"p"}.
#' @param p2 The name of the column containing p-values in \env{gwas2}. Defaults to \samp{"p"}.
#' @param z1 The name of the column containing z-values in \env{gwas1}. Defaults to \samp{"NULL"}.
#' @param z2 The name of the column containing z-values in \env{gwas2}. Defaults to \samp{"NULL"}.
#' @param sug The threshold for suggestive significance, plotted as a light grey dashed line.
#' @param sig The threshold for genome-wide significance, plotted as a dark grey dashed line.
#' @param pcutoff The p-value threshold below which SNPs will be ignored. Defaults to 0.1. It is not recommended to set this higher as it will narrow the central gap between the two plot region where the chromosome number is plotted.
#' @param topcols A vector of two colours to plot alternating chromosomes in for the upper plot. Defaults to green3 and darkgreen.
#' @param botcols A vector of two colours to plot alternating chromosomes in for the lower plot. Defaults to royalblue1 and navy.
#' @param yAxisInterval The interval between tick marks on the y-axis. Defaults to 5, 2 may be more suitable for plots with larger minimum p-values.
#' @export
#'
#' @return In addition to creading a Miami plot, the function returns a data frame containing x coordinates for chromosome start positions (required for \code{\link[gap]{labelManhattan}})
#'
#' @author Jonathan Marten
#' @note Extended to handle extreme P values.
#' @keywords GWAS Miami Manhattan
#' @examples
#' \dontrun{
#' # miamiplot2(gwas1, gwas2)
#' # chrmaxpos <- miamiplot2(gwas1, gwas2)
#' gwas <- within(mhtdata[c("chr","pos","p")], {z=qnorm(p/2)})
#' chrmaxpos <- miamiplot2(gwas,gwas,name1="Batch 2",name2="Batch 1",z1="z",z2="z")
#' labelManhattan(chr=c(2,16),pos=c(226814165,52373776),name=c("AnonymousGene","FTO"),
#'                gwas,gwasZLab="z",chrmaxpos=chrmaxpos)
#' }

miamiplot2 <- function(gwas1,gwas2,name1="GWAS 1",name2="GWAS 2",chr1="chr",chr2="chr",pos1="pos",pos2="pos",p1="p",p2="p",z1=NULL,z2=NULL,
                       sug=1e-5,sig=5e-8,pcutoff=0.1,topcols=c("green3","darkgreen"),botcols=c("royalblue1","navy"),yAxisInterval=5)
{
	## plots two manhattan plots top to bottom
	## TO DO: allow highlighting regions?

   #Input checks:
   if(length(which(names(gwas1)==chr1))==0){
      stop(paste0("Could not find column ",chr1," in gwas1"))
   }
   if(length(which(names(gwas2)==chr2))==0){
      stop(paste0("Could not find column ",chr2," in gwas2"))
   }
   if(length(which(names(gwas1)==pos1))==0){
      stop(paste0("Could not find column ",pos1," in gwas1"))
   }
   if(length(which(names(gwas2)==pos2))==0){
      stop(paste0("Could not find column ",pos2," in gwas2"))
   }
   if(length(which(names(gwas1)==p1))==0){
      stop(paste0("Could not find column ",p1," in gwas1"))
   }
   if(length(which(names(gwas2)==p2))==0){
      stop(paste0("Could not find column ",p2," in gwas2"))
   }

	# make data frames with just the required cols from the input. Probably v. inefficient but hey.
        if(is.null(z1))
	dat1 <- within(data.frame("chr"=gwas1[,which(names(gwas1)==chr1)], "pos"=gwas1[,which(names(gwas1)==pos1)], "p"=gwas1[,which(names(gwas1)==p1)]),
                {log10p=-log10(p)})
        else
	dat1 <- within(data.frame("chr"=gwas1[,which(names(gwas1)==chr1)], "pos"=gwas1[,which(names(gwas1)==pos1)], "z"=gwas1[,which(names(gwas1)==z1)]),
                {log10p=-log10p(z)})
        if(is.null(z2))
	dat2 <- within(data.frame("chr"=gwas2[,which(names(gwas2)==chr2)], "pos"=gwas2[,which(names(gwas2)==pos2)], "p"=gwas2[,which(names(gwas2)==p2)]),
                {log10p=-log10(p)})
        else
	dat2 <- within(data.frame("chr"=gwas2[,which(names(gwas2)==chr2)], "pos"=gwas2[,which(names(gwas2)==pos2)], "z"=gwas2[,which(names(gwas2)==z2)]),
                {log10p=-log10p(z)})

	# make object with cumulative chr start and end positions - to allow plotting of GWASes with different numbers of SNPs without things being misaligned.
	chrmaxpos <- data.frame("chr"=1:22,"maxpos"=NA,"genomestartpos"=NA)
	for(i in 1:22){
		chrmaxpos$maxpos[i] <- max(c(dat1$pos[dat1$chr==i],dat2$pos[dat2$chr==i]),na.rm=TRUE)
		if(i==1){
			chrmaxpos$genomestartpos[i] <- 0
		} else {
			chrmaxpos$genomestartpos[i] <- chrmaxpos$genomestartpos[i-1] + chrmaxpos$maxpos[i-1]
		}
	}
	# Positions for plotting chr labels in centre of chromosomes
	chrmaxpos$labpos <- chrmaxpos$genomestartpos + 0.5*chrmaxpos$maxpos

	# Get 'genome' positions of each SNP
	dat1$gpos <- NA
	dat2$gpos <- NA

	for(ch in 1:22){
		vec <- which(dat1$chr==ch)
		dat1$gpos[vec] <- dat1$pos[vec]+chrmaxpos$genomestartpos[ch]
		vec2 <- which(dat2$chr==ch)
		dat2$gpos[vec2] <- dat2$pos[vec2]+chrmaxpos$genomestartpos[ch]
	}

	# Parameters for graph - ylim, vector of SNPs to plot and colours
	#maxp <- max(-log10(min(dat1$p)),-log10(min(dat2$p))) +2
	maxp <- max(max(dat1$log10p,na.rm=TRUE),max(dat2$log10p,na.rm=TRUE))
	plotvec1 <- which(dat1$log10p>=-log10(pcutoff))
	plotvec2 <- which(dat2$log10p>=-log10(pcutoff))
	col1 <- rep(topcols,22)[dat1$chr[plotvec1]]
	col2 <- rep(botcols,22)[dat2$chr[plotvec2]]

	#Plot graph
	plot(dat1$gpos[plotvec1], dat1$log10p[plotvec1], pch=20, cex=0.8,col=col1,ylim=c((-maxp-2), (maxp+8)),xaxt="n",yaxt="n",ylab="-log10(P-value)",xlab="",bty="n")
	points(dat2$gpos[plotvec2],-dat2$log10p[plotvec2], pch=20, cex=0.8,col=col2)

	# Add significance lines and chromosome numbers
	abline(h=-log10(sug),col="gray70",lty=2)
	abline(h=log10(sug),col="gray70",lty=2)
	abline(h=-log10(sig),col="gray50",lty=2)
	abline(h=log10(sig),col="gray50",lty=2)
	#axis(1,at=chrmaxpos$genomestartpos,labels=rep("",22),pos=0)
	text(chrmaxpos$labpos,0,as.character(1:22),cex=0.8)
	axisLabNum <- floor(maxp/yAxisInterval) # number of labels needed for each side of axis
	axisLabAt <-  seq(-yAxisInterval*axisLabNum, yAxisInterval*axisLabNum, yAxisInterval) # position of each axis tick mark
	axisLabels <- c(seq(axisLabNum*yAxisInterval,0,-yAxisInterval),seq(yAxisInterval,axisLabNum*yAxisInterval,yAxisInterval))

	axis(2, at=axisLabAt, labels=axisLabels)

	# Add titles to each side of the plot
	mtext(name1,side=3,font=2)
	mtext(name2,side=1,font=2)

	# Returns the position-modifier values for each chromosome - required for labelManhattan function
	return(chrmaxpos)
}
