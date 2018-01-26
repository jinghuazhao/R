\name{plot.pedigree}
\alias{plot.pedigree}
\title{
  plot pedigrees
}
\description{
plot objects created with the function pedigree
}
\usage{
\method{plot}{pedigree}(x, id=x$id, sex=x$sex, status=x$status,
affected=x$affected, cex=1, col=rep(1, length(x$id)),
symbolsize=1, branch=0.6, packed=T, align=packed, width=8, density = c(5, 0, 100, 50),
mar=c(4.1,1,4.1,1), angle=c(90, 70, 50, 0), keep.par=F, \ldots)
}
\arguments{
\item{x}{
object created by the function pedigree.
}
\item{id}{
id variable - used for labelling.  
}
\item{sex}{
sex variable - used to determine which symbols are plotted. 
}
\item{status}{
can be missing.  If it exists, 0=alive/missing and 1=death.
}
\item{affected}{
variable, or matrix, of up to 4 columns representing 4 different
affected statuses.
}
\item{cex}{
controls text size.  Default=1.
}
\item{col}{
color for each id.  Default assigns the same color to everyone.
}
\item{symbolsize}{
controls symbolsize. Default=1.
}
\item{branch}{
defines how much angle is used to connect various levels of nuclear families.
}
\item{packed}{
default=T.  If T, uniform distance between all individuals at a given level.
}
\item{align}{
default=T.  Indicates extra effort should be spent trying to align parents
and children.  Set to F to speed up plotting.
}
\item{width}{
default=8.  For a packed pedigree, the minimum width allowed in the
realignment of pedigrees.
}
\item{density}{
defines density used in the symbols.  Takes up to 4 different values.
}
\item{mar}{
no comment(s)
}
\item{angle}{
defines angle used in the symbols.  Takes up to 4 different values.
}
\item{keep.par }{
Default = F, allows user to keep the parameter settings the same as they
were for plotting (useful for adding extras to the plot)
}
\item{\dots }{
Extra options that feed into the plot function.
}
}
\value{
returns points for each plot plus original pedigree.
}
\section{Side Effects}{
creates plot on current plotting device.
}
\seealso{
\code{\link{pedigree}}
}
\examples{
\dontrun{
# The original example
ptemp <- pedigree(id=d10$upn, dadid=d10$dadid, momid=d10$momid, 
                  sex=d10$sex, affected=d10$affect)

plot(ptemp)

col.founder <- rep(1,length(d10$affect))
col.founder[d10$id==1] <- 2

plot(ptemp, affected=cbind(d10$affect,d10$affect2,d10$affect3,d10$affect4),
     col=col.founder, id=paste(ptemp$id,'\n','(',d10$post,')',sep=''),
     angle=c(90,80,70,60), density=c(-1,90,70,50))

# The example to R News and Bioinformatics
# pedigree 10081 in GAW14
p1 <- scan(nlines=16,what=list(0,0,0,0,0,"",""))
 1   2   3  2  2  7/7  7/10 
 2   0   0  1  1  -/-  -/-  
 3   0   0  2  2  7/9  3/10 
 4   2   3  2  2  7/9  3/7  
 5   2   3  2  1  7/7  7/10 
 6   2   3  1  1  7/7  7/10 
 7   2   3  2  1  7/7  7/10 
 8   0   0  1  1  -/-  -/-  
 9   8   4  1  1  7/9  3/10 
10   0   0  2  1  -/-  -/- 
11   2  10  2  1  7/7  7/7 
12   2  10  2  2  6/7  7/7 
13   0   0  1  1  -/-  -/- 
14  13  11  1  1  7/8  7/8 
15   0   0  1  1  -/-  -/- 
16  15  12  2  1  6/6  7/7 

p2 <- as.data.frame(p1)
names(p2) <- c("id","fid","mid","sex","aff","GABRB1","D4S1645")
attach(p2)
sex <- sex-1
p3 <- pedigree(id,fid,mid,sex,aff)
pdf("10081.pdf")
par(xpd=TRUE)
plot(p3)
dev.off()
}
}
\keyword{hplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
