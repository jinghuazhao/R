.dat<-read.table('_ped2.csv',sep=',',row.names=NULL)
.dat$ped<-.dat$V1
.pnums<-unique(.dat$ped)
# hsrcolor()
pedplot<-function(pn) {
#### Legend Configuration
  par(fig=c(0,1,0,1/15),new=F)
  .leg<-read.table('_leg2.csv',sep=',',row.names=NULL)
  col.avail<-rep(1,length(.leg$V1))
  col.avail[.leg$V5==2] <- 2
  col.avail[.leg$V5==3] <- 3
  col.avail[.leg$V5==4] <- 4
  col.avail[.leg$V5==5] <- 10
  .legd<-pedigree(id=.leg$V1,dadid=.leg$V2,momid=.leg$V3,sex=.leg$V4)
  aff <- cbind(.leg$V7, .leg$V8, .leg$V9, .leg$V10)
  den <- rep(-1,7)
  ang <- rep(90,7)
  id <- paste(.leg$V6)
  plot(.legd,affected=aff,density=den,angle=and,col=col.avail,id=id,symbolsize=1.5,packed=T,cex=.5,ylim=c(-2.75,-1.98))
#### End Legend Configuration
  .in<-.dat[.dat$ped==pn,]
  V13<-.in$V13;V13[V13=='_']<-' '
  V14<-.in$V14;V14[V14=='_']<-' '
  V15<-.in$V15;V15[V15=='_']<-' '
  V16<-.in$V16;V16[V16=='_']<-' '
  V17<-.in$V17;V17[V17=='_']<-' '
  V18<-.in$V18;V18[V18=='_']<-' '
  V19<-.in$V19;V19[V19=='_']<-' '
  V20<-.in$V20;V20[V20=='_']<-' '
  col.avail<-rep(1,length(.in$ped))
  col.avail[.in$V7==2] <- 2
  col.avail[.in$V7==3] <- 3
  col.avail[.in$V7==4] <- 4
  col.avail[.in$V7==5] <- 10
  .twn<-read.table('_twin2.csv',sep=',',row.names=NULL,header=T)
  .twn<-.twn[.twn$ped==pn,2:4]
  ntwn<-nrow(.twn)
  browser()
  if (ntwn > 0) { .ptemp<-pedigree(id=.in$V2,dadid=.in$V3,momid=.in$V4,sex=.in$V5,relations=.twn)}
  if (ntwn == 0) {.ptemp<-pedigree(id=.in$V2,dadid=.in$V3,momid=.in$V4,sex=.in$V5)}
  par(fig=c(0,1,1/50,1))
  aff <- cbind(.in$V9,.in$V10,.in$V11,.in$V12)
  den <- c(-1,-1,-1,-1)
  ang <- c(90,90,90,90)
  id <- paste(.ptemp$id,'\n',V13,'\n',V14,'\n',V15,'\n',V16,'\n',V17,'\n',V18,'\n',V19,'\n',V20,sep='')
  plot(.ptemp,affected=aff,density=den,angle=ang,col=col.avail,id=id,status=.in$V6,symbolsize=.5,packed=F,cex=.4,keep.par=T)
  title(paste(.in$V8,'Pedigree',pn,sep=' '))
  pstamp()
}

postscript(file='allpeds.ps',horizontal=T)
library(kinship)
library(Hmisc) # pstamp()
for(pn in(76)) {pedplot(pn)}
dev.off()
