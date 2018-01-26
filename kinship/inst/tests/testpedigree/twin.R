library(kinship)

#source(tilde('~atkinson/research/kinship/plot.pedigree.s'))
#source('../plot.pedigree.s')

postscript(file='TwinTest.ps')

id<-c(1:5)
sex<-c(2,rep(1,4))
fa<-c(0,0,rep(2,3))
mo<-c(0,0,rep(1,3))

### Identical Twins
 twin<-matrix(c(3,4,4,3,1,1),2,3)

 my.ped<-pedigree(id=id,dadid=fa,momid=mo,sex=sex,relations=twin)

 plot(my.ped,packed=F,symbolsize=.5,cex=.5)
 
### Fraternal Twins
 twin<-matrix(c(3,4,4,3,2,2),2,3)

 my.ped<-pedigree(id=id,dadid=fa,momid=mo,sex=sex,relations=twin)

 plot(my.ped,packed=F,symbolsize=.5,cex=.5)

### 3 Identical Triplets
 twin<-matrix(c(3,4,5,5,3,4,1,1,1),3,3)

 my.ped<-pedigree(id=id,dadid=fa,momid=mo,sex=sex,relations=twin)

 plot(my.ped,packed=F,symbolsize=.5,cex=.5)
 
### 3 Fraternal Triplets
 twin<-matrix(c(3,4,5,5,3,4,2,2,2),3,3)

 my.ped<-pedigree(id=id,dadid=fa,momid=mo,sex=sex,relations=twin)

 plot(my.ped,packed=F,symbolsize=.5,cex=.5)

dev.off()

 

## use relation to modify plot to encorporate twins
## mztwins have a code of 1, dztwins have a code of 2

test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
		    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
		    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
		    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))

test3 <- rbind(test1,
               data.frame(id=c(15,16,17,18,19), mom=c(0,15,15,15,14),
                          dad=c(0,9,9,9,17), sex=c(1, 0,0,1,0)))

test3.ped <- pedigree(test3$id, test3$dad, test3$mom, test3$sex,
                      relation=matrix(c(16,17,17,18,1,2), nrow=2))
plot(test3.ped)
title('Test1.ped: Add in MZ twins (16 & 17) and an inbred marriage')


## check Curt's latest puzzle

test4 <- read.table('jhu.csv',sep=',',row.names=NULL,header=T)
twins <- read.table('_twin.csv',sep=',',row.names=NULL,header=T)

## twins on the ends of the family 
ok <- test4$ped==1736
ptemp1736 <- pedigree(id=test4$id[ok],dadid=test4$dad[ok],momid=test4$mom[ok],
                  sex=test4$sex[ok],relations=twins[twins$ped==1736,2:4])
plot(ptemp1736,packed=F)

## one twin married 3 times
ok <- test4$ped==2312
ptemp2312 <- pedigree(id=test4$id[ok],dadid=test4$dad[ok],momid=test4$mom[ok],
                  sex=test4$sex[ok],relations=twins[twins$ped==2312,2:4])
plot(ptemp2312)


## add in twins in 2 different nuclear families
ok <- test4$ped==1736

twins2 <- rbind(twins,matrix(c(1736,12,13,2,
                               1736,47,48,2,
                               1736,16,4,1),ncol=4,byrow=T))

ptemp1736 <- pedigree(id=test4$id[ok],dadid=test4$dad[ok],momid=test4$mom[ok],
                  sex=test4$sex[ok],relations=twins2[twins2$ped==1736,2:4])
plot(ptemp1736,packed=F)
