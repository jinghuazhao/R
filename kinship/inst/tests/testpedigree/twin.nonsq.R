id <-c(1:10)
dad<-c(rep(0,2),rep(1,8))
mom<-c(rep(0,2),rep(2,8))
sex<-c(1,2,rep(1,8))
postscript(file='Twins.ps',horizontal=T)
par(mfrow=c(3,2))

### Structure

 p<-pedigree(id,dad,mom,sex)
 plot(p)
 
### Sequential Twin ID

 p<-pedigree(id,dad,mom,sex,relations=matrix(c(3,4,2,4,3,2),2,3,byrow=T))
 plot(p)
 
### Twin Spaced 1 apart 
 
 p<-pedigree(id,dad,mom,sex,relations=matrix(c(3,5,2,5,3,2),2,3,byrow=T))
 plot(p)

### Twin Spaced 2 apart 
 
 p<-pedigree(id,dad,mom,sex,relations=matrix(c(3,6,2,6,3,2),2,3,byrow=T))
 plot(p)

### Twin Spaced 3 apart #### Here is where it mis-twinned!!!
 
 p<-pedigree(id,dad,mom,sex,relations=matrix(c(3,7,2,7,3,2),2,3,byrow=T))
 plot(p)

## try to fix it by changing the hints matrix
 p$hints[3:4] <- c(3,2)
 plot(p);title('swapped people 3 and 4')

##NOTE: Corrected Sept 22, 2003
dev.off()

id <-c(1:12)
dad<-c(rep(0,2),rep(1,8),0,3)
mom<-c(rep(0,2),rep(2,8),0,11)
sex<-c(1,2,rep(1,8),2,2)
p<-pedigree(id,dad,mom,sex,relations=matrix(c(3,7,2,7,3,2),2,3,byrow=T))

### Twin Spaced 4 apart #### Compare to 3 apart
 
p<-pedigree(id,dad,mom,sex,relations=matrix(c(3,8,2,8,3,2),2,3,byrow=T))
plot(p)


