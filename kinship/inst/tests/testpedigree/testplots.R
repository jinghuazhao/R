## A SIMPLE PLOT ##

postscript('testplots.ps')

data1 <- scan(what=list(id=0, dadid=0, momid=0, affected=0, status=0, sex=0),
             file='test1.dat')

ped1 <- pedigree(data1$id, dadid=data1$dadid, momid=data1$momid, 
		 sex=2-data1$sex, status=data1$status, affected=data1$affected)

plot(ped1)
title('Test1 Data - Default Plot')

## SOME PEDIGREES THAT HAVE GIVEN US PROBLEMS IN THE PAST
## (DATA HAS BEEN ALTERED FROM THE ORIGINAL)

data2 <- read.table("test2.dat",
                    col.names=c("pedid", "id", "father", 
                      "mother", "sex", "status", "dead", "proband"))

fun2 <- function(i, ...) {
    temp <- data2[data2$pedid==i,]
    ptemp <- pedigree(temp$id, dadid=temp$father, momid=temp$mother, temp$sex,
		      status=temp$dead, affected=temp$status==2)
    plot(ptemp, ...)
    title(paste("Test 2 Data: Pedigree", i))
    invisible(ptemp)
    }

# Some interesting ones (broke the routines at times)
fun2(1)  #double marriage
fun2(2)  #should be able to straighten the connections more
         # (The springs are oscillating, plot(..., width=7) is perfect)
fun2(3)  # autohint not quite smart enough
fun2(4)  

## CREATE SOME INTERESTING PEDIGREES TO TRY AND STUMP THE FUNCTIONS

# A test pedigree, with odd cross-ties
#  This originally put kindepth() into an infinite recursion

test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
		    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
		    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
		    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))

cat('Print depth information\n')
print(kindepth(test1$id, test1$dad, test1$mom))
print(kindepth(test1$id, test1$dad, test1$mom, T))
cat('\n')

test1.ped <- pedigree(test1$id, test1$dad, test1$mom, test1$sex)
plot(test1.ped)
title('Test1.ped: odd cross-ties')

# Add a spouse with no children, and force them to be plotted
test2 <- rbind(test1, data.frame(id= 15, mom=0, dad=0, sex=1))

test2.ped <- pedigree(test2$id, test2$dad, test2$mom, test2$sex,
                      relation=matrix(c(9,15,4), nrow=1))
plot(test2.ped)
title('Test1.ped: add a spouse w/ no kids & force them to be plotted')

# Add in a pair of twins, and an inbred marriage
test3 <- rbind(test1,
               data.frame(id=c(15,16,17,18,19), mom=c(0,15,15,15,14),
                          dad=c(0,9,9,9,17), sex=c(1, 0,0,1,0)))
test3.ped <- pedigree(test3$id, test3$dad, test3$mom, test3$sex,
                      relation=matrix(c(16,17,17,18,1,2), nrow=2))
plot(test3.ped)
title('Test1.ped: Add in MZ twins (16 & 17) and an inbred marriage')

## MORE PROBLEM PEDIGREES

# ped1 has a sister marrying two brothers.  An early version of the
#  routines didn't plot this at all.  Autohint gets it right sometimes
#  (depending on the release), but only by luck.

ped1 <- pedigree(id=1:18, 
		 momid=c(0,0,0,2,2,0,4,4,4,5,5,7,7,7,7,11,14,16),
		 dadid=c(0,0,0,1,1,0,3,3,3,6,6,8,8,9,9, 9,13,15),
		 sex  =c(1,2,1,2,2,1,2,1,1,1,2,2,1,2,1, 2, 1, 1))
plot(ped1)
title('Ped1: sister marrying her own 2 brothers - try to confuse functions')

# A person marries 4 times, all have kids
ped2 <- pedigree(id=1:11, dadid=c(0,0,0,0,1,0,0,5,5,5,5),
		          momid=c(0,0,0,0,2,0,0,3,4,6,7),
			  sex=  c(1,2,2,2,1,2,2,1,2,1,1))
plot(ped2)
title('Person marries 4 times, all with kids')

# Same as ped2, but with parents for all of the marriages
# Here, hints don't help a lot.  The way in which the routine does its
#  processing precludes the joining that we want
ped3 <- pedigree(id=1:19, dadid=c(0,0,14,16,1,18,12,5,5,5,5,0,0,0,0,0,0,0,0),
		          momid=c(0,0,15,17,2,19,13,3,4,6,7,0,0,0,0,0,0,0,0),
		          sex=  c(1,2, 2, 2,1, 2, 2,1,2,1,1,1,2,1,2,1,2,1,2))

ped3$hints[c(2,13,15,17,19), 1] <- c(2,5,1,3,4)
ped3$hints[3,2] <- 5
plot(ped3)
title('Person marries 4 times, all with kids: add in all parents')

## SWAP PEOPLE AROUND

data4 <- read.table('test4.dat')
names(data4) <- c("pedid","id","father","mother","sex","status")
data4$sex <- data4$sex -1

ped1 <- pedigree(data4$id, data4$father, data4$mother, data4$sex, 
	        affected=data4$status ==2)
plot(ped1)
title('Test4.dat: Basic Plot')

ped1$hints[7,1] <- .1  #make subject 7 leftmost sib in the family
ped1$hints[5,1] <- .2  # make subject 5 be second from the left
plot(ped1)
title('Test4.dat: Move subjects 5 & 7')

# (The default hints in a pedigree are always between 1 and #siblings).

ped1$hints[ped1$id==15,1] <- .1   # Now make #15 leftmost in his family
ped1$hints[ped1$id==13,1] <- 100  # and #13 rightmost
plot(ped1)
title('Test4.dat: Move subjects 15 & 13')

data5 <- read.table('test5.dat',col.names=c('famid','upn','dadid','momid','sex'))

ok <- data5$famid==27
ped5.1 <- pedigree(data5$upn[ok],dadid=data5$dadid[ok],momid=data5$momid[ok],
                   sex=data5$sex[ok])
plot(ped5.1)
title('Test5.dat, family 27: 2nd husband not directly related to main pedigree') 

ok <- data5$famid==71
ped5.2 <- pedigree(data5$upn[ok],dadid=data5$dadid[ok],momid=data5$momid[ok],
                   sex=data5$sex[ok])
plot(ped5.2)
title('Test5.dat, family 71: 2-level loop (different generations')

dev.off()
