#
# This is a test for lmekin, using data sets created from the GAW
#   simulations.  The IBD matrix is read in from a data set created
#   by SOLAR.  Becuase Solar mixes up the order of subjects in the
#   families, this is a good test of bdsmatix.reconcile
#
temp1 <- read.table('data.gaw1', header=F, skip=1,
                   row.names=NULL, sep=',',  
                   col.names=c('famid', 'id', 'father', 'mother', 'sex'))

temp2 <- read.table('data.gaw2', header=F, skip=1, 
                    row.names=NULL, sep=',', 
                    col.names=c('famid', 'id', 'age', 'ef1', 'ef2', 
                                'q1', 'q2', 'q3', 'q4', 'q5'))

aeq(temp1$id, temp2$id)
aeq(temp1$famid, temp2$famid)
adata <- cbind(temp1, temp2[,-(1:2)])
rm(temp1, temp2)

newfam <- makefamid(adata$id, adata$father, adata$mother)
kmat <- makekinship(newfam, adata$id, adata$father, adata$mother)

#
# Read an ibd file
#   It corresponds to chromosome 6, position 90
# The pedindex file is also created by Solar, and contains the
#   original identifiers for each "solar" id.
#
temp <-matrix( scan("pedindex.out"), ncol=9, byrow=T)
pedindex <- temp[,c(1,9)]

temp <- read.table("data.solar", row.names=NULL, 
                       col.names=c("id1", "id2", "x", "dummy"))
ibd6.90 <- bdsmatrix.ibd(temp$id1, temp$id2, temp$x, pedindex)

#
# First test, see if I get the same answer as lme
#
all(diff(adata$famid) >=0)  #make sure the data set is in famid order
fit1 <- lme(age ~ q1, adata[!is.na(adata$age),], random=~1|famid, method='ML')

fmat <- bdsBlock(adata$id, adata$famid)
fit2 <- lmekin(age ~ q1, adata, random= ~1|id, varlist=list(fmat))
# this will fail 
#   fit3 <- lmekin(age ~ q1, adata, random= ~1|famid)
# as lmekin still requires exactly one random effect per subject

aeq((fit1$coefficients)$fixed, coef(fit2)$fixed)
aeq(fit1$sigma^2, fit2$theta[2])
aeq(fit1$logLik, fit2$loglik)

#
# Now for a real kinship fit, and a real ibd one
#
kfit1 <- lmekin(age ~1, adata, random=~1|id, varlist=list(kmat))
kfit2 <- lmekin(age ~1, adata, random=~1|id, varlist=list(kmat, ibd6.90))

kfit1
kfit2
