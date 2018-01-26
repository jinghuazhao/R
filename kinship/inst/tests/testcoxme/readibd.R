#
# See if we can create a bdsmatrix from the ibd style of file.
#
# data.restore('cdata.sdump')   # cdata already created by another test
# kmat <- makekinship(cdata$famid, cdata$gid, cdata$dadid, cdata$momid) # ditto

tmat <- as.matrix(kmat)
who <- (row(tmat) >= col(tmat) & tmat >0)
ibddata <- cbind((row(tmat))[who], (col(tmat))[who], tmat[who])

remake <- bdsmatrix.ibd(ibddata, idmap=cbind(1:1000,dimnames(kmat)[[1]]))

# The remake and kmat matrices won't look alike, since bdsmatrix.ibd has
#   in essence done what makefamid does: find singletons
newfam <- makefamid(cdata$gid, cdata$dadid, cdata$momid)
kmat2 <- makekinship(newfam, cdata$gid, cdata$dadid, cdata$momid)

# The kmat2 matrix still won't match, since it is in a different order!
temp2 <- bdsmatrix.reconcile(list(remake, kmat2), dimnames(kmat)[[1]])
aeq(temp2[[1]], temp2[[2]])


# Do some checks on kmat vs kmat2
xx <- match(cdata$momid, cdata$gid[newfam==0], nomatch=0)
all(xx==0) # should be true -- no one with kids has newfam set to 0
all(match(cdata$dadid, cdata$gid[newfam==0], nomatch=0) ==0)
all(cdata$momid[newfam==0] =='')
all(cdata$dadid[newfam==0] =='')


id2 <- match(cdata$famid, unique(cdata$famid)) -1
aeq(newfam[newfam!=0], id2[newfam!=0])  #other than zeros, everyone the same?
#
# kmat2 is about 2 times the size of kmat, in terms of storage, because of
#  the extra zeros for marry-ins.  See if they are essentially the same,
#  however.
# This isn't so easy: we need to expand out a submatrix and check it.

temp <- rep(T,10)
for (i in 1:10) {  # ten families
    idlist <- cdata$gid[newfam==i]
    row1 <- match(idlist, dimnames(kmat)[[1]])
    row2 <- match(idlist, dimnames(kmat2)[[1]])
    temp[i] <- aeq(as.matrix(kmat[row1,row1]), as.matrix(kmat2[row2,row2]))
    }
temp

rm(temp, id2, xx, row1, row2)
rm(tmat, who, ibddata, remake, kmat2, temp2, newfam)
