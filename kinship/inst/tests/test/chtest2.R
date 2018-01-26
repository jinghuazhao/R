#
# Inverse of the matrix:  
#
inv1 <- solve(smat)
inv2 <- as.matrix(solve(tmat)) # the result is a full, non-sparse matrix
aeq(inv1, inv2)

inv3 <- solve(gchol(tmat))         #sparse version, not all parts will be there
inherits(inv3, 'bdsmatrix')         #This should be true
aeq(inv3@blocksize, tmat@blocksize) # Should be the same shape at tmat
inv3 <- as.matrix(inv3)		    # What is returned should be correct
aeq(inv1[1:3,1:3], inv3[1:3, 1:3])
aeq(inv1[4:5,4:5], inv3[4:5, 4:5])
aeq(inv1[6:7,6:7], inv3[6:7, 6:7])
aeq(inv1[8:11,8:11], inv3[8:11, 8:11])
aeq(inv1[,12:13], inv3[, 12:13])    # and rmat the same too
