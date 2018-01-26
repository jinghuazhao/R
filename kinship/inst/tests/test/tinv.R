# look at inverses more closely
#  (I needed this when some of the other tests weren't being passed,
#  to figure out where in the decomposition/inversion/multiply process
#  the flaw was).

ch1 <- gchol(smat)
ch2 <- gchol(tmat)

inv1 <- solve(as.matrix(ch1))
inv2 <- solve(ch2,full=F)  #inverse of the cholesky, not of tmat
aeq(inv1, as.matrix(inv2))


#
# Now test the solution to a partial solve
#  We want to be able to transform a matrix to uncorrelated form
#  If tmat= LDL', and A is general, I want (D^{-1/2}) L^{-1} A
#
amat <- matrix(runif(5*nrow(tmat)), nrow=nrow(tmat))
xx1 <- diag(1/sqrt(diag(ch1))) %*% solve(as.matrix(ch1), amat) 
xx2 <- solve(ch2, amat, full=F)
aeq(xx1, xx2)

xx1 <- diag(1/sqrt(diag(ch1))) %*% solve(as.matrix(ch1), yy)
xx2 <- solve(ch2, yy, full=F)
aeq(xx1, xx2)
