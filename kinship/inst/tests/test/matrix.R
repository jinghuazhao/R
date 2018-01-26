# matrix multiplication
zz <- runif(13)
aeq(zz%*% smat, zz%*% tmat)
aeq(smat%*%zz, tmat%*% zz)

xx <- matrix(1:39, ncol=3)
aeq(smat %*% zz, tmat %*% zz)
aeq(t(xx) %*% smat, t(xx) %*% tmat)


amat <- tmat
amat@offdiag <- pi
bmat <- as.matrix(amat)

aeq(zz%*% amat, zz%*% bmat)
aeq(amat%*%zz, bmat%*% zz)


# Solve the right-hand side wrt a matrix
yy2 <- cbind(yy, -yy, yy+3)
zz1 <- solve(smat, yy2)
zz2 <- solve(tmat, yy2)
aeq(zz1, zz2)
aeq(zz2[,1], solve(tmat, yy))
