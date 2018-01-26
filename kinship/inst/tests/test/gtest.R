# Create a matrix that is symmetric, but not positive definite
#   The first one, temp, has column 6 redundant with cols 1-5
temp <- smat[c(1:5, 5:10), c(1:5, 5:10)]
ch1  <- gchol(temp)
aeq(diag(ch1)[6], 0)  # Check that it has a zero in the proper place
ginv <- solve(ch1)    # see if I get a generalized inverse
aeq(temp %*% ginv %*% temp, temp)
aeq(ginv %*% temp %*% ginv, ginv)

# Now create one that is negative definite 
ch2 <- gchol(smat)
temp2 <- as.matrix(ch2)
temp3 <- diag(ch2) * rep(c(1, -1), length=nrow(smat))
xmat  <- temp2 %*% diag(temp3) %*% t(temp2)
xmat  <- (xmat + t(xmat))/2  #work out round-off errors
ch3 <- gchol(xmat)

aeq(diag(ch3), temp3)
aeq(as.matrix(ch3), temp2)
