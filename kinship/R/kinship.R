# $Id: kinship.s,v 1.5 2003/01/04 19:07:53 Therneau Exp $
#
# Create the kinship matrix, using the algorithm of K Lange,
#  Mathematical and Statistical Methods for Genetic Analysis,
#  Springer, 1997, p 71-72.
#
# The rows (cols) of founders are just .5 * identity matrix, no further
#    processing is needed for them.
#  Parents must be processed before their children, and then a child's
#    kinship is just a sum of the kinship's for his/her parents.
#
kinship <- function(id, father.id, mother.id) {
    n <- length(id)
    if (any(duplicated(id))) stop("All id values must be unique")
    kmat <- diag(n+1) /2
    kmat[n+1,n+1]    <- 0  #if A and B both have "unknown" dad, this ensures
                           # that they won't end up 'related' in the matrix

    pdepth <- kindepth(id, father.id, mother.id)
    # id number "n+1" is a placeholder for missing parents
    mrow <- match(mother.id, id, nomatch=n+1) #row number of the mother
    drow <- match(father.id, id, nomatch=n+1) #row number of the dad 

    # Those at depth==0 don't need to be processed
    # Subjects with depth=i must be processed before those at depth i+1.
    # Any parent is guarranteed to be at a lower depth than their children
    #  The inner loop on "i" can NOT be replaced with a vectorized expression:
    # sibs' effect on each other is cumulative.
    for (depth in 1:max(pdepth)) {
	indx <- (1:n)[pdepth==depth]
	for (i in indx) {
	    mom <- mrow[i]
	    dad <- drow[i]
	    kmat[i,]  <- kmat[,i] <- (kmat[mom,] + kmat[dad,])/2
	    kmat[i,i] <- (1+ kmat[mom,dad])/2
	    }
	}
    
    kmat <- kmat[1:n,1:n]
    dimnames(kmat) <- list(id, id)
    kmat
    }
