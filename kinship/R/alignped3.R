# $Id: alignped3.s,v 1.1 2003/02/05 00:08:44 atkinson Exp $ 

# The real work is in three routines:
#   alignped1: called with a single subject, returns the subtree founded
#      on this subject, as though it were the only subtree.
#   alignped2:  called with a set of sibs: calls align multiple times and
#      merges the result
#   alignped3: glue two results together
#
# The work here is keeping positions up to date,
#   and recognizing left/right people that join up
#

alignped3 <- function(x1, x2, packed, space=1) {
    maxcol <- max(x1$n + x2$n)
    maxlev <- length(x1$n)
    n1 <- max(x1$n)   # These are always >1
    n  <- x1$n + x2$n

    nid <- matrix(0, maxlev, maxcol)
    nid[,1:n1] <- x1$nid
    
    pos <- matrix(0.0, maxlev, maxcol)
    pos[,1:n1] <- x1$pos

    fam <- matrix(0, maxlev, maxcol)
    fam[,1:n1] <- x1$fam
    fam2 <- x2$fam

    if (!packed) {
	# The aligned case: instead of sliding each row separately
	#   we think of the two trees as solids, slide till the solids touch
	# First task : compute the amount of the slide
	slide <- 0
        for (i in 1:maxlev) {
            n1 <- x1$n[i]
            n2 <- x2$n[i]
	    if (n1 >0 & n2 >0) {
		if (nid[i,n1] == x2$nid[i,1])
			temp <- pos[i, n1] - x2$pos[i,1]
		else    temp <- space + pos[i, n1] - x2$pos[i,1]
		if (temp > slide) slide <- temp
		}
	    }
	}

    # Now, merge the two trees
    for (i in 1:maxlev) {
	n1 <- x1$n[i]
	n2 <- x2$n[i]
	if (n2 >0) {   # If there is someone to add to the row
	    if (n1>0 && (nid[i,n1] == floor(x2$nid[i,1]))) {
		# same person in the two subtrees, touching: combine
		n[i] <- n[i] -1
		if (n2==1) {  #only 1 person: easy case
		    # whichever instance is the 'spouse' will have fam=0, 
		    #   keep the `sib' value, which is >0
		    fam[i,n1] <- max(fam[i,n1], fam2[i,1])
		    if (!packed) pos[i,n1] <- x2$pos[i,1] + slide
		    next;  
		    }
		nid[i, n1] <- max(nid[i,n1], x2$nid[i,1]) #preserve a ".5"
		zz <- 2:n2      
		nid[i, n1 + zz -1] <- x2$nid[i, zz]
		fam[i, n1 + zz -1] <- fam2[i,zz]
		fam[i, n1] <- max(fam[i,n1], fam2[i,1]) #parent of overlap
		if (i<maxlev) {
		    # adjust the pointers of any children (look ahead)
		    temp <- fam2[i+1,]
		    fam2[i+1,] <- ifelse(temp==0, 0, temp + n1 -1)
		    }
		if (packed) pos[i, n1 + zz -1] <- x2$pos[i,zz] + pos[i,n1]
		else {
                    pos[i, n1 + zz -1] <- x2$pos[i,zz] + slide
		    # put the person in common under their parents
		    #  or, if a double spouse, midway between
		    if (fam[i,n1]==0) {
			if (fam2[i,1]==0)
				pos[i, n1]<- (pos[i,n1] + x2$pos[i,1]+ slide)/2
			else    pos[i, n1]<- x2$pos[i,1] + slide
			}
		    }
		}
            else {  #no overlapping subject
                zz <- 1:n2      
		nid[i, n1 + zz] <- x2$nid[i, zz]
		fam[i, n1 + zz] <- fam2[i, zz]
		if (i<maxlev) {
		    temp <- fam2[i+1,]
		    fam2[i+1,] <- ifelse(temp==0, 0, temp + n1)
		    }
		if (packed) {
		    if (n1 >0)
			    pos[i, n1 + zz] <- x2$pos[i,zz] + pos[i,n1] + space
		    else    pos[i, zz]      <- x2$pos[i,zz]
		    }
		else  pos[i, n1 + zz] <- x2$pos[i,zz] + slide
		}
	    }
	}

    # sometimes, the merging of an id makes the result matrix narrower
    if (max(n) < maxcol) {
	maxcol <- max(n)
	nid <- nid[,1:maxcol]
	pos <- pos[,1:maxcol]
	fam <- fam[,1:maxcol]
	}

    list(n=n, nid=nid, pos=pos, fam=fam)
    }
