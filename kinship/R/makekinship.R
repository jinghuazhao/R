# $Id: makekinship.s,v 1.6 2003/03/20 13:59:27 Therneau Exp $
#
# Make a sparse (bdsmatrix) kinship matrix
#  The families are processed one by one.
#  The final row order of the matrix will not necessarily be that of the
#    input data set -- pay attention to the final dimnames in order to 
#    select subsets!
#  Subjects with a family id of "unrelated" are assumed to be unrelated
#    uniques (often marry-ins with no children).  This allows for efficiency
#    in computation, by making them each a small block, without having to 
#    create a unique family id to each.  They are always put first in 
#    the matrix.
#
makekinship <- function(famid, id, father.id, mother.id, unrelated=0) {
    n <- length(famid)
    if (length(id)    != n) stop("Mismatched lengths: famid and id")
    if (length(mother.id) != n) stop("Mismatched lengths: famid and mother.id")
    if (length(father.id) != n) stop("Mismatched lengths: famid and father.id")
    if (any(is.na(famid)))  stop("One or more subjects with missing family id")
    if (any(is.na(id)))     stop("One or more subjects with a missing id")
    if (any(famid <0))      stop("Invalid family id, must be >0")
    if (any(duplicated(id))) stop("Subject id values must be distinct")

    famlist <- sort(unique(famid))  #same order as the counts table
    idlist <- id            # will be overwritten, but this makes it the
                            #  correct data type and length
    counts <- table(famid)
    cumcount <- cumsum(counts)    
    
    if (any(famid==unrelated)) {
	# Assume that those with famid of 0 are unrelated uniques
	#   (usually the marry-ins)
        temp <- match(unrelated, names(counts))
	nzero <- counts[temp]    
	counts <- counts[-temp]
	famlist <- famlist[famlist != unrelated]
	idlist[1:nzero] <- id[famid== unrelated]
	cumcount <- cumsum(counts) + nzero
	}
    else nzero <- 0
    
    blockn <- counts*(counts+1)/2   #size of storage for each block
    n2 <- sum(blockn) 	    # total amount needed
    bdata <- double(n2)
    j <- cumsum(blockn)     
    for (i in 1:length(counts)) {
	who <- (famid == famlist[i])
        if (sum(who) ==1) bdata[j[i]] <- 0.5  # family of size 1
        else {
            temp <- kinship(id[who], mother.id[who], father.id[who])
            bdata[seq(to=j[i], length=blockn[i])] <- temp[row(temp)>=col(temp)]
            }
	idlist[seq(to=cumcount[i], length=counts[i])] <- id[who]
	}

    bdsmatrix(blocksize=c(rep(1,nzero), counts),
	      blocks =  c(rep(.5,nzero), bdata),
	      dimnames=list(idlist, idlist))
    }
    
