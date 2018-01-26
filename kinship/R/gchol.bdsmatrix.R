# $Id: gchol.bdsmatrix.s,v 1.3 2003/02/05 19:29:12 Therneau Exp $
#
# Cholesky decompostition for block-diagonal square matrices
#
setClass('gchol.bdsmatrix',
	 representation(blocksize = 'integer',
			blocks    = 'numeric',
			rmat      = 'matrix',
			rank      = 'integer',
			.Dim      = 'integer',	
			.Dimnames = 'list or NULL'))

gchol.bdsmatrix <- function(x, tolerance=1e-9) {
    if (class(x) != 'bdsmatrix') stop ("Bad argument")
    if (x@offdiag !=0) return(gchol(as.matrix(x)))
    dd <- x@.Dim
    if (length(x@rmat) >0) {
        nc <- ncol(x@rmat)
        temp <- .C("gchol_bds", as.integer(length(x@blocksize)),
		                as.integer(x@blocksize),
                                as.integer(dd),
                                dmat= as.double(x@blocks),
                                rmat= as.double(x@rmat),
                                flag= as.double(tolerance), PACKAGE="kinship")
        newr <- matrix(temp$rmat, ncol=nc)
        if (nc>1) {  
	    # The C-routine doesn't zero out t(r) above the diagonal
	    #   (the lower right corner)
	    d3 <- sum(x@blocksize)
            for (i in 1:(nc-1)) newr[(1+d3+i):dd[1],i] <- 0
            }
        new('gchol.bdsmatrix', blocksize=as.integer(x@blocksize), blocks=temp$dmat, 
	                rmat=as.matrix(newr), .Dim=as.integer(x@.Dim), 
                        rank=as.integer(temp$flag),
	                .Dimnames=x@.Dimnames)
        }
    else {
        temp <- .C("gchol_bds", as.integer(length(x@blocksize)),
                                as.integer(x@blocksize),
                                as.integer(dd),
                                blocks =as.double(x@blocks),
                                as.double(0),
                                flag= as.double(tolerance), PACKAGE="kinship")

        new('gchol.bdsmatrix', blocksize=as.integer(x@blocksize), blocks=temp$blocks, 
	                rmat=as.matrix(numeric(0)), .Dim=as.integer(x@.Dim), 
                        rank=as.integer(temp$flag),
	                .Dimnames=x@.Dimnames)
        }
    }

# 
#  return L, from the LDL' decompostion
#
as.matrix.gchol.bdsmatrix <- function(x, ones=F, ...){
    dd <- x@.Dim
    n <- dd[1]
    newmat <- matrix(0., n, n, dimnames=x@.Dimnames)
    temp <- .C('bdsmatrix_index2',
	       as.integer(length(x@blocksize)),
	       as.integer(x@blocksize),
	       rows= integer(length(x@blocks)),
	       cols= integer(length(x@blocks)), PACKAGE="kinship")
    rindex <- match(temp$rows, 1:n, nomatch=0)
    cindex <- match(temp$cols, 1:n, nomatch=0)
    newmat[cbind(rindex, cindex)] <- x@blocks
    if (length(x@rmat)){
	d3 <- sum(x@blocksize)
	newmat[(d3+1):n, ]<-  t(x@rmat)
	}
    if (ones) diag(newmat) <- 1
    newmat
    } 

setAs('gchol.bdsmatrix', 'matrix', 
      function (from) as.matrix.gchol.bdsmatrix(from))

setMethod('diag', signature(x='gchol.bdsmatrix'),
    function(x,nrow=1,ncol=1) {
	d <- x@.Dim
	d3 <- sum(x@blocksize)
	temp <- .C('bdsmatrix_index1',
		   as.integer(length(x@blocksize)),
		   as.integer(x@blocksize),
		   as.integer(c(0,1,0)),
		   as.integer(d3),
		   as.integer(1:d3 -1),
		   integer(1),
		   indexb = integer(d3),
		   integer(1), PACKAGE="kinship")$indexb

	if (length(x@rmat) > 0) {
	    temp2 <- seq(from=d3+1, by= d[2]+1, length= d[1] - d3)
	    c(x@blocks[temp], x@rmat[temp2])
	    }
	    else x@blocks[temp]
	})

setMethod('dim', 'gchol.bdsmatrix', function(x) x@.Dim)

setMethod('show', 'gchol.bdsmatrix',
          function(object) show(as.matrix(object)))
    

# The subscript method is almost identical to that for bdsmatix,
#   the main difference being that the bdsmatrix method fills in symmetry
#   when the result is not sparse
setMethod('[', 'gchol.bdsmatrix', 
 function(x, i, j, ..., drop=T) {
    if (class(x) != 'gchol.bdsmatrix') stop("Must be a gchol.bdsmatrix object")
    allargs <- list(...)
#   if (nDotArgs(...) !=2) stop("Two subscripts are required")
#   if (length(allargs)!=2) stop("Two subscripts are required")

    nblock <- length(x@blocksize)
    d <- x@.Dim
    d3 <- sum(x@blocksize)
    d4 <- length(x@blocks)
#   if (any(..1 > d[1]))
    if ((i > d[1]))
        stop(paste("Array subscript (", max(i), 
                   ") out of bounds, should be at most ", d[1], sep=''))

    if (any(j > d[2]))
        stop(paste("Array subscript (", max(j), 
                   ") out of bounds, should be at most ", d[2], sep=''))

#   rows <- (1:d[1])[..1]
#   cols <- (1:d[2])[..2]
    rows <- (1:d[1])[i]
    cols <- (1:d[2])[j]
    
    # The only case where the result is still a Cholesky is if you grab the
    #  first k rows/cols
    if (length(rows)==length(cols) && all(rows==cols) && 
                                      all(rows== 1:(length(rows)))) {
        brows <- rows[rows <= d3]  #the block-diagonal portion
        brlen <- length(brows)
        # The result will be block-diagonal symmetric
        # Note: we don't allow for reordering the row/col indices: too hard
        #   to keep track of what's happening
        temp <- .C('bdsmatrix_index1',
                   as.integer(nblock),
                   bsize = as.integer(x@blocksize),
                   as.integer(c(0,0,1)),
                   as.integer(brlen),
                   as.integer(brows -1),
                   integer(1),
                   integer(1),
                   indexc = integer(d4), PACKAGE="kinship")

        x@blocksize <- temp$bsize[temp$bsize>0]
        x@blocks <- x@blocks[temp$indexc]
        if (length(x@rmat)) {
            if (any(cols>d3)) x@rmat <- x@rmat[rows, cols[cols>d3]-d3, drop=F]
            else              x@rmat <- numeric(0)
            }
        if (!is.null(x@.Dimnames)) 
                x@.Dimnames <- list((x@.Dimnames[[1]])[rows],
                                    (x@.Dimnames[[2]])[cols])
        x@.Dim <- rep(length(rows),2)
        dd <- diag(x)
        x@rank <- sum(dd!=0)
        x
        }

    else { # The result is not a gchol.bdsmatrix object
        brows <- rows[rows <= d3]  #the block-diagonal portion
        brlen <- length(brows)
        bcols <- cols[cols <= d3]
        bclen <- length(bcols)
        if (brlen>0 && bclen>0) {
            bdsize <- .Options$bdsmatrixsize
            if (is.null(bdsize)) bdsize <- 1000
            if (prod(x@.Dim) > bdsize )
                  stop("Automatic conversion would create too large a matrix")
            # I need to subscript the block diagonal portion
            #  index2 is the rows() and cols() function for the block portion
            temp <- .C('bdsmatrix_index2',
                       as.integer(nblock),
                       as.integer(x@blocksize),
                       rows= integer(length(x@blocks)),
                       cols= integer(length(x@blocks)), PACKAGE="kinship")
            newmat <- matrix(x@offdiag, brlen, bclen)
            rindex <- match(temp$rows, brows, nomatch=0)
            cindex <- match(temp$cols, bcols, nomatch=0)
            keep <- (rindex>0 & cindex >0)  #the row/col is one we want to keep
            if (any(keep)) 
                newmat[cbind(rindex[keep], cindex[keep])] <- x@blocks[keep]

            if (length(x@rmat)) {
                if (any(rows > d3)) {
                    newmat <- rbind(newmat, t(x@rmat[bcols, rows[rows>d3]-d3]))
                    }
                if (any(cols > d3)) {
                    newmat <- cbind(newmat, x@rmat[rows, cols[cols>d3]-d3])
                    }
                }
            }
        else newmat <-x@rmat[rows, cols[cols>d3]-d3, drop=F]

        if (!is.null(x@.Dimnames)) 
            dimnames(newmat)<- list((x@.Dimnames[[1]])[rows],
                                    (x@.Dimnames[[2]])[cols])
        newmat[,,drop=drop]
        }
    })


#
# This follows bdsmatrix for multiplication in the other direction
# added by JH Zhao 19/3/2004
#

setMethod("%*%", signature(x='matrix', y='gchol.bdsmatrix'),
    function(x, y) {
        t(y%*% t(x))
        })
    
setMethod("%*%", signature(x='numeric', y='gchol.bdsmatrix'),
    function(x, y) {
        t(y%*% x)
        })
