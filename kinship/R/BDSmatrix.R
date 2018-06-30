# $Id: bdsmatrix.s,v 1.3 2002/12/26 22:54:48 Therneau Exp $
rmat <- NULL
if(version$major>=2) setClassUnion("list or NULL",c("list","NULL"))
if(version$major==1 && version$minor>=8.0) setClassUnion("list or NULL",c("list","NULL"))
if(version$major==1 && version$minor <8.0) {
  setClass("list or NULL")
  setIs("list","list or NULL")
  setIs("NULL","list or NULL")
}
setClass('bdsmatrix',
	 representation(blocksize = 'integer',
			blocks    = 'numeric',
			rmat      = 'matrix',
			offdiag   = 'numeric',
			.Dim='integer', .Dimnames="list or NULL"))

setMethod('Math', signature(x='bdsmatrix'),
	  function(x) {
	      x@offdiag <- callGeneric(x@offdiag)
	      x@blocks  <- callGeneric(x@blocks)
	      x@rmat    <- callGeneric(x@rmat)
	      x })

setMethod('Math2', signature(x='bdsmatrix'),
	  function(x, digits) {
	      x@offdiag <- callGeneric(x@offdiag, digits)
	      x@blocks  <- callGeneric(x@blocks, digits)
	      x@rmat    <- callGeneric(x@rmat, digits)
	      x })

# For the summary method, we need to count the number of zeros (the off
#  diagonal elements of the block portion) that are not stored, and put them
#  into the computation.  This is trivial for min, max, and etc, but for
#  means and products we have written them out as weighted computations.
#  (The number of off-diagonal elements can be in the billions, rep() would
#   not be wise).
# Per a note from Bill Dunlap, max(c(x1,x2,x3)) is faster than max(x1, x2, x3),
#  when x1, x2, etc are all numeric.  (Up to 50 times faster!)
#
setMethod("max", signature("bdsmatrix"),
	  function(x, ..., na.rm=F) {
	      if (length(x@rmat))
	           max(c(x@offdiag, x@blocks, x@rmat), na.rm=na.rm)
	      else
	           max(c(x@offdiag, x@blocks), na.rm=na.rm)
	      })

setMethod("min", signature("bdsmatrix"),
	  function(x, ..., na.rm=F ) {
	      if (length(rmat))
	           min(c(x@offdiag, x@blocks, x@rmat), na.rm=na.rm)
	      else
	           min(c(x@offdiag, x@blocks), na.rm=na.rm)
	      }, valueClass="numeric")

setMethod('range',signature('bdsmatrix'),
	  function(x, ..., na.rm=F) {
	      if (length(x@rmat))
	           range(c(x@offdiag, x@blocks, x@rmat), na.rm=na.rm)
	      else
	           range(c(x@offdiag, x@blocks), na.rm=na.rm)
	      })

setMethod('any', signature('bdsmatrix'),
	  function(x, ..., na.rm=F) {
	      if (length(x@rmat))
	           any(c(x@offdiag, x@blocks, x@rmat), na.rm=na.rm)
	      else
	           any(c(x@offdiag, x@blocks), na.rm=na.rm)
	      })

setMethod('all', signature('bdsmatrix'),
	  function(x, ..., na.rm=F) {
	      if (length(x@rmat))
	           all(c(x@offdiag, x@blocks, x@rmat), na.rm=na.rm)
	      else
	           all(c(x@offdiag, x@blocks), na.rm=na.rm)
	      })

setMethod('sum', signature('bdsmatrix'),
	  function(x, ..., na.rm=F) {
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
			 integer(1), PACKAGE="kinship")$indexb     #index of diagonal elements
	      
	      n2 <- length(x@blocks)
	      nz <- d3^2 - sum(x@blocksize^2) #number of "offdiag" elements
	      wts <- rep(2, n2)
	      wts[temp] <- 1      # the diagonal elements
	      tsum <- sum(c(nz *x@offdiag, wts*x@blocks), na.rm=na.rm)

	      if (length(x@rmat)) {
		  wt2 <- rep(2, length(x@rmat))
		  wt2[row(x@rmat) > d3] <- 1
		  tsum <- tsum + sum(wt2*x@rmat, na.rm=na.rm)
		  }
	      tsum
	      })

setMethod('prod', signature('bdsmatrix'),
	  function(x, ..., na.rm=F) {
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
			 integer(1), PACKAGE="kinship")$indexb  #index of diagonal elements

	      n2 <- length(x@blocks)
	      nz <- d3^2 - sum(x@blocksize^2) #number of "offdiag" 
	      tprod <- 1
	      if (nz>0) {
		  if (x@offdiag==0) return(x@offdiag)
		  if (!is.na(x@offdiag) || na.rm==F)  tprod<- x@offdiag^nz
		  }

	      wts <- rep(2, n2)
	      wts[temp] <- 1      # the diagonal elements
	      tprod <- tprod * prod(x@blocks^wts, na.rm=na.rm)
	      if (length(x@rmat)) {
		  wt2 <- rep(2, length(x@rmat))
		  wt2[row(x@rmat) > d3] <- 1
		  tprod <- tprod * prod(x@rmat^wt2, na.rm=na.rm)
		  }
	      tprod
	      })


#
# For arithmetic operations, adding a single number preserves the structure
#  of the matrix, but adding a vector creates a matrix result which is
#  not block-diagonal.  Ditto for *, -, etc
#
setMethod('Ops', signature(e1='bdsmatrix', e2='numeric'),
	  function(e1, e2) {
	      if (length(e2)==1) {
		  e1@offdiag <- callGeneric(e1@offdiag, e2)
		  e1@blocks  <- callGeneric(e1@blocks,  e2)
		  if (length(e1@rmat))
			  e1@rmat    <- callGeneric(e1@rmat, e2)
		  e1
		  }
	      else {
		  bdsize <- .Options$bdsmatrixsize
		  if (is.null(bdsize)) bdsize <- 1000
		  if (prod(e1@.Dim) > bdsize)
		      stop("Automatic conversion would too large a matrix")
		  else callGeneric(as(e1, 'matrix'), e2)
		  }
	      } )

setMethod('Ops', signature(e1='numeric', e2='bdsmatrix'),
	  function(e1, e2) {
	      if (length(e1)==1) {
		  e2@offdiag <- callGeneric(e1, e2@offdiag)
		  e2@blocks  <- callGeneric(e1, e2@blocks)
		  if (length(e2@rmat))
			  e2@rmat    <- callGeneric(e1, e2@rmat)
		  e2
		  }
	      else {
		  bdsize <- .Options$bdsmatrixsize
		  if (is.null(bdsize)) bdsize <- 1000
		  if (prod(e2@.Dim) > bdsize)
		   stop("Automatic conversion would create too large a matrix")
		  else callGeneric(e1, as(e2, 'matrix'))
		  }
	      } )

setMethod('Ops', signature(e1='bdsmatrix', e2='bdsmatrix'), 
	  function(e1, e2) {
	      if (all(e1@.Dim == e2@.Dim) && 
		  (length(e1@blocksize) == length(e2@blocksize)) &&
		  all(e1@blocksize== e2@blocksize)) {
		  e1@offdiag <- callGeneric(e1@offdiag, e2@offdiag)
		  e1@blocks  <- callGeneric(e1@blocks,  e2@blocks)
		  if (length(e1@rmat))
			  e1@rmat <- callGeneric(e1@rmat, e2@rmat)
		  e1
		  }
	      else {
		  bdsize <- .Options$bdsmatrixsize
		  if (is.null(bdsize)) bdsize <- 1000
		  if (prod(e1@.Dim) > bdsize || prod(e2@.Dim) > bdsize)
		   stop("Automatic conversion would create too large a matrix")
	          else callGeneric(as(e1, 'matrix'), as(e2, 'matrix'))
		  }
	      })

setMethod('unique', c('bdsmatrix','missing'),
	  function(x, ...) unique(c(x@offdiag, x@blocks, x@rmat)))

bdsmatrix <- function(blocksize, blocks, rmat, dimnames=NULL) {
    nblock <- length(blocksize)
    if (any(blocksize <=0)) stop("Block sizes must be >0")
    if (any(as.integer(blocksize) != blocksize)) 
            stop("Block sizes must be integers")
    n1 <- sum(blocksize)
    n2 <- sum(blocksize^2)
    n3 <- sum(blocksize * (blocksize+1))/2
    if (length(blocks) == n2) {
	# Assume that they gave the full blocks, we only want the bottom
	#  half
        temp <- .C("bdsmatrix_index3",
		   as.integer(nblock),
		   as.integer(blocksize),
		   index=integer(n3), PACKAGE="kinship")$index
	blocks <- blocks[temp]
	}
    else if (length(blocks) != n3) 
	    stop("Length mismatch between blocks and blocksize")
    
    if (missing(rmat) || length(rmat)==0) {
	rmat <- numeric(0)
	n2 <- n1
	}
    else {
        rmat <- as.matrix(rmat)
	n2 <- n1 + ncol(rmat)           
	if (nrow(rmat) != n2) stop("Incompatible dimension for rmat")
	}

    if (!missing(dimnames) && !is.null(dimnames)) {
	if (is.list(dimnames) && length(dimnames)==2) {
	    if (length(dimnames[[1]])==0) val1 <- NULL
	    else { 
		val1 <- dimnames[[1]]
		if (length(val1) != n2) 
			stop("Invalid length for row dimnames")
		}
	    if (length(dimnames[[2]])==0) val2 <- NULL
	    else { 
		val2 <- dimnames[[2]]
		if (length(val2) != n2) 
			stop("Invalid length for column dimnames")
		}

	    dimnames <- list(val1, val2)
	    }
	else stop("dimnames must be a list of length 2")
	}

    new('bdsmatrix', .Dim=as.integer(c(n2,n2)), 
                      blocksize=as.integer(blocksize), 
                      blocks=as.numeric(blocks),
                      rmat=as.matrix(rmat), offdiag=as.numeric(0), 
                      .Dimnames=dimnames)
    }

setGeneric("bdsmatrix")
setMethod('[', 'bdsmatrix', 
 function(x, i, j, ..., drop=T) {
    if (class(x) != 'bdsmatrix') stop("Must be a bdsmatrix object")
    allargs <- list(...)
#   if (nDotArgs(...) !=2) stop("Two subscripts are required")
#   if (length(allargs)!=2) stop("Two subscripts are required")

    nblock <- length(x@blocksize)
    d <- x@.Dim
    d3 <- sum(x@blocksize)
    d4 <- length(x@blocks)
#   if (any(..1 > d[1]))
    if (any(i > d[1]))
        stop(paste("Array subscript (", max(i), 
		   ") out of bounds, should be at most ", d[1], sep=''))

#   if (any(..2 > d[2]))
    if (any(j > d[2]))
        stop(paste("Array subscript (", max(j), 
		   ") out of bounds, should be at most ", d[2], sep=''))

#   rows <- (1:d[1])[..1]
#   cols <- (1:d[2])[..2]
    rows <- (1:d[1])[i]
    cols <- (1:d[2])[j]
    
    brows <- rows[rows <= d3]  #the block-diagonal portion
    bcols <- cols[cols <= d3]
    brlen <- length(brows)
    bclen <- length(bcols)
    if (brlen>1 && (length(rows)==length(cols)) && all(rows==cols) &&
        all(diff(rows)>0)) {
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
	x		
	}

    else {
	# Now if brows==bcols, I would still have a bdmatrix object (the
	#   only asymmetry was in columns/rows of rmat),
	#   but the case is rare enough that I'm ignoring it.  Otherwise...
	# The result will not be block diagonal!
	if (brlen>0 && bclen>0) {
	    bdsize <- .Options$bdsmatrixsize
	    if (is.null(bdsize)) bdsize <- 1000
	    if (length(rows)*length(cols) > bdsize )
		  stop("Automatic conversion would create too large a matrix")
	    # I need to subscript the block diagonal portion
	    #  index2 is the rows() and cols() function for the block portion
	    temp <- .C('bdsmatrix_index2',
		       as.integer(nblock),
		       as.integer(x@blocksize),
		       rows= integer(length(x@blocks)),
		       cols= integer(length(x@blocks)),PACKAGE="kinship")
	    newmat <- matrix(x@offdiag, brlen, bclen)
	    rindex <- match(temp$rows, brows, nomatch=0)
	    cindex <- match(temp$cols, bcols, nomatch=0)
	    keep <- (rindex>0 & cindex >0)  #the row/col is one we want to keep
	    if (any(keep)) 
		newmat[cbind(rindex[keep], cindex[keep])] <- x@blocks[keep]

	    # the above has snatched and inserted all of the below the
	    #  diagonal parts.  For above diagonal, realize that I can just
	    #  swap the temp$rows, temp$cols for the 'upper trianglar'
	    #  stored version of blocks
	    rindex <- match(temp$cols, brows, nomatch=0) 
	    cindex <- match(temp$rows, bcols, nomatch=0)
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

## the following is similar to gchol and added by JH Zhao 14/03/04
## in fact they are already defined in diag.bdsmatrix.R so comment them
#setMethod('diag', signature(x='bdsmatrix'),
#    function(x,nrow=1,ncol=1) {
#        diag(as.matrix(x))
#        })

#setMethod('diag<-', 'bdsmatrix',
#function (x,value) {
#   if (class(x)!="bdsmatrix") stop("argument 1 should be a bdsmatrix!")
#   if (length(value)!=dim(x)[1]) stop("argument 2 has invalid length!")
#   blocksize <- x@blocksize
#   # upper left is blocks
#   ul <- sum(blocksize)
#   blocks <- x@blocks
#   nc <- 0
#   for(i in 1:length(blocksize))
#   { 
#      n <- blocksize[i]
#      z <- rep(0,n)
#      s <- 1
#      z[1] <- s
#      for(j in 1:(n-1))
#      {
#         s <- s + (n-j+1)
#         z[j+1] <- s
#      }
#      nall <- n*(n+1)/2
#      d <- blocks[(nc+1):(nc+nall)]
#      d[z] <- value[(nc+1):(nc+n)]
#      if (i==1) newblocks <- d
#      else newblocks <- c(newblocks,d)
#      nc <- nc + nall
#   }
#   if (length(x@rmat)>0) {
#      nn <- nrow(x@rmat)
#      # right lower is part of rmat
#      rl <- x@rmat[(ul+1):nn,]
#      diag(rl) <- value[(ul+1):nn]
#      rmat <- rbind(x@rmat[1:ul,],rl)
#   }
#   else rmat <- x@rmat
#   newmat <- bdsmatrix(blocksize=x@blocksize, blocks=newblocks,
#             rmat=rmat, dimnames=x@.Dimnames)
#   newmat
#})

is.matrix.bdsmatrix <- function(x, ...) is(x, "bdsmatrix")
is.matrix.gchol.bdsmatrix <- function(x, ...) is(x,"gchol.bdsmatrix")
is.list.bdsmatrix <- function(x, ...) !is(x, "bdsmatrix")
is.list.gchol.bdsmatrix <- function(x, ...) !is(x, "gchol.bdsmatrix")
is.bdsmatrix <- function(x) inherits(x, "bdsmatrix")
is.gchol.bdsmatrix <- function(x) inherits(x, "gchol.bdsmatrix")
