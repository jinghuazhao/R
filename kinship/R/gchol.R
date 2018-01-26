# $Id: gchol.s,v 1.2 2002/12/26 22:54:52 Therneau Exp $
#
# Code for the generalized cholesky  A = LDL', where L is lower triangular
#   with 1's on the diagonal, and D is diagonal.
# The decompostions exists for any square symmetric matrix.  
# If A is positive definite, then all elements of D will be positve.
# If A is not full rank, then 0's on the diagonal of D signal the redundant 
#   columns.  
#

#The original form adding setGeneric does not generate gchol.bdsmatrix class
#gchol <- function(x, tolerance=1e-10) UseMethod("gchol")
#package='kinship',valueClass=c('gchol.bdsmatrix',"gchol"))

# JH Zhao 27/3/2004
setGeneric("gchol", function(x,tolerance=1e-10) if(class(x) == 'matrix')
         standardGeneric("gchol") else gchol.bdsmatrix(x,tolerance))

setClass('gchol', 
	 representation(.Data= 'numeric',
			.Dim = 'integer',
			.Dimnames = 'list or NULL',
			rank = 'integer'))
			
as.matrix.gchol <- function(x, ones=T, ...) {
    temp <- matrix(x@.Data, x@.Dim[1], dimnames=x@.Dimnames, byrow=T)
    if (ones) diag(temp) <- 1
    temp
    }

setAs('gchol', 'matrix', function(from) as.matrix.gchol(from))

setMethod('gchol', signature(x='matrix'),
    function(x,  tolerance) {
	d <- dim(x)
	if (d[1] != d[2]) 
		stop("Cholesky decomposition requires a square matrix")
	if (!is.logical(all.equal(as.vector(x), as.vector(t(x)))))
		stop("Cholesky decomposition requires a symmetric matrix")
	temp <- .C("gchol", as.integer(d[1]),
		   x =   as.double(x),
		   rank= as.double(tolerance), PACKAGE="kinship")
	new('gchol', .Data= temp$x  , .Dim=as.integer(d), 
	    .Dimnames= dimnames(x), rank=as.integer(temp$rank))
	})

setMethod('diag', signature(x='gchol'),
    function(x,nrow=1,ncol=1) {
	d <- x@.Dim[1]
	x@.Data[ seq(1, length=d, by=d+1)]
	})

setMethod('show', 'gchol', function(object) show(as.matrix(object, FALSE)))
