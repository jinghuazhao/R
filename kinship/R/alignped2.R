#  $Id: alignped2.s,v 1.1 2003/02/05 00:08:44 atkinson Exp $

# The real work is in three routines:
#   alignped1: called with a single subject, returns the subtree founded
#      on this subject, as though it were the only subtree.
#   alignped2:  called with a set of sibs: calls align multiple times and
#      merges the result
#   alignped3: glue two results together
#

# Input: x -- the list of sibs
#        id, dad, mom, level, horder, sorder: as defined at start of alignped
#	 packed: what result do I want
# Return: a return list for the subtree, of the form of align's final one
#
# Pedigrees are small: efficiency isn't the issue, rather ease of writing

alignped2 <- function(x, id, dad, mom, level, horder, sorder, packed) {
    x <- x[order(horder[x])]  # Use the hints
    rval <- alignped1(x[1], id, dad, mom, level, horder, sorder, packed) 

    if (length(x) >1) {
	mylev <- level[x[1]]
	for (i in 2:length(x)) {
	    rval2 <-  alignped1(x[i], id, dad, mom, level,
				horder, sorder, packed)
	    
	    # Deal with a very unusual special case:
	    #   Subject x married one of his/her sibs, so has already
	    #       appeared somewhere to the right (in rval)
	    #   This instance added no spouses/kids  (rval2$n =1)
	    #   So we don't want to append this second instance of x
	    # In all other cases (99.99%) we do want to append
	    if ((rval2$n[mylev] > 1) ||
		          (is.na(match(x[i], floor(rval$nid[mylev,])))))
		rval <- alignped3(rval, rval2, packed)
	    }
	}
    rval
    }
