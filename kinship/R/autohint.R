# $Id: autohint.s,v 1.4 2007/12/12 23:13:24 therneau Exp $ 
#
# Try to guess a decent hint matrix at the start of a pedigree
# The algorithm is pretty simple minded:
#   Get the plot order of the pedigree, and see if there are any "arcs"
#     joining spouses that are not next to each other
#   Move said spouses to be the leftmost/rightmost among their siblings
#   Make an entry in hints[,2] to capture the left/right side.
#   Iterate this for level 1, then 2, etc.
# For many pedigrees, founders need to be moved to get the "best" picture.
#   This routine isn't smart enough to do that.
#
autohint <- function(ped) {
    spousid <- NULL
    n <- length(ped$depth)
    #
    # First, find out if we have any twins
    #  if so, the vector "twinset" will be 0 for non-twins, and
    #  a unique integer for each unique set of twins.  "twinrel" is
    #  the relationships matrix, with only twin information retained
    #  twinord is the implied order of twins within a set
    #
    if (!is.null(ped$relation) && any(ped$relation[,3] <4)) {
	temp <- (ped$relation[,3] < 4)
	twinlist <- unique(c(ped$relation[temp,1:2]))  #list of twin id's 
	twinrel  <- ped$relation[temp,,drop=F]
	
	twinset <- rep(0,n)
	twinset[twinlist] <- 1:length(twinlist)  #give every twin a unique id
	twinord <- rep(1,n)
	for (i in 2:length(twinlist)) {
	    # Now, for any pair of twins on a line of twinrel, give both
	    #  of them the minimum of the two ids
	    # For a set of triplets, it might take two iterations for the
	    #  smallest of the 3 numbers to "march" across the threesome.
	    #  For quads, up to 3 iterations, for quints, up to 4, ....
	    newid <- pmin(twinset[twinrel[,1]], twinset[twinrel[,2]])
	    twinset[twinrel[,1]] <- newid
	    twinset[twinrel[,2]] <- newid
	    twinord[twinrel[,2]] <- pmax(twinord[twinrel[,2]], 
					 twinord[twinrel[,1]]+1)
	    }	
	}
    else {
	twinset <- rep(0,n)
	twinrel <- NULL
	}

    #
    # Here is an internal function which rearranges someone to be
    #  the leftmost or rightmost of his/her siblings.  The only
    #  real complication is twins -- if one of them moves
    #  the other has to move too.  And we need to keep the
    #  monozygotics together within a band of triplets....
    # Algorithm: first move all the twins to the left end (or right
    #  as the case may be), then move all the monozygotes to the
    #  left, then move the subject himself to the left.
    #
    shift <- function(id, sibs, goleft, hint, twinrel, twinset) {
	if (twinset[id]> 0)  { 
	    shift.amt <- 1 + diff(range(hint[sibs]))  # enough to avoid overlap
            twins <- sibs[twinset[sibs]==twinset[id]]
	    if (goleft) 
		 hint[twins] <- hint[twins] - shift.amt
	    else hint[twins] <- hint[twins] + shift.amt
		    
	    mono  <- any(twinrel[c(match(id, twinrel[,1], nomatch=0),
				   match(id, twinrel[,2], nomatch=0)),3]==1)
	    if (mono) {
		#
		# ok, we have to worry about keeping the monozygotics
		#  together within the set of twins.
		# first, decide who they are, by finding those monozygotic
                #  with me, then those monozygotic with the results of that
                #  iteration, then ....  If I were the leftmost, this could
                #  take (#twins -1) iterations to get us all
                #
		monoset <- id
		rel2 <- twinrel[twinrel[,3]==1, 1:2, drop=F]
		for (i in 2:length(twins)) {
		    newid1 <- rel2[match(monoset, rel2[,1], nomatch=0),2]
		    newid2 <- rel2[match(monoset, rel2[,2], nomatch=0),1]
		    monoset <- unique(c(monoset, newid1, newid2))
		    }
		if (goleft) 
		       hint[monoset]<- hint[monoset] - shift.amt
		else   hint[monoset]<- hint[monoset] + shift.amt
		}
	    }

	#finally, move the subject himself
	if (goleft) hint[id] <- min(hint[sibs]) -1   
	else	    hint[id] <- max(hint[sibs]) +1

        hint[sibs] <- rank(hint[sibs])  # aesthetics -- no negative hints
	hint
	}

    #
    # Now, get an ordering of the pedigree to "look at"
    #    we don't need the final "prettify" step, hence align=F
    # If there is a hints matrix entered, we retain it's non-zero entries,
    #   otherwise people are put into the order of the data set.  Twins are
    #   then further reordered
    if (is.null(ped$hints)) temp <- integer(n)
    else                    temp <- ped$hints[,1]
    
    for (i in unique(ped$depth)) {
	who <- (ped$depth==i & temp==0)
	if (any(who)) temp[who] <- 1:sum(who)
	}
    if (any(twinset>0)) {
	# First, make any set of twins a cluster: 6.01, 6.02, ...
	#  By using fractions, I don't have to worry about other sib's values
	for (i in unique(twinset)) {
	    if (i==0) next
	    who <- (twinset==i)
	    temp[who] <- mean(temp[who]) + twinord[who]/100
	    }

	# Then reset to integers
	for (i in unique(ped$depth)) {
	    who <- (ped$depth==i)
	    temp[who] <- rank(temp[who])
	    }
	}

    hints <- cbind(temp, 0)
    plist <- align.pedigree(ped, packed=T, align=F, hints=hints)


    #
    # Now, find the list of duplicate numbers on each level, and
    #   arrange the children/spouse hints to move them closer
    # By default, put the "marriage direction hint" onto the person who
    #   is connected to parents
    #
    maxlev <- nrow(plist$nid)
    for (lev in 1:maxlev) {
	idlist <- plist$nid[lev,1:plist$n[lev]] 
	dups <- duplicated(idlist)
	while (any(dups)) {
	    id <- (idlist[dups])[1] # the id to be fixed
	    dups[idlist==id] <- F   # avoid redoing this subject
	    xpos <- (1:plist$n[lev])[idlist==id] #positions at which he/she is
		
	    # Look at the first instance of the person.  My "doppleganger"
	    #  will be to the right, so make sure that I am listed as the
	    #  right hand part of any marriage, and of any set of sibs.  That
            #  way there is a chance that the two symbols for me will overlap.
	    # Is "id" the left or right part of this union?  
	    # Use that to get my spouse's id
	    x <- xpos[1]   # index, in the row, of this plot symbol
            if (plist$spouse[lev,x]) spouse.x <- x+1
            else                     spouse.x <- x-1
            spouseid <- idlist[spouse.x]

	    family <- plist$fam[lev,x]                
	    if (family > 0) {
		# This particular instance of me in the tree is the one that
                #  is attached to my parents (and sibs if I have them).  Make
		#  me the rightmost of the sibs.
		sibs <- match(plist$fam[lev,], family, nomatch=0)
		sibs <- idlist[sibs>0]
		hints[ ,1] <- shift(id, sibs, F, hints[,1], twinrel, twinset)

		# I have a dotted line connection to the right.  This means
		#  that I'm the spouse of someone further to the right in the
		#  pedigree (second listing of me is the redundant one).  The
		#  correct marking for the order of that marriage will be
		#  taken care of further below.
		# If I also have a solid line connection to the right this may
		#  be an insoluble case (both these spouses need to be to the
		#  immediate right).  Assume the simplifying case, i.e., that
                #  I displace the current "right spouse" to without
                #  harming things (unless there is already one to the left).
		# 
		if (plist$spouse[lev,x] && 
		          (x==1 || !plist$spouse[lev,x-1])){
		    hints[id,2] <- spouseid
		    hints[spouseid,2] <- 0   #don't let them hint at me
		    }
		}
	    else {
		# My family is somewhere to the right, but my spouse is here
		#  Make sure that I am listed as his/her "right" marriage
		#  If spouse has sibs at this point, make him/her rightmost
                # If both spouse slots are taken, then don't do anything.
                #   (You sometimes can't get multiple families correct all
                #   at the same time, especially with animals).
                if (hints[spouseid,2]== 0 || hints[id,2]== 0){
                    if (hints[spouseid,2] == 0) hints[spouseid,2] <-  id
                    else                        hints[id,2] <- -spouseid
		
                    family <- plist$fam[lev,spouse.x]
                    if (family > 0) {
                        sibs <- match(plist$fam[lev,], family, nomatch=0)
                        sibs <- idlist[sibs>0]
                        hints[,1] <- shift(spouseid, sibs, F, hints[,1],
				                twinrel, twinset)
                        }
                    }
		}

	    #
	    # Block 2
	    # Now care for the second instance of my symbol.  My "doppleganger"
	    #  will be to the left, so make sure that I am listed as the
	    #  left hand part of any marriage, and of any set of sibs.
            # 
	    x <- xpos[2]  
            if (plist$spouse[lev,x]) spouse.x <- x+1
            else                     spouse.x <- x-1
            spouseid <- idlist[spouse.x]
	    
            family <- plist$fam[lev,x]
	    if (family > 0) {
		# I have sibs at this spot
		sibs <- match(plist$fam[lev,], family, nomatch=0)
		sibs <- idlist[sibs>0]
		hints[,1] <- shift(id, sibs, T, hints[,1], twinrel, twinset)

		# I have a dotted line connection to the left.  If I also
		#  have a solid line connection, and no current right-hand
		#  spouse, we can displace that spouse.
		if (plist$spouse[lev, x-1] && !plist$spouse[lev,x]) {
		    hints[id,2] <- spouseid
		    hints[spouseid, 2] <- 0
		    }
		}
	    else {
		# Must be a marriage at this point, make me leftmost
		#  If spouse has sibs at this point, make him/her leftmost
                # If both hint slots are taken, or these two subjects are 
                #  already mapped to each other, leave things alone.
                if ((hints[id]==0 || hints[spouseid]==0) && 
                        (abs(hints[id]) != spouseid)  &&
                        (abs(hints[spousid]) != id)) {
                    if (hints[spouseid,2] == 0) hints[spouseid,2] <- -id
                    else                        hints[id,2] <- spouseid
		
                    family <- plist$fam[lev,spouse.x]
                    if (family > 0) {
                        sibs <- match(plist$fam[lev,], family, nomatch=0)
                        sibs <- idlist[sibs>0]
                        hints[,1] <- shift(spouseid, sibs, T, hints[,1],
				                twinrel, twinset)
                        }
                    }
		}
            # Block 3 - take care of the third, fourth, etc instances of "id"
            #   on the same line.  
            # No block 3 code has been written!  Why?  This case only tends to
            #   occur in very complex families, e.g., animals, and in pedigrees
            #   with that level of interconnect the resultant plot is going to
            #   look like spagetti no matter what we do.
	    }
        #
        # Recompute, since this shifts things on levels below
        #
	plist <- align.pedigree(ped, packed=T, align=F, hints=hints)
	}
    dimnames(hints) <- list(NULL, NULL)

    # Cleanup
    check.hint(hints, ped$sex)
}
