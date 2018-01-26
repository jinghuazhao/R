#  $Id: alignped4.s,v 1.1 2003/02/05 00:08:45 atkinson Exp $ 
# Do what we can to line it up, sliding unfilled rows.
# The algorithm is to try to mimic a physical system: imagine that
#   each "block" of the pedigree (children and their spouses) was on
#   rollers and could move left-right, the whole pedigree was in a box
#   of width "width", and the connections from parent to children were
#   each a spring.
# This version of the algorithm has 2 steps: first get a "perfect"
#   lineup by making each row fixed in turn, and moving it's neighbors
#   to line up with it.  Then, push in on the box.
#
alignped4 <- function(rval, spouse, level, iter, width=2) {
    n <- length(level)
    maxlev <- max(level)

    pos <- rval$pos
    width <- max(width-1, rval$n)   #actually, the max x-coord

    # Make the indicator of "blocks", start with sibs, i.e., rval$fam, 
    #    then add spouses
    # Then, give a fake famid to each of the founder families
    fam2   <- rval$fam      # to save typing 'rval$fam' all the time
    blocks <- rval$fam
    for (i in 1:maxlev) {
	if (any(spouse[i,])) {
	    has.spouse <- (1:ncol(spouse))[spouse[i,]]
	    for (j in 1:3) {  # A max of 3 spouses to either side
		temp <- pmax(blocks[i, has.spouse], blocks[i, has.spouse+1])
		if (any(temp==0)) { 
		    # founder marriages
		    temp[temp==0] <- n + 1:(sum(temp==0))
		    }
		blocks[i, has.spouse] <- blocks[i, has.spouse+1] <- temp
		}
	    
	    # Last step -- some families may have merged 
	    #  (both husband and wife have sibs): all need to move as
	    #  a block. Key off people who "changed" their family id 
	    changed.fam <- (fam2[i,]>0 & (fam2[i,] != blocks[i,]))
	    if (any(changed.fam)) {
		oldfam <- fam2[i, changed.fam]
		newfam <- blocks[i, changed.fam]
		for (j in (length(newfam)):1) { # go high to low
		    blocks[i, blocks[i,]== oldfam[j]] <- newfam[j]
		    }
		}
	    }
	}
    # 
    # At this point block contains the grouping of the subjects
    # For instance, the 3rd row might be 2,4,4,7,15,15,0 
    #   This means:  6 plotted symbols at this level (0's are a filler)
    #                the widest level has 7 plotted symbols
    #                objects 2 and 3 must stick together, as must 6 and 7
    #                parents of the first subject are positions 2 and 3 of
    #          blocks[2,], the parents of one or both of 2 and 3 are
    #          subjects 4 and 5 in blocks[2,], of the third are cols 7&8
    #		 The subjects plotted at postitions 7 and 8 have no
    #          parents in the pedigree.
    #  If fam2[3,] were 2, 4,4,... then subjects 2-3 are sibs, if it were
    #                   2, 0,4,... then subject 3 is the child, 2 a spouse
    #  The vector pull[] will contain the computed left-right shift
    #    for each block of subjects, i.e., 4 elements in row 3 for the
    #    sample blocks[] given above
    #

    pfun <- function(blocks, fam, pos) {
	# Given a subset (2 rows) of the blocks, fam and pos matrices,
	#  return a 3 column matrix with one row for each connection ("spring")
	#  and colums of top row attachment (block number), bottom row
	#  attachment block, and (top row pos) - (bottom row pos)
	# The lower attachment point of each spring is kept in "centers",
	#  and is the mean of the children's positions (but not spouses of
	#  children.  The upper attachment point is centered between spouses.
	ff <- sort(unique(fam[2,]))  # the unique family id's in lower row
	centers <- tapply(pos[2,], fam[2,], 
			      function(x) mean(range(x))) 
	if (any(ff==0)) {  # the 0's were spouses of kids, not a family of 0's
	    ff <- ff[-1]
	    centers <- centers[-1]
	    }
	pcenter <- (pos[1,ff] + pos[1,ff+1]) /2 #upper attachment point

	blockid <- blocks[2, match(ff, fam[2,])]
	res <- cbind(blocks[1,ff], blockid, pcenter - centers)
	dimnames(res) <- NULL
	res
	}

    # the "touching" function
    tfun <- function(pos, pull, ns, block, width) {
	# Adjust the pulls so that no one overlaps, and no one sticks out
	#  past the edges.  If width==0, ignore the edge constraint
	fixed <- 0*pos       # =1 for cells that can't move any more
	while (T) {
	    nopush <- T      #if this is still true at the bottom, we're done
	    pos2 <- pos + pull
	    if (width>0) {   # width==0 means "don't check width"
		if (min(pos2) <0) {
		    nopush <- F
		    i <- (unique(block[pos2<0]))[1]  
		    # The [1] above: I'm not sure that it's possible to have
		    #  two blocks hanging over the edge at once -- but maybe
		    who <- (block==i)
		    pull[who] <- pull[who] - min(pos2)
		    fixed[who] <- 1
		    }
		if (max(pos2) > width) {
		    nopush <- F
		    i <- (rev(unique(block[pos2 > width])))[1] 
		    #see [1] comment above
		    # if there were 2, I want the last one in the list
		    who <- (block==i)
		    pull[who] <- pull[who] - (max(pos2) - width)
		    fixed[who] <- 1
		    }
		pos2 <- pos + pull
		}

	    space<- diff(pos2)  
	    rspace <- round(space,2)  #stop infinite loops due to roundoff
	    if (length(pos2)>1 && any(rspace <1)) {
		# line above: round() stops infinite loops due to round-off
		nopush <- F
	        # Because multiple blocks may have hit, I need to do this
	        #  one at at time.  Blocks get "merged" so that ones that
	        #  have hit slide all together from that point on.
	        # Suppose that there are lots of overlaps.  Then the leftmost
	        #  block might get pushed and pushed and pushed and...
	        #  So it is necessary to use "nspring" as a weight
		i <- (1:length(space))[rspace<1] #the left member of collisions
		i <- i[1]
		lblock <- block[i];  lwho <- (block==lblock)
		rblock <- block[i+1];rwho <- (block==rblock)

		if (any(fixed[lwho]==1)) {
		    # left one can't be  moved
		    if (any(fixed[rwho]==1)) stop("Logic bug 1")
		    pull[rwho] <- pull[rwho] + (1 - space[i])
		    block[rwho] <- lblock
		    }
		else if (any(fixed[rwho]==1)) {
		    # right one can't move
		    pull[lwho] <- pull[lwho] - (1-space[i])
		    block[lwho] <- rblock
		    }
		else {
		    nspring <- ns[i] + ns[i+1]
		    slide <- (1-space[i])/nspring
		    pull[lwho] <- pull[lwho] - slide*ns[i+1]
		    pull[rwho] <- pull[rwho] + slide*ns[i]
		    block[rwho] <- lblock;
		    ns[rwho | lwho] <- nspring;
		    }
		}

	    if (nopush) break
	    }
	as.vector(pull)  #no names
	}

    #
    # Do k iterations of the algorithm, with the first pass ignoring
    #   boundary constraints.  
    # At iteration 3, level 1 is only adjusted from bottom up, and level
    #   maxlev from the top down.  At iter 4, levels 1-2 only get adjusted
    #   from the bottom, and etc.  Thus, there are at most 2 + floor(maxlev/2)
    #   iterations (at most 2 if maxlev==2)
    #
    maxiter <- 2 + floor(maxlev/2)
    if (maxlev==2) maxiter<-2
    if (is.numeric(iter)) icount <- min(max(1, round(iter)), maxiter)
    else icount <- maxiter
    for (iter in 1:icount) {
	realval <- (col(pos) <= rval$n[row(pos)])     # not structural zeros
	if (max(pos[realval]) > width) {
	    # try to pre-center so that the smallest number of blocks 
	    #   sits outside the boundaries
	    temp1 <- sort(unique(pos))
	    temp2 <- table(pos[realval])
	    temp3 <- temp1
	    for (i in 1:length(temp3)) {
		lower.limit <- temp1[i]
###############
##Begin patch##
###############
tmp_var1 <- (temp1 < lower.limit)
tmp_var2 <- (temp1 > (lower.limit + width))
        if( length(tmp_var1) > 0 )
        {
        tmp_summand_1 <- sum(temp2[tmp_var1], na.rm=TRUE)
        }
        else
        {
        tmp_summand_1 <- 0
        }
        if( length( tmp_var2) > 0 )
        {
        tmp_summand_2 <- sum(temp2[tmp_var2], na.rm=TRUE)
        }
        else
        {
        tmp_summand_2 <- 0
        }
                        #temp3[i] <- sum(temp2[temp1 < lower.limit]) + sum(temp2[temp1 > (lower.limit + width)])
                        temp3[i] <- tmp_summand_1 + tmp_summand_2
#############
##End Patch##
#############
		}
	    pos[realval] <- pos[realval] - mean(temp1[temp3==min(temp3)])
	    }

	# Go from the top down, treating row lev as fixed and moving row 
	#   'lev+1'.
	for (lev in max(1,iter-1):(maxlev-1)) {
	    temp <- pfun(blocks[lev+0:1,], fam2[lev +0:1,], pos[lev+ 0:1,])

	    utemp <- sort(unique(temp[,2]))
	    pull <- tapply(temp[,3], temp[,2], mean)  #one number per block
	    pull <- pull[match(blocks[lev+1,], utemp)] #one number per subject
	    nspring <- tapply(temp[,3], temp[,2], length)
	    nspring <- nspring[match(blocks[lev+1,], utemp)] 
	    np <- is.na(pull)   # those with no parents have no attachments
	    pull[np]   <- 0
	    nspring[np]<- 0

	    # Center the child under the parent 
	    nn <- rval$n[lev+1]
	    pull2 <- tfun(pos[lev+1, 1:nn], pull[1:nn], nspring[1:nn], 
			  blocks[lev+1, 1:nn], width= width*(iter>1)) 
	    pos[lev+1, 1:nn] <- pos[lev+1, 1:nn] + pull2
	    }

	# Now repeat the process from bottom up, centering parents above
	# their children
	# In this case the pull vector is for the parents, above it was 
	# for kids
	for (lev in min(maxlev-1, 1+maxlev-iter):1) {
	    temp <- pfun(blocks[lev+0:1,], fam2[lev +0:1,], pos[lev+ 0:1,])

	    utemp <- sort(unique(temp[,1]))
	    pull <- tapply(-temp[,3], temp[,1], mean)  #one number per block
	    pull <- pull[match(blocks[lev, ], utemp)]  #one number per subject
	    nspring <- tapply(temp[,3], temp[,1], length)
 	    nspring <- nspring[match(blocks[lev,], utemp)] 
	    np <- is.na(pull)   # those with no children have no attachments
	    pull[np]   <- 0
	    nspring[np]<- 0

	    # Ensure no overlaps
	    nn <- rval$n[lev]
	    pull2 <- tfun(pos[lev, 1:nn], pull[1:nn], nspring[1:nn], 
			  blocks[lev, 1:nn], width= width*(iter>1)) 
	    pos[lev, 1:nn] <- pos[lev, 1:nn] + pull2
	    }

	# When subtracting a global constant, we want to leave the "structural
	#  zeros" (values of the pos matrix that correspond to no one) alone.
	#
	pos[realval] <- pos[realval] - min(pos[realval])  # set left edge =0
	}

    pos
    }
