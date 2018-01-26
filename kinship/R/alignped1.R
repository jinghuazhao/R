# $Id: alignped1.s,v 1.1 2003/02/05 00:08:43 atkinson Exp $ 
# The real work is in three routines:
#   alignped1: called with a single subject, returns the subtree founded
#      on this subject, as though it were the only subtree.
#   alignped2:  called with a set of sibs: calls alignped1 multiple times and
#      merges the result
#   alignped3: glue two results together
#
# In these routines, the nid array = "final nid" + (1/2)"final spouse array"
#  align.pedigree takes them back apart
#
#
# Input: x -- the nid of the subject in question
#        id, dad, mom, level, horder, sorder: as defined at start of alignped
#	 packed: what result do I want
# Return: a return list for the subtree, of the form of align's final one

alignped1 <- function(x, id, dad, mom, level, horder, sorder, packed) {
    # Set a few constants
    maxlev <- max(level)
    lev <- level[x]
    n <- integer(maxlev)

    # Get list of not-yet-processed spouses
    spouselist <- spouselist       # copy this from frame 0
    if (length(spouselist)==0)  spouse <- NULL
    else {
	if (any(spouselist[1,]==x)){
	    sex <- 1		                  # I'm female
	    which.col <- (spouselist[1,]==x)
	    spouse <- spouselist[2, which.col]
	    spouse <- spouse[level[spouse]<= lev]  # cross-level marriages are
	                                           # plotted on the lower level
	    }
	else {
	    sex <- 2
	    which.col <- (spouselist[2,]==x)   #could be none at all
	    if (any(which.col)) {
		spouse <- spouselist[1, which.col]
		spouse <- spouse[level[spouse] <= lev]  
		}
	    else spouse <- NULL
	    }
	}
    nspouse <- length(spouse)  # Almost always 0, 1 or 2

    # Create the return structures, of a width (1 + nspouse), which is
    #  the minimum necessary.  If there are wider levels below due to
    #  children, they will be widened below
    nid <- fam <- matrix(0, maxlev, nspouse+1)
    pos <- matrix(0.0, maxlev, nspouse +1)
    n[lev] <- nspouse +1       
    pos[lev,] <- 0:nspouse
    if (nspouse ==0) {          #No children, or already processed them
        # Easy case: the "tree rooted at x" is --- just x
	nid[lev,1] <- x
	return(list(nid=nid, pos=pos, fam=fam, n=n))
	}

    # Ok, we have a list of spouses.  Now, look at the spouse hints (sorder)
    #  to see if any of them are specified as being plotted to the left
    #  or to the right. 
    lspouse <- rspouse <- NULL
    temp <- (abs(sorder) ==x)
    if (any(temp)) { # I am someone else's hint
        who <- (1:length(temp))[temp]
        indx <- match(who, spouse, nomatch=0)  #will only be zero with
        who <- who[indx>0]                     # certain multiple marriages
                                               # (one spouse already processed)
        if (length(who)) {
            lspouse <- (spouse[indx])[sorder[who]>0]
            rspouse <- (spouse[indx])[sorder[who]<0]
            }
        }
    if (sorder[x]!=0) {
        indx <- match(abs(sorder[x]), spouse)  #which spouse is in the hint
	if (is.na(indx)) {
	    # This can happen -- there are multiple spouses, the one mentioned
	    #  in the hint is someone's child, and this marriage has 
	    #   already been processed and thus removed from the spouselist.
	    # There is already one marriage to the left, so put more of these
	    #  to the right
            bias <- .5  #the variable not yet used -- till some smarter version
            }
	else {
	    if (sorder[x] <0)   #hinted spouse to the left
                 lspouse <- unique(c(lspouse, spouse[indx]))
	    else rspouse <- unique(c(rspouse, spouse[indx]))
            # unique() is needed because of double hints: spouse 2 shows #7
            #  on thier right, #7 shows #2 on thier left, giving repeats in
            #  the lspouse or rspouse vector
	    }  
	}
    #
    # At this point some future routine may put in a little more smarts:
    #  if #spouses >2, one of them has parents, and those parents are to the
    #  right (not yet processed): perhaps we should not processes that union,
    #  leave them on the spouselist, so that the marriage and children are
    #  drawn further to the right, where the pedigree is less crowded.
    # (Yes, if I have 3 spouses, but my #3 has been wed more times yet, we are
    #  better off here.  Will this actually happen?)
    #
    
    # If there are more spouses still, divide them to the left & right.
    # For 3 spouses, we put 2 on the right, 1 on the left (perhaps this
    #  person already has one more, dotted line spouse, to the left!).
    ntemp <- length(c(lspouse, rspouse))
    if (ntemp==0) {
	# There is no spouse hint: decide how many are plotted to the left
        # Unhinted spouses print in the order found in spouselist.
	if (nspouse==1) {
	    #plot male on the left by default
	    if (sex==1) nid[lev,1:2] <- c(spouse+.5, x) 
	    else        nid[lev,1:2] <- c(x+.5, spouse)
	    }
	else {
	    nleft <- floor(nspouse/2)  #number plotted to the left
	    temp <- c(spouse[1:nleft], x, spouse[(nleft+1):nspouse])
	    nid[lev, ] <- temp + c(rep(.5, nspouse), 0)
	    }
	}
    else {
        # Put them in the printing array, with lspouse closest to x on the
        #   left, rspouse closest on the right, and any other ones split 
        #   between "farther left" and "farther right".
        # At the end, put spouselist into this order, so as to properly
        #   attach children.  (A lot of work for the few triple+ marriages).
        if (nspouse > ntemp) {  #more spouses to do!
            tspouse <- spouse[spouse!=c(lspouse, rspouse)]  # those remaining
            nleft <- floor(nspouse/2) - length(lspouse)
            temp <- rep(T, length(tspouse))
            if (nleft>0) temp[1:nleft] <- F
            nid[lev,] <- c(tspouse[!temp], lspouse, x, rspouse, tspouse[temp])+
                         c(rep(.5, nspouse), 0)
            spouse <- c(tspouse[!temp], lspouse, rspouse, tspouse[temp])
            }
        else { #only the left and right ones, as determined by hints
            nid[lev,] <- c(lspouse, x, rspouse) +  c(rep(.5, nspouse), 0)
            spouse <- c(lspouse, rspouse)
            }
	}

    # Remove the ones we'll take care of from the spouselist
    toss <- (spouselist[sex,]==x  & !is.na(match(spouselist[3-sex,], spouse))) 
    assign('spouselist', spouselist[,!toss, drop=F], envir=as.environment(sys.frame()))
	
    # Now to work  
    # For each spouse, find the children.  If there are any, call alignped2
    #  to create their tree, then mark them with this parent, and call
    #  alignped3 to join the trees.
    nokids <- T   #haven't found any kids yet
    for (i in 1:nspouse) {
	ispouse <- spouse[i]
	children <- id[(dad==x | mom==x) & (dad==ispouse | mom==ispouse)]
	if (length(children) > 0) {
	    children <- children[order(horder[children])]
	    rval1 <- alignped2(children, id, dad, mom, level, horder, sorder, 
			      packed)
	    # set the parentage for any kids
	    #  a nuisance: it's possible to have a child appear twice, when
	    #  via inbreeding two children marry --- makes the "indx" line
	    #  below more complicated
	    temp <- floor(rval1$nid[lev+1,])  # cut off the .5's for matching
	    indx <- (1:length(temp))[match(temp,children, nomatch=0) >0]
	    rval1$fam[lev+1,indx] <- i   #set the kids parentage
	    if (!packed) {
		# line the kids up below the parents
		# The advantage at this point: we know that there is 
		#   nothing to the right that has to be cared for
		kidmean <- mean(rval1$pos[lev+1, indx])
		parmean <- mean(pos[lev, i + 0:1])
		if (kidmean > parmean) {
		    # kids to the right of parents: move the parents
		    indx <- i:(nspouse+1)
		    pos[lev, indx] <- pos[lev, indx] + (kidmean - parmean)
		    }
		else {
		    # move the kids and thier spouses and all below
		    shift <- parmean - kidmean
		    for (j in (lev+1):maxlev) {
			jn <- rval1$n[j]
			if (jn>0) 
			    rval1$pos[j, 1:jn] <- rval1$pos[j, 1:jn] +shift
			}
		    }
		}
	    if (nokids) {
		rval <- rval1
		nokids <- F
		}
	    else {
		rval <- alignped3(rval, rval1, packed)
		}
	    }
	}
    #
    # Add x and the spouses onto the pedigree just created
    #
    if (nokids) {
	return(list(nid=nid, pos=pos, fam=fam, n=n))
	}

    if (max(rval$n) >= 1+nspouse) {
	# The rval list has room for me!
	rval$n[lev] <- n[lev]
	indx <- 1:(nspouse+1)
	rval$nid[lev, indx] <- nid[lev,]
	rval$pos[lev, indx] <- pos[lev,]
	}
    else {
	#my structure has room for them
	indx <- 1:(max(rval$n))
	rows <- (lev+1):maxlev
	n[rows] <- rval$n[rows]
	nid[rows,indx] <- rval$nid[rows,]
	pos[rows,indx] <- rval$pos[rows,]
	fam[rows,indx] <- rval$fam[rows,]
	rval <- list(nid=nid, pos=pos, fam=fam, n=n)
	}
    rval
    }
