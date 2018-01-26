# $Id: align.pedigree.s,v 1.3 2007/12/12 23:13:23 therneau Exp $ 
#
# This routine returns a matrix of positions for either a packed or
#   spread out pedigree.
# The structure it builds is 1. a set of matrices, each of dimension
#   maxdepth rows and maxwidth columns.  
#	  nid: left to right, the index number of the subject to be plotted
#	  pos: the x position of said subject
#         fam: the index of the attachment point, upward.  Ties in this
#	       row demark the sibships.  If the people in col 5&6 in the
#              level above are the parents, then all sibs have "5" (we index
#	       to the leftmost parent)
#         spouse: If True, the person to my right is my spouse
#
#  and 2. a vector "n" giving the number of subjects per row
#
# The max possible in a row, worst case of arrangement, is the number
#   of subjects at that level + # marriages at that level.      
#

# Arguments: ped     a pedigree
#	     packed  return a packed or "compressed" view
#            align   extra effort to align parents and children
#            hints   helps in laying out the order
#	     width   for a packed pedigree, the minimum width we have to
#			work with in the realignment stage
#
#  Final width of a packed pedigree is max(width, max(subjects on any one row))
#
align.pedigree <- function(ped, packed=T, hints=ped$hints, width=6, align=T) {
    n <- length(ped$depth)
    if (ncol(hints) !=2 || nrow(hints) != n) stop("Invalid hints matrix")
    hints <- check.hint(hints, ped$sex)
    
    # all the work will be done with "row number" as the id variable
    nid <- 1:n
    dad<- match(ped$dadid, ped$id, nomatch=0)
    mom<- match(ped$momid, ped$id, nomatch=0) 
    level <- ped$depth +1

    horder <- hints[,1]   # relative order of siblings within a family
    sorder <- hints[,2]   # -x = id of left spouse, +x = id of right

    # Currently, the alignment routine requires that you have either
    #  0 parents or 2 parents, not 1.
    if (any(dad==0 & mom>0) || any(dad>0 & mom==0))
	    stop("Everyone must have 0 parents or 2 parents, not just one")

    #
    # Compute the list of spouses for each subject
    #   Mothers are row 1, fathers row 2, each column a spouse pair
    #   Either having children, having a hint, or a spousal
    #   relationship gets you on the list
    # Start with a hashlist: fatherid * max(id) + mother
    hash <- (dad*n + mom)[mom>0 & dad>0]  #those with children
    temp1 <- abs(sorder)
    temp1 <- ifelse(ped$sex=='male', (1:n)*n + temp1, 1:n + n*temp1)
    hash <- c(hash, temp1[sorder!=0]) #those with a spouse hint
    if (!is.null(ped$relation) && any(ped$relation[,3]==4)) {
	who <- (ped$relation[,3]==4)  # add spouses from relationship list
	indx <- ped$relation[who,1]   # id of the first spouse
	temp1 <- ifelse(ped$sex[indx]=='male', n*indx + ped$relation[who,2],
			                       indx + n*ped$relation[who,2])
	hash <- c(temp1, hash) #being first is important -- it controls plot
	                       # order per the documentation
	}

    hash <- hash[!duplicated(hash)]   #eliminate duplicates
    hash <- hash -1                   #change to range of 0 to n-1 (for %%)
    spouselist <- rbind(1+ hash%%n, floor(hash/n))  # uncompress it
    assign('spouselist', spouselist, envir=as.environment(sys.frame()))

    #
    # Process the founders, 1 by 1 (choose females, with males to the
    #   left as our default hint).  
    # Those with no kids aren't shown in the plot, so skip em
    # To be a founder, neither the wife nor the husband can have parents
    #   (you must have both a dad and a mom, or neither)

    noparents <- (dad[spouselist[1,]]==0 & dad[spouselist[2,]]==0)
    founders <-  unique(spouselist[1,noparents])
    founders <-  founders[order(horder[founders])]  #use the hints
    rval <- alignped1(founders[1], nid, dad, mom, level, horder, sorder,
		              packed=packed)

    if (length(founders)>1) {
	for (i in 2:length(founders)) {
	    rval <- alignped3(rval, alignped1(founders[i], nid, dad, mom,
					      level, horder, sorder, packed), 
				      packed)
	    }
	}

    #
    # Unhash out the spouse and nid arrays
    #
    nid    <- floor(rval$nid)
    spouse <- 1*(rval$nid != nid)
    maxdepth <- nrow(nid)
    # For each spouse pair, find out if it should be connected with
    #  a double line.  This is the case if they have a common ancestor
    ancestor <- function(me, momid, dadid) {
	alist <- me
	repeat {
	    newlist <- c(alist, momid[alist], dadid[alist])
	    newlist <- sort(unique(newlist[newlist>0]))
	    if (length(newlist)==length(alist)) break
	    alist <- newlist
	    }
	alist[alist!=me]
	}
    for (i in (1:length(spouse))[spouse>0]) {
        a1 <- ancestor(nid[i], mom, dad)
        a2 <- ancestor(nid[i+maxdepth],mom, dad)  #matrices are in column order
	if (any(duplicated(c(a1, a2)))) spouse[i] <- 2
	}


    # Create the twins array, if neccessary
    #  Its values are: 0 = nothing, 1= the sib to my right is a
    #   monzygotic twin, 2= the sib to my right is a dizygote,
    #   3= the sib to my right is a twin, unknown zyogosity.
    #
    if (!is.null(ped$relation) && any(ped$relation[,3] < 4)) {
	twins <- 0* nid
	who  <- (ped$relation[,3] <4)
	ltwin <- ped$relation[who,1]
	rtwin <- ped$relation[who,2]
	ttype <- ped$relation[who,3]
	
	# find where each of them is plotted (any twin only appears
	#   once with a family id, i.e., under their parents)
        ntemp <- ifelse(rval$fam>0, nid,0) # matix of connected-to-parent ids
	ltemp <- (1:length(ntemp))[match(ltwin, ntemp, nomatch=0)]
	rtemp <- (1:length(ntemp))[match(rtwin, ntemp, nomatch=0)]
	twins[pmin(ltemp, rtemp)] <- ttype
	}
    else twins <- NULL

    # At this point the pedigree has been arranged, with the positions
    #  in each row going from 1 to (# subjects in the row).
    # Now pretty up the alignment (for plotting).  Trees that are
    #  not packed need only one iteration of clean up.
    # Note that alignped4 wants a T/F spouse matrix: it doesn't care
    #  about your degree of relationship to the spouse.
    if (align && max(level) >1) {
	if (packed || is.numeric(align))
		pos <- alignped4(rval, spouse>0, level, align, width)
	else    pos <- alignped4(rval, spouse>0, level, 1, width)
	}
    else pos <- rval$pos

    remove('spouselist', envir = as.environment(sys.frame()))
    if (is.null(twins))
         list(n=rval$n, nid=nid, pos=pos, fam=rval$fam, spouse=spouse)
    else list(n=rval$n, nid=nid, pos=pos, fam=rval$fam, spouse=spouse, 
              twins=twins)
    }

