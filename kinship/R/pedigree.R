#  $Id: pedigree.s,v 1.2 2003/10/22 18:14:31 Atkinson Exp $ 
pedigree <- function(id, dadid, momid, sex, affected, status, relations) {
    # Do some data checks
    n <- length(id)
    if (length(momid) != n) stop("Mismatched lengths, id and momid")
    if (length(dadid) != n) stop("Mismatched lengths, id and momid")
    if (length(sex  ) != n) stop("Mismatched lengths, id and sex")

    # Don't allow missing id values
    if (any(is.na(id))) stop("Missing value for the id variable")
    if (!is.numeric(id) && any(as.character(id) == ''))
	stop("Empty string is not allowed as the id variable")

    # Allow for character/numeric/factor in the sex variable
    if(is.factor(sex))
	    sex <- as.character(sex)
    codes <- c("male","female", "unknown", "terminated")
    if(is.character(sex)) sex<- charmatch(casefold(sex, upper = F), codes, 
					  nomatch = 3)	

    # assume either 0/1/2/4 =  female/male/unknown/term, or 1/2/3/4
    #  if only 1/2 assume no unknowns
    if(min(sex) == 0)
	    sex <- sex + 1
    sex <- ifelse(sex < 1 | sex > 4, 3, sex)
    if(all(sex > 2))
	    stop("Invalid values for 'sex'")
	else if(mean(sex == 3) > 0.25)
		warning("More than 25% of the gender values are 'unknown'")
    sex <- factor(sex, 1:4, labels = codes)

    father <- match(dadid, id, nomatch = 0)
    if(any(sex[father] != "male")) {
	who <- unique((id[father])[sex[father] != "male"])
	stop(paste("Id not male, but is a father:", paste(who, collapse = " ")))
	}

    mother <- match(momid, id, nomatch = 0)
    if(any(sex[mother] != "female")) {
	who <- unique((id[mother])[sex[mother] != "female"])
	stop(paste("Id not female, but is a mother:", paste(who, collapse = " ")))
	}

    # Paste the data together into the skeleton of a pedigree
    depth <- kindepth(id, momid, dadid, align = T)
    temp <- list(id = id, momid = momid, dadid = dadid, sex = sex, depth = depth)

    # More data checks -- do these if affected/status/relations are filled in
    #   (the optional parts of a pedigree structure)
    if (!missing(affected)) {
	if (is.matrix(affected)){
	    if (nrow(affected) != n) stop("Wrong number of rows in affected")
	    if (is.logical(affected)) affected <- 1* affected
	    } 
	else {
	    if (length(affected) != n)
		stop("Wrong length for affected")

	    if (is.logical(affected)) affected <- as.numeric(affected)
	    if (is.factor(affected))  affected <- as.numeric(affected) - 1
	    }
	if(max(affected) > min(affected)) affected <- affected - min(affected)
        else if (max(affected) == min(affected)) affected <- affected - 1
	if (!all(affected==0 | affected==1 | affected==2))
		stop("Invalid code for affected status")
	temp$affected <- affected
	}

    if(!missing(status)) {
	if(length(status) != n)
		stop("Wrong length for status")
	if(any(status != 0 & status != 1))
		stop("Invalid status code")
	temp$status <- status
	}

    if (!missing(relations)) {
	if (is.matrix(relations)) {
	    if (ncol(relations) != 3) 
		    stop("Relations matrix must have 3 columns")
	    id1 <- relations[,1]
	    id2 <- relations[,2]
	    code <- relations[,3]
	    }
	else if (is.list(relations)) {
	    id1 <- relations$id1
	    id2 <- relations$id2
	    code <- relations$code
	    if (is.null(id1) || is.null(id2) || is.null(code)) 
		stop("Relations list must have id1, id2, and code")
	    if (length(id1) != length(id2) || (length(id1) != length(code)))
		stop("Id1, id2 and code in the relations list are different lengths")
	    }
        else stop("Relations argument must be a matrix or a list")
	
	if (any(code != floor(code)) || min(code) <1 || max(code) >4)
	    stop("Invalid relationship code")
     
	# Is everyone in this relationship in the pedigree?
	temp1 <- match(id1, id, nomatch=0)
	temp2 <- match(id2, id, nomatch=0)
	if (any(temp1==0 | temp2==0))
		stop("Subjects in relationships that are not in the pedigree")
	if (any(temp1==temp2)) {
	    who <- temp1[temp1==temp2]
	    stop(paste("Subject", id[who], "is their own spouse or twin"))
	    }

	# Check, are the twins really twins?
	if (any(code<3)) {
	    twins <- (code<3)
	    if (any(momid[temp1[twins]] != momid[temp2[twins]]))
		stop("Twins with non-identical parents")
	    if (any(dadid[temp1[twins]] != dadid[temp2[twins]]))
		stop("Twins with non-identical parents")
	    }
        
	# Check, are the monozygote twins the same gender?
	if (any(code==1)) {
	    mztwins <- (code==1)
	    if (any(sex[mztwins] != sex[mztwins]))
		stop("MZ Twins with different genders")
	    }
	
	temp$relation <- cbind(temp1, temp2, code)
	}

    # autohint fills in the hints array -- our first guess about how to plot it
    temp <- c(temp, list(hints = autohint(temp)))
    oldClass(temp) <- "pedigree"
    temp
    }



"[.pedigree" <- function(x, ..., drop=F) {
    allargs <- list(...)
#   if (nDotArgs() != 1) stop ("Only 1 subscript allowed")
    if (length(allargs) != 1) stop ("Only 1 subscript allowed")
    old.n <- length(x$id)
    i <- (1:old.n)[..1]

    z <- list(id=x$id[i],momid=x$momid[i],dadid=x$dadid[i],
	     sex=x$sex[i],depth=x$depth[i],affected=x$affected[i],
	     hints=x$hints[i,])

    # For the relationship matrix it's a bit more work.
    #  It's first and second columns contain the row number (in the
    #  id list) of the two who are related.  I need to remap those to
    #  the new row numbers
    if (!is.null(x$relation)) {
	indx1 <- match(x$relation[,1], i, nomatch=0)
	indx2 <- match(x$relation[,2], i, nomatch=0)
	keep <- (indx1 >0 & indx2 >0)  #keep only if both id's are kept
	if (any(keep))
	     z$relation <- cbind(indx1[keep], indx2[keep], x$relation[keep,3])
	}
	
    oldClass(z) <- 'pedigree'
    z
    }
