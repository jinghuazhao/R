#  $Id: plot.pedigree.s,v 1.10 2003/10/24 14:13:34 Atkinson Exp $

plot.pedigree <- function(x, id = x$id, sex = x$sex, status = x$status, 
			  affected = x$affected, 
			  cex = 1, col = rep(1, length(x$id)), 
			  symbolsize = 1, branch = 0.6, 
			  packed = T, align = packed, width = 8, 
			  density=c(5, 0, 100, 50), mar=c(4.1, 1, 4.1, 1),
			  angle=c(90, 70, 50, 0), keep.par=F, ...) # density = c(-1, 50, 70, 90)
{
    maxlev <- max(x$depth) + 1
    n <- length(x$depth)	


    ## DATA CHECKS 
    sex <- as.numeric(sex)
    sex[is.na(sex)] <- 3
    if(any(sex < 1 | sex > 4))
      stop("Invalid sex code")
    if(length(sex) != n)
      stop("Wrong length for sex")

    if(is.null(status))
      status <- rep(0, n)
    else {
	if(!all(status == 0 | status == 1 | status == 2))
	  stop("Invalid status code")
	if(length(status) != n)
	  stop("Wrong length for status")
    }
    if(!is.null(id)) {
	if(length(id) != n)
	  stop("Wrong length for id")
    }

    ## if someone has an unknown affected status, assign them 
    ## a value of 0 (not affected)

    # More data checks -- do these if affected/status/relations are filled in
    #   (the optional parts of a pedigree structure)

    if(is.null(affected)){
      affected <- matrix(0,nrow=n)
    }
    else {
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
	if(max(affected) > min(affected)) 
              affected <- matrix(affected - min(affected),nrow=n)
        else
              affected <- matrix(affected,nrow=n)
	if (!all(affected==0 | affected==1 | affected==2))
		stop("Invalid code for affected status")
	}
  
    if(!missing(status)) {
	if(length(status) != n)
		stop("Wrong length for affected")
	if(any(status != 0 & status != 1))
		stop("Invalid status code")
	}

    
    ## Define plotting symbols for sex by affected
    ## in S, 0=square, 1= circle, 5=diamond, 2=triangle
    ## 15, 16, 18 and 17 are filled versions of them

    symbol <-  c(0, 1, 5, 2)[sex]

    ## CREATE SYMBOLS/PLOTS

    points.sym <- function(x, y, symbol, status, affected, size, col, aspect,
			   angle, density, adj)
      {

	  circle <- function(cx, cy, r, code=0, col=1, angle, density, ...)
	    {

		#   cx, cy, coordinates for center; r is radius

		z <- (0:360 * pi)/180
		pin <- par()$pin
		usr <- par()$usr
		adj <- (pin[2]/pin[1])/((usr[4] - usr[3])/(usr[2] - usr[1]))
		x <- sin(z) * r + cx
		y <- cos(z) * r * 1/adj + cy

		if(sum(code,na.rm=T)==0) polygon(x, y, border=T, 
			density=density[1],col=col, ...) #density=0
		else {

		    if(length(code)==1) polygon(x,y,border=T,col=col,
			       density=density[code[1]+1], angle=angle[1],...) #density[1]

		    if(length(code)==2) {
			polygon(x,y,border=T,density=0, ...)
			z <- (0:180 * pi)/180
			x <- sin(z) * r + cx
			y <- cos(z) * r * 1/adj + cy

			polygon(x,y,border=T,col=col,density=density[code[2]+1],
				angle=angle[2], ...) # density[2]*code[2]
			z <- (180:360 * pi)/180
			x <- sin(z) * r + cx
			y <- cos(z) * r * 1/adj + cy

			polygon(x,y,border=T,col=col,density=density[code[1]+1],
				angle=angle[1], ...) # code[1]*density[1]
		    }

		    if(length(code)==3) {
			polygon(x,y,border=T,density=0, ...)
			
			z <- (0:90 * pi)/180
			x <- c(cx,sin(z) * r + cx)
			y <- c(cy,cos(z) * r * 1/adj + cy)
			polygon(x,y,border=T,col=col,density=density[code[3]+1],
				angle=angle[3],...) # code[3]*density[3]

			z <- (180:270 * pi)/180
			x <- c(cx,sin(z) * r + cx)
			y <- c(cy,cos(z) * r * 1/adj + cy)
			polygon(x,y,border=T,col=col,density=density[code[1]+1],
				angle=angle[1], ...) # code[1]*density[1]

			z <- (270:360 * pi)/180
			x <- c(cx, sin(z) * r + cx)
			y <- c(cy, cos(z) * r * 1/adj + cy)
			polygon(x,y,border=T,col=col,density=density[code[2]+1],
				angle=angle[2], ...) # code[2]*density[2]
		    }

		    if(length(code)==4) {
			polygon(x,y,border=T,density=0, ...)
			
			z <- (0:90 * pi)/180
			x <- c(cx,sin(z) * r + cx)
			y <- c(cy,cos(z) * r * 1/adj + cy)
			polygon(x,y,border=T,col=col,density=density[code[3]+1],
				angle=angle[3], ...) # code[3]*density[3]

			z <- (90:180 * pi)/180
			x <- c(cx,sin(z) * r + cx)
			y <- c(cy,cos(z) * r * 1/adj + cy)
			polygon(x,y,border=T,col=col,density=density[code[4]+1], 
				angle=angle[4],...) # code[4]*density[4]

			z <- (180:270 * pi)/180
			x <- c(cx,sin(z) * r + cx)
			y <- c(cy,cos(z) * r * 1/adj + cy)
			polygon(x,y,border=T,col=col,density=density[code[1]+1],
				angle=angle[1], ...) # code[1]*density[1]

			z <- (270:360 * pi)/180
			x <- c(cx, sin(z) * r + cx)
			y <- c(cy, cos(z) * r * 1/adj + cy)
			polygon(x,y,border=T,col=col,density=density[code[2]+1],
				angle=angle[2], ...) # code[2]*density[2]
		    }

		    if(length(code)>4) stop('Can only plot up to 4 levels of codes')
		}  
		invisible()
	    }

	  square <- function(cx, cy, r, code=0, col=1, angle, density, ...) {

	      #  cx, cy, coordinates for center; r is radius

	      pin <- par()$pin
	      usr <- par()$usr
	      adj <- (pin[2]/pin[1])/((usr[4] - usr[3])/(usr[2] - usr[1]))
	      
	      x <- cx + c(-r,-r,r,r)
	      y <- cy + (1/adj)*c(-r,r,r,-r)

	      if(sum(code,na.rm=T)==0) polygon(x,y,border=T,density=density[1], col=col,...) # density=0
	      else {
		  if(length(code)==1) polygon(x,y,border=T,col=col,
			     density=density[code[1]+1], angle=angle[1], ...) # density[1]

		  if(length(code)==2) {
		      polygon(x,y,border=T,density=0, ...)
		      
		      x <- cx + c(-r,-r,0,0)
		      polygon(x,y,border=T,col=col,density=density[code[1]+1],
			      angle=angle[1], ...) # code[1]*density[1]
		      x <- cx + c(0,0,r,r)
		      polygon(x,y,border=T,col=col,density=density[code[2]+1],
			      angle=angle[2], ...) # density[2]*code[2]
		  }
		  if(length(code)==3) {
		      polygon(x,y,border=T,density=0, ...)
		      
		      x <- cx + c(-r,-r,0,0)
		      y <- cy + (1/adj)*c(-r,0,0,-r)
		      polygon(x,y,border=T,col=col,density=density[code[1]+1],
			      angle=angle[1], ...) # code[1]*density[1]
		      x <- cx + c(-r,-r,0,0)
		      y <- cy + (1/adj)*c(0,r,r,0)
		      polygon(x,y,border=T,col=col,density=density[code[2]+1],
			      angle=angle[2], ...) # code[2]*density[2]
		      x <- cx + c(0,0,r,r)
		      y <- cy + (1/adj)*c(0,r,r,0)
		      polygon(x,y,border=T,col=col,density=density[code[3]+1],
			      angle[3], ...) # code[3]*density[3]
		  }
		  if(length(code)==4) {
		      polygon(x,y,border=T,density=0, ...)
		      
		      x <- cx + c(-r,-r,0,0)
		      y <- cy + (1/adj)*c(-r,0,0,-r)
		      polygon(x,y,border=T,col=col,density=density[code[1]+1],
			      angle=angle[1], ...) # code[1]*density[1]
		      x <- cx + c(-r,-r,0,0)
		      y <- cy + (1/adj)*c(0,r,r,0)
		      polygon(x,y,border=T,col=col,density=density[code[2]+1],
			      angle=angle[2], ...) # code[2]*density[2]
		      x <- cx + c(0,0,r,r)
		      y <- cy + (1/adj)*c(0,r,r,0)
		      polygon(x,y,border=T,col=col,density=density[code[3]+1],
			      angle=angle[3], ...) # code[3]*density[3]
		      x <- cx + c(0,0,r,r)
		      y <- cy + (1/adj)*c(-r,0,0,-r)
		      polygon(x,y,border=T,col=col,density=density[code[4]+1],
			      angle[4], ...) # code[4]*density[4]
		  }
		  if(length(code)>4) stop('Can only plot up to 4 levels of codes')
	      }  
	      invisible()
	  }

	  diamond <- function(cx, cy, r, code=0, col=1, angle, density, ...) {

	      #  cx, cy, coordinates for center; r is radius

	      pin <- par()$pin
	      usr <- par()$usr
	      adj <- (pin[2]/pin[1])/((usr[4] - usr[3])/(usr[2] - usr[1]))
	      
	      x <- cx + c(-r,0,r,0)
	      y <- cy + (1/adj)*c(0,r,0,-r)

	      if(sum(code,na.rm=T)==0) polygon(x,y,border=T,density=density[1],col=col, ...) # density=0
	      else {
		  if(length(code)==1) polygon(x,y,border=T,col=col,
			     density=density[code[1]+1], angle=angle[1], ...) # density[1]

		  if(length(code)==2) {
		      polygon(x,y,border=T,density=0, ...)
		      
		      x <- cx + c(-r,0,0)
		      y <- cy + (1/adj)*c(0,r,-r)
		      polygon(x,y,border=T,col=col,density=density[code[1]+1],
			      angle=angle[1], ...) # code[1]*density[1]
		      x <- cx + c(0,0,r)
		      y <- cy + (1/adj)*c(r,-r,0)
		      polygon(x,y,border=T,col=col,density=density[code[2]+1],
			      angle=angle[2], ...) # code[2]*density[2]
		  }
		  if(length(code)==3) {
		      polygon(x,y,border=T,density=0, ...)
		      
		      x <- cx + c(-r,0,0)
		      y <- cy + (1/adj)*c(0,0,-r)
		      polygon(x,y,border=T,col=col,density=density[code[1]+1],
			      angle=angle[1], ...) # code[1]*density[1]
		      x <- cx + c(-r,0,0)
		      y <- cy + (1/adj)*c(0,r,0)
		      polygon(x,y,border=T,col=col,density=density[code[2]+1],
			      angle=angle[2], ...) # code[2]*density[2]
		      x <- cx + c(0,0,r)
		      y <- cy + (1/adj)*c(0,r,0)
		      polygon(x,y,border=T,col=col,density=density[code[3]+1],
			      angle=angle[3], ...) # code[3]*density[3]
		  }
		  if(length(code)==4) {
		      polygon(x,y,border=T,density=0, ...)
		      
		      x <- cx + c(-r,0,0)
		      y <- cy + (1/adj)*c(0,0,-r)
		      polygon(x,y,border=T,col=col,density=density[code[1]+1],
			      angle=angle[1], ...) # code[1]*density[1]
		      x <- cx + c(-r,0,0)
		      y <- cy + (1/adj)*c(0,r,0)
		      polygon(x,y,border=T,col=col,density=density[code[2]+1],
			      angle=angle[2], ...) # code[2]*density[2]
		      x <- cx + c(0,0,r)
		      y <- cy + (1/adj)*c(0,r,0)
		      polygon(x,y,border=T,col=col,density=density[code[3]+1],
			      angle=angle[3], ...) # code[3]*density[3]
		      x <- cx + c(0,0,r)
		      y <- cy + (1/adj)*c(-r,0,0)
		      polygon(x,y,border=T,col=col,density=density[code[4]+1],
			      angle=angle[4], ...) # code[4]*density[4]
		  }
		  if(length(code)>4) stop('Can only plot up to 4 levels of codes')
	      }  
	      invisible()
	  }

	  triangle <- function(cx, cy, r, code=0, col=1, angle, density, ...) {

	      #  cx, cy, coordinates for center; r is radius

	      pin <- par()$pin
	      usr <- par()$usr
	      adj <- (pin[2]/pin[1])/((usr[4] - usr[3])/(usr[2] - usr[1]))
	      
	      r <- r*1.25

	      a <- 3*r/sqrt(3)
	      b <- r/2

	      x <- cx + c((-1/2)*a, (1/2)*a, 0)
	      y <- cy + (1/adj)*c(-b, -b, r)

	      if(sum(code,na.rm=T)==0) polygon(x,y,border=T,density=density[1],col=col, ...) # density=0
	      else {
		  if(length(code)==1) polygon(x,y,border=T,col=col,
			     density=density[code[1]+1],angle=angle[1], ...) # density[1]

		  if(length(code)==2) {
		      polygon(x,y,border=T,density=0, ...)
		      
		      x <- cx + c((-1/2)*a, 0, 0)
		      y <- cy + (1/adj)*c(-b, -b, r)
		      polygon(x,y,border=T,col=col,density=density[code[1]+1],
			      angle=angle[1], ...) # code[1]*density[1]
		      
		      x <- cx + c(0, (1/2)*a, 0)
		      y <- cy + (1/adj)*c(-b, -b, r)
		      polygon(x,y,border=T,col=col,density=density[code[2]+1],
			      angle=angle[2], ...) # code[2]*density[2]
		  }

		  if(length(code)==3) {
		      polygon(x,y,border=T,density=0, ...)
		      midx <- (r*(.5)*a)/(b+r)

		      x <- cx + c((-1/2)*a,-midx, 0, 0)
		      y <- cy + (1/adj)*c(-b, 0, 0, -b)
		      polygon(x,y,border=T,col=col,density=density[code[1]+1],
			      angle=angle[1], ...) # code[1]*density[1]
		      x <- cx + c(-midx, 0,0)
		      y <- cy + (1/adj)*c(0,r,0)
		      polygon(x,y,border=T,col=col,density=density[code[2]+1],
			      angle=angle[2], ...) # code[2]*density[2]
		      x <- cx + c(0,0,midx)
		      y <- cy + (1/adj)*c(0,r,0)
		      polygon(x,y,border=T,col=col,density=density[code[3]+1],
			      angle=angle[3],...) # code[3]*density[3]
		  }
		  if(length(code)==4) {
		      polygon(x,y,border=T,density=0, ...)
		      midx <- (r*(.5)*a)/(b+r)
		      
		      x <- cx + c((-1/2)*a, -midx , 0, 0)
		      y <- cy + (1/adj)*c(-b, 0, 0, -b)
		      polygon(x,y,border=T,col=col,density=density[code[1]+1],
			      angle=angle[1], ...) # code[1]*density[1]
		      x <- cx + c(-midx, 0,0)
		      y <- cy + (1/adj)*c(0,r,0)
		      polygon(x,y,border=T,col=col,density=density[code[2]+2],
			      angle=angle[2], ...) # code[2]*density[2]
		      x <- cx + c(0,0,midx)
		      y <- cy + (1/adj)*c(0,r,0)
		      polygon(x,y,border=T,col=col,density=density[code[3]+1],
			      angle=angle[3], ...) # code[3]*density[3]
		      x <- cx + c(0,midx,(1/2)*a,0)
		      y <- cy + (1/adj)*c(0,0,-b,-b)
		      polygon(x,y,border=T,col=col,density=density[code[4]+1],
			      angle=angle[4], ...) # code[4]*density[4]
		  }
		  if(length(code)>4) stop('Can only plot up to 4 levels of codes')
	      }  
	      invisible()
	  }

	  #   in S, 0=square, 1= circle, 5=diamond, 2=triangle
	  #   15, 16, 18 and 17 are filled versions of them

	  for(i in (1:length(x))[!is.na(match(symbol,c(0,15)))]){
	      square(x[i],y[i],size, code=(affected[i,,drop=F]), 
		     col=col[i], angle=angle, density=density)
	  }

	  for(i in (1:length(x))[!is.na(match(symbol,c(1,16)))]){
	      circle(x[i],y[i],size, code=(affected[i,,drop=F]), 
		     col=col[i],angle=angle, density=density)
	  }

	  for(i in (1:length(x))[!is.na(match(symbol,c(5,18)))]){
	      diamond(x[i],y[i],size, code=(affected[i,,drop=F]), 
		      col=col[i],angle=angle, density=density)
	  }

	  for(i in (1:length(x))[!is.na(match(symbol,c(2,17)))]){
	      triangle(x[i],y[i],size, code=(affected[i,,drop=F]), 
		       col=col[i],angle=angle, density=density)
	  }


          # Draw slash if status=1
          who <- (status == 1)
	  if(any(who)) {
	      deltax <- size*cos(pi/4) + aspect
	      deltay <- (1/adj)*(size*sin(pi/4) + aspect)
	      segments(x[who] - deltax, y[who] - deltay,
		       x[who] + deltax, y[who] + deltay)
	  }
      }
      
      #
      # Now, get the structure of the plot, and lay out the main region
      #
	
    plist <- align.pedigree(x, packed = packed, width = width, align = 
			    align)
    who <- (plist$nid > 0)
    xp <- plist$pos[who]
    yp <-  - (row(plist$nid))[who]
    np <- plist$nid[who]
    oldpar <- par(mar = mar, err=-1)
    par(usr = c(range(xp) + c(-0.5, 0.5), - (maxlev + 0.5), -0.5))

    pin <- par()$pin
    usr <- par()$usr
    adj <- (pin[2]/pin[1])/((usr[4] - usr[3])/(usr[2] - usr[1]))

    symbolsize <- symbolsize * adj
    radius <- 0.08 * symbolsize	

    delta <- .1
    aspect <- 0.05 * symbolsize  #used for status variable	

#   textoff <- - (1/adj)*(radius) - .5*cex*(par()$"1em"[2])

    textoff <- - (1/adj)*(radius) - .5*cex*(par()$"cin"[2]/7)

    plot(xp,yp,axes=F,type='n',xlab='',ylab='',...)

    points.sym(xp, yp, symbol[np], status[np], affected[np, ,drop=F], 
	       radius, col[np], aspect, angle=angle, density=density, adj=adj)
    if(!is.null(id)) {
	text(xp, yp + textoff, id[np], cex = cex, col = col[np])
	
    }

    # Draw in the family connections
    for(i in 1:maxlev) {
	# between spouses
	if(any(plist$spouse[i,  ]>0)) {
	    temp <- (1:ncol(plist$spouse))[plist$spouse[i,  ]>0]
	    segments(plist$pos[i, temp] + radius, rep( - i, length(temp)), 
		     plist$pos[i, temp + 1] - radius, rep(-i, length(temp)))
            # to be added -- double segment for plist$spouse==2
	}
    }
    for(i in 2:maxlev) {

	# Children to parents, and across sibs
	zed <- unique(plist$fam[i,  ])
	zed <- zed[zed > 0]

	for(fam in zed) {
	    xx <- plist$pos[i - 1, fam + 0:1]
	    parentx <- mean(xx)
	    # Draw uplines from each subject
            if(!is.null(plist$twins)){
            tw.left <- (plist$twins[i, ] > 0 & plist$twins[i, ] < 4) & plist$fam[i,]==fam
            mz.left <- (plist$twins[i, ] == 1) & plist$fam[i,]==fam

            un.left <- (plist$twins[i, ] == 3) & plist$fam[i,]==fam

            ## figure out twins IN FAMILY (find next person in family)
            tw.right <- rep(F,length(tw.left))
            famlst <- plist$fam[i,] 
            twloc <- (1:length(tw.right))[tw.left]
            famloc <- (1:length(tw.right))[famlst==fam]
            flag <- NULL
            for(j in twloc) flag <- c(flag,min(famloc[j<famloc]))
            tw.right[flag] <- T

            mz.right <- rep(F,length(mz.left))
            mzloc <- (1:length(mz.right))[mz.left]
            flag <- NULL
            for(j in mzloc) flag <- c(flag,min(famloc[j<famloc]))
            mz.right[flag] <- T

            ## unknown zygosity of the twin (type=3)
            un.right <- rep(F,length(un.left))
            unloc <- (1:length(un.right))[un.left]
            flag <- NULL
            for(j in unloc) flag <- c(flag,min(famloc[j<famloc]))
            un.right[flag] <- T

          }
            
            if(is.null(plist$twins)){
              twn <- length(plist$fam[i,]==fam)
              tw.left <- rep(F,twn)
              tw.right <- rep(F,twn)
              mz.left <- rep(F,twn)
              mz.right <- rep(F,twn)

              un.left <- rep(F,twn)
              un.right <- rep(F,twn)
}

            tw <- tw.left|tw.right
            mz <- mz.left|mz.right          
            un <- un.left|un.right          
	    who <- (plist$fam[i,  ] == fam) & !tw

	    xx <- plist$pos[i, who]
	    yy <- rep( - i, length = sum(who))
	    ww <- plist$nid[i, who]

            xx.l <- plist$pos[i, tw.left]
            yy.l <- rep( - i, length = sum(tw.left))
	    ww.l <- plist$nid[i, tw.left]

            xx.r <- plist$pos[i, tw.right]
            yy.r <- rep( - i, length = sum(tw.right))
	    ww.r <- plist$nid[i, tw.right]
            
            segments(xx, yy + (1/adj) * radius, xx, yy + 3 * delta, col
		     = col[ww])

            #####################################################
            # draw diag line for twins
             ## left side of twin
	    segments(xx.l, yy.l + (1/adj) * radius, .5*(xx.l+xx.r), yy.l + 3 * delta, 
		     col = col[ww.l])

             ## right side of twin
	    segments(xx.r, yy.r + (1/adj) * radius, .5*(xx.l+xx.r), yy.r + 3 * delta, 
		     col = col[ww.r])

             ## draw midpoint MZ twin line

            xx.lm <- plist$pos[i, mz.left]
            yy.lm <- rep( - i, length = sum(mz.left))

            xx.rm <- plist$pos[i, mz.right]
            yy.rm <- rep( - i, length = sum(mz.right))

            segments(.5*(xx.lm+.5*(xx.rm+xx.lm)),
                     .5*(yy.lm+(1/adj)*radius + yy.lm+3*delta),
                     .5*(xx.rm+.5*(xx.rm+xx.lm)),
                     .5*(yy.rm+(1/adj)*radius + yy.lm+3*delta))
            
             ## add question mark for unknown zygosity

            xx.lu <- plist$pos[i, un.left]
            yy.lu <- rep( - i, length = sum(un.left))

            xx.ru <- plist$pos[i, un.right]
            yy.ru <- rep( - i, length = sum(un.right))

            ## text('?',.5*(xx.lu+.5*(xx.ru+xx.lu)),
            ##      .5*(yy.lu+(1/adj)*radius + yy.lu+3*delta), cex=adj*.3)
                    
            #####################################################

            
	    who <- (plist$fam[i,  ] == fam) 
	    xx <- plist$pos[i, who]
	    yy <- rep( - i, length = sum(who))
	    ww <- plist$nid[i, who]

            # add connector
	    xx2 <- plist$pos[i,]
            xx2 <- xx2[who]

            ## check to see if twins on either end of family
            if(sum(tw)>=2){
              tw.who <- tw[who]
              n <- length(tw.who)
              n2 <- sum(tw.who)
              flagL <- sum(tw.who[1:2]) == 2
              flagR <- sum(tw.who[c(n,n-1)]) == 2
              xx2.save <- xx2

              if(flagL) xx2[tw.who][1:2] <- rep(mean(xx2.save[tw.who][1:2]),2)
              if(flagR) xx2[tw.who][c(n2,n2-1)] <- rep(mean(xx2.save[tw.who][c(n2,n2-1)]),2)
              }                                   
           
	    segments(min(xx2), 3 * delta - i, max(xx2), 3 * delta - i)	
            
  	    # Draw line to parents
	    x1 <- mean(range(xx))
	    y1 <- 3 * delta - i

	    if(branch == 0)
	      segments(x1, y1, parentx,  - (i - 1))
	    else {
		y2 <- 1 - i
		x2 <- parentx
		ydelta <- ((y2 - y1) * branch)/2
		segments(c(x1, x1, x2), c(y1, y1 + ydelta, y2 - ydelta), 
			 c(x1, x2, x2), c(y1 + ydelta, y2 - ydelta, y2))
	    }



          }
    }

    # Connect up remote spouses
    arcconnect <- function(x, y)
      {
        # draw an approx arc whose midpoint is 1/2 unit higher than its chord, 
	#  with x[1], y[1], and x[2], y[2] as the ends of the chord
	  xx <- seq(x[1], x[2], length = 15)
	  yy <- seq(y[1], y[2], length = 15) + 0.5 - (seq(-7, 7))^2/98
	  lines(xx, yy, lty = 2)
      }

    xx <- table(np)
    xx <- xx[xx > 1]	

    #those who appear multiple times
    if(length(xx) > 0) {
	multi <- as.numeric(names(xx))
	for(i in 1:length(xx)) {
	    x2 <- xp[np == multi[i]]
	    y2 <- yp[np == multi[i]]
	    nn <- xx[i]
	    for(j in 1:(nn - 1))
	      arcconnect(x2[j + 0:1], y2[j + 0:1])
	}
    }

    ### print message if not plotting everyone

    ckall <- x$id[is.na(match(x$id,x$id[plist$nid]))]
    if(length(ckall>0)) cat('Did not plot the following people:',ckall,'\n')
    
    if(!keep.par) par(oldpar)

    ## reorder xp and yp so that they are in the same order as plist
    tmp <- plist$nid[plist$nid!=0]
    xp2 <- xp[order(tmp)]
    yp2 <- yp[order(tmp)]

    invisible(list(plist=plist,object=x,x=xp2,y=yp2,textoff=textoff,
                   symbolsize=symbolsize))
}
