allele.recode <- function (a1, a2, miss.val = NA)
{
    n <- length(a1)
    if (is.factor(a1))
        a1 <- as.character(a1)
    if (is.factor(a2))
        a2 <- as.character(a2)
    is.ch <- is.character(a1) | is.character(a2)
    if (is.ch) {
        t <- factor(c(a1, a2), exclude = miss.val)
    }
    if (!is.ch) {
        lev <- sort(unique(c(a1, a2)))
        t <- factor(c(a1, a2), levels = lev, exclude = miss.val)
    }
    allele.label <- levels(t)
    t <- as.numeric(t)
    a1 <- t[1:n]
    a2 <- t[(n + 1):(2 * n)]
    return(list(a1 = a1, a2 = a2, allele.label = allele.label))
}

geno.recode <- function (geno, miss.val = 0)
{
    n.loci <- ncol(geno)/2
    alist <- vector("list", n.loci)
    grec <- NULL
    for (i in 1:n.loci) {
        t <- (i - 1) * 2 + 1
        tmp <- allele.recode(geno[, t], geno[, t + 1], miss.val = miss.val)
        grec <- cbind(grec, tmp$a1, tmp$a2)
        alist[[i]] <- list(allele = tmp$allele.label)
    }
    return(list(grec = grec, alist = alist))
}

a2g <- function(a1,a2)
{
  i <- ifelse(a1 < a2,a2,a1)
  j <- ifelse(a1 < a2,a1,a2)
  genocode <- ifelse (j==0, 0, i*(i-1)/2 + j)
  invisible (genocode)

}

g2a.c <- function (g)
{
    d <- 1 + 8 * (g - 1)
    u <- 1 + ((1 + (sqrt(d) - 1) - 1) / 2)
    u <- ceiling(u)
    l = g - u * (u - 1) / 2
    list (l=l,u=u)
}

g2a <- function (g)
{
    i <- 1 + floor((sqrt(8 * g + 1) - 1)/2)
    j <- g - i * (i - 1)/2
    i <- ifelse(j == 0, i - 1, i)
    j <- ifelse(j == 0, i, j)
    return(cbind(j,i))
}

is.miss <- function(data,data.int,miss.val=0)
{
   if (data.int==2) # genotype
   {
      id <- array(FALSE, length(data))
      for (i in 1:length(miss.val))
          id <- id | data==miss.val[i]
   }
   else # allele
   {
      id <- array(FALSE, length(data[,1]))  
      for (i in 1:length(miss.val))
          id <- id | apply(data==miss.val[i],1,any)
   }
   return (id)  
}

# adapted from gcp.c on 24/9/2004
# 17-6-2004 JH Zhao

revhap <- function(loci,hapid)
{
   nloci <- length(loci)
   nalleles <- vector("numeric",nloci)
   nalleles[nloci] <- 1
   m <- nloci
   for(k in nloci:1)
   {
      m <- m - 1
      nalleles[m] <- nalleles[m+1] * loci[k]
   }
  n <- length(hapid)
  hap <- matrix(0,n,nloci)
  for (i in 1:n)
  {
    l <- hapid[i] - 1
    for (j in 1:nloci)
    {
      hap[i,j] <- floor(l/nalleles[j])
      if (j==nloci) hap[i,j] <- l
      else l <- l%%nalleles[j]
      hap[i,j] <- hap[i,j] + 1
    }
  }
  invisible(hap)

}

revhap.i <- function(loci,hapid)
{
   nloci <- length(loci)
   nalleles <- vector("numeric",nloci)
   nalleles[nloci] <- 1
   m <- nloci
   for(k in nloci:1)
   {
      m <- m - 1
      nalleles[m] <- nalleles[m+1] * loci[k]
   }

  hap <- vector("numeric",nloci)
  l <- hapid - 1
  for (j in 1:nloci)
  {
    hap[j] <- floor(l/nalleles[j])
    if (j==nloci) hap[j] <- l
    else l <- l%%nalleles[j]
    hap[j] <- hap[j] + 1
  }
  invisible(hap)

}

gcode <- function(a1,a2) {

  i <- ifelse(a1 < a2,a2,a1)
  j <- ifelse(a1 < a2,a1,a2)
  genocode <- i*(i-1)/2 + j
  return(genocode)

}

ungcode <- function(g) {

  i <- 1 + floor((sqrt(8*g+1)-1)/2)
  j <- g - i*(i-1)/2
  
  # the following 2 lines were added as a patch to make this ungcode work:
  i <- ifelse(j==0,i-1,i)
  j <- ifelse(j==0,i,j)

  return(cbind(j,i))
}

grec2g <- function (h, n, t)
{
  hh <- h
  for (i in 1:n)
  {
    hh[,i] <- t$alist[[i]]$allele[h[,i]]
  }
  invisible(hh)
}

# a simple scheme to represent SNPs
# similar to a2g() and not unlike Mike Weale's Matlab function
# 10-2-2005 JHZ
m2plem <- function(a1,a2)
{
  a <- vector("numeric")
# 1=A, 2=a
  for (i in 1:length(a1))
  {
      if(a1[i]==1&a2[i]==2) {
         a[i] <- 0
      }
      else if (a1[i]==1&a2[i]==1) {
         a[i] <- 1
      }
      else if (a1[i]==2&a2[i]==2) {
         a[i] <- 2
      }
  }
  invisible(a)
}

plem2m <- function(a)
{
  a1 <- vector("numeric")
  a2 <- vector("numeric")
# 1=A, 2=a
  for (i in 1:length(a))
  {
      if(a[i]==0) {
         a1[i] <- 1
         a2[i] <- 2
      }
      else if (a[i]==1) {
         a1[i] <- 1
         a2[i] <- 1
      }
      else if (a[i]==2) {
         a1[i] <- 2
         a2[i] <- 2
      }
  }
  invisible(list(a1,a2))
}

micombine <- function (est, std.err, confidence = 0.95)
{
    qstar <- est[[1]]
    for (i in 2:length(est)) {
        qstar <- cbind(qstar, est[[i]])
    }
    qbar <- apply(qstar, 1, mean)
    u <- std.err[[1]]
    for (i in 2:length(std.err)) {
        u <- cbind(u, std.err[[i]])
    }
    if (!is.null(dimnames(qstar)[[1]]))
        dimnames(u)[[1]] <- dimnames(qstar)[[1]]
    u <- u^2
    ubar <- apply(u, 1, mean)
    bm <- apply(qstar, 1, var)
    m <- dim(qstar)[2]
    tm <- ubar + ((1 + (1/m)) * bm)
    rem <- (1 + (1/m)) * bm/ubar
    nu <- (m - 1) * (1 + (1/rem))^2
    alpha <- 1 - (1 - confidence)/2
    low <- qbar - qt(alpha, nu) * sqrt(tm)
    up <- qbar + qt(alpha, nu) * sqrt(tm)
    pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
    fminf <- (rem + 2/(nu + 3))/(rem + 1)
    result <- list(est = qbar, std.err = sqrt(tm), df = nu, signif = pval,
        lower = low, upper = up, r = rem, fminf = fminf)
    result
}

VR <- function(v1,vv1,v2,vv2,c12)
{
  nv <- v2^2*vv1 + v1^2*vv2 - 2*v1*v2*c12
  dv <- v2^4
  nv/dv
}

h2G <- function(V,VCOV,verbose=TRUE)
{
  VG <- V[1]
  Ve <- V[2]
  Vp <- VG + Ve
  VVG <- VCOV[1,1]
  VVe <- VCOV[2,2]
  cVGVe <- VCOV[2,1]
  h2G <- VG / Vp
  VVp <- VVG + VVe + 2 * cVGVe
  cVpVG <- VVG + cVGVe
  Varh2G <- VR(VG,VVG,Vp,VVp,cVpVG)
  if (verbose) {
     cat("Vp =",Vp,"SE =",sqrt(VVp),"\n")
     cat("h2G =",h2G,"SE =",sqrt(Varh2G),"\n")
  }
  invisible(list(Vp=Vp,VVp=VVp,h2G=h2G,Varh2G=Varh2G))
}

h2GE <- function(V, VCOV, verbose=TRUE)
 {
   n <- length(V)
   if (n==3)
   {
     VG <- V[1]
     VGE <- V[2]
     Ve <- V[3]
     Vp <- VG + VGE + Ve
     VVG <- VCOV[1, 1]
     VVGE <- VCOV[2, 2]
     VVe <- VCOV[3, 3]
     cVGVGE <- VCOV[2, 1]
     cVGVe <- VCOV[3, 1]
     cVGEVe <- VCOV[3, 2]
     s <- VVG + VVGE + VVe
     s12 <- 2 * cVGVGE
     s13 <- 2 * cVGVe
     s23 <- 2 * cVGEVe
     VVp <- s + s12 + s13 + s23
     cVpVG <- VVG + cVGVGE + cVGVe
     cVpVGE <- cVGVGE + VVGE + cVGEVe
     h2G <- VG/Vp
     Varh2G <- VR(VG, VVG, Vp, VVp, cVpVG)
     h2GE <- VGE/Vp
     Varh2GE <- VR(VGE, VVGE, Vp, VVp, cVpVGE)
     if (verbose) {
         cat("The old function results: \n",
             "Vp =", Vp, "SE =", sqrt(VVp), "\n",
             "h2G =", h2G, "SE =", sqrt(Varh2G), "\n",
             "h2GE =", h2GE, "SE =", sqrt(Varh2GE), "\n")
     }
   }
   VarV <- diag(VCOV)
   # phenotypic variance
   P <- sum(V)
   VarP <- sum(VCOV)
   # heritabilities for individual components
   GGE <- V[-n]
   s <- vector()
   for(j in 1:(n-1)) s[j] <- VR(GGE[j],VarV[j],P,VarP,sum(VCOV[j,]))
   G <- GGE[1]
   h2G <- G/P
   Varh2G <- s[1]
   GE <- GGE[-1]
   h2GE <- GE/P
   Varh2GE <- s[-1]
   # total interaction variances/heritability
   m <- -c(1,n)
   sGE <- sum(GE)
   VarsGE <- sum(VCOV[m,m])
   CovsGEP <- sum(VCOV[1,m])+VarsGE+sum(VCOV[n,m])
   h2sGE <- sGE/P
   Varh2sGE <- VR(sGE,VarsGE,P,VarP,CovsGEP)
   # genotypic variance/heritability including interactions
   sGGE <- sum(GGE)
   VarsGGE <- sum(VCOV[-n,-n])
   CovsGGE <- sum(VCOV[n,-n])
   h2sGGE <- sGGE/P
   Varh2sGGE <- VR(sGGE,VarsGGE,P,VarP,VarsGGE+CovsGGE)
   if (verbose) {
      cat("\nVariance, SE\n")
      source <- c("G",paste("GxE",1:(n-2),sep=""),"E")
      print(data.frame(V,SE=sqrt(diag(VCOV)),row.names=source))
      cat("\nV(G) (SE) =", sGGE, sqrt(VarsGGE),
          "including V(GE) (SE) =", sGE, sqrt(VarsGE),
          "\nV(P) (SE) =", P, sqrt(VarP), "\n")
      cat("\nHeritability, SE\n")
      h2 <- data.frame(h2=c(h2G,h2GE),SE=sqrt(s),row.names=source[-n])
      print(h2)
      cat("\nsum of h2GE (SE) =", h2sGE, sqrt(Varh2sGE),
          "\nsum of h2G and h2GE (SE) =", h2sGGE, sqrt(Varh2sGGE), "\n")
   }
   invisible(list(G=G, VarG=VarV[1],
                  GE=GE, VarGE=VarV[m],
                  P=P, VarP=VarP,
                  sGE=sGE, VarsGE=VarsGE,
                  sGGE=sGGE, VarsGGE=VarsGGE,
                  h2G=h2G, Varh2G=Varh2G,
                  h2GE=h2GE, Varh2GE=Varh2GE,
                  h2sGE=h2sGE, Varh2sGE=Varh2sGE,
                  h2sGGE=h2sGGE, Varh2sGGE=Varh2sGGE))
}

h2l <- function(K=0.05,P=0.5,h2,se,verbose=TRUE)
{
  x <- qnorm(1-K)
  z <- dnorm(x)
  1/sqrt(2*pi)*exp(-x^2/2)
  fK <- (K*(1-K)/z)^2
  fP <- P*(1-P)
  f <- fK/fP
  h2l <- f*h2
  sel <- f*se
  z2 <- K^2*(1-K)^2/(f*fP)
  x2 <- -log(2*pi*z2)
  if (verbose) {
     cat("K = ", K, "P = ", P, "\n")
     cat("h2 =",h2,"SE =",se,"h2l =",h2l,"SE =",sel,"\n")
  }
  invisible(list(h2=h2,se=se,h2l=h2l,sel=sel,cc=f,z=sqrt(x2)))
}

# R script to read the GRM binary file

ReadGRMBin <- function(prefix, AllN=FALSE, size=4)
{
  BinFileName <- paste(prefix,".grm.bin",sep="")
  NFileName <- paste(prefix,".grm.N.bin",sep="")
  IDFileName <- paste(prefix,".grm.id",sep="")
  id <- read.table(IDFileName)
  n <- dim(id)[1]
  BinFile <- file(BinFileName, "rb");
  grm <- readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  close(BinFile)
  NFile <- file(NFileName, "rb");
  if(AllN) N <- readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  else N <- readBin(NFile, n=1, what=numeric(0), size=size)
  close(NFile)
  i <- sapply(1:n, function(i) i*(i+1)/2)
  GRM <- matrix(NA,n,n)
  GRM[upper.tri(GRM,diag=TRUE)] <- grm
  GRM[lower.tri(GRM)] <- t(GRM)[lower.tri(GRM)]
  invisible(list(grm=grm, id=id, N=N, GRM=GRM))
}

WriteGRMBin <- function(prefix, grm, N, id, size=4)
{
  BinFileName <- paste(prefix,".grm.bin",sep="")
  NFileName <- paste(prefix,".grm.N.bin",sep="")
  IDFileName <- paste(prefix,".grm.id",sep="")
  grm.bin <- file(BinFileName, "wb")
  writeBin(grm,grm.bin,size=size)
  close(grm.bin)
  grm.N.bin <- file(NFileName, "wb")
  writeBin(N,grm.N.bin,size=size)
  close(grm.N.bin)
  write.table(id,IDFileName,col.names=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
}

ReadGRM <- function(prefix=51)
{
  idfile <- paste(prefix,".grm.id",sep="")
  id <- read.delim(idfile,header=FALSE,sep="\t")
  N <- dim(id)[1]
  L <- N * (N + 1) /2
  grmfile <- paste(prefix,".grm.gz",sep="")
  gz <- gzfile(grmfile)
  grm_lines <- readLines(gz)
  close(gz)
  GRM_list <- sapply(grm_lines,strsplit,"\t")
  M <- as.integer(lapply(GRM_list,"[[",3))
  grm <- as.numeric(lapply(GRM_list,"[[",4))
  GRM <- matrix(NA,N,N)
  GRM[upper.tri(GRM,diag=TRUE)] <- grm
  GRM[lower.tri(GRM)] <- t(GRM)[lower.tri(GRM)]
  invisible(list(GRM=GRM,N=M,id=id,grm=grm))
}

WriteGRM <- function(prefix=51,id,N,GRM)
{
  M <- N
  N <- dim(GRM)[1]
  k2l <- matrix(NA,N*(N+1)/2,4)
  L <- 1
  for(i in 1:N)
  {
    for(j in 1:i)
    {
      k2l[L,] <- c(i,j,M[L],GRM[i,j])
      L <- L + 1
    }
  }
  idfile <- paste(prefix,".grm.id",sep="")
  write(t(id),file=idfile,2,sep="\t")
  grmfile <- paste(prefix,".grm.gz",sep="")
  gz <- gzfile(grmfile,"w")
  write.table(k2l,file=gz,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  close(gz)
}

ReadGRMPLINK <- function(prefix,diag=1)
{
  f <- paste(prefix,".genome",sep="")
  g <- read.table(f,header=TRUE,as.is=TRUE)
  L <- dim(g)[1]
  N <- (1+sqrt(1+8*L))/2
  pi_hat <- g["PI_HAT"]
  PIHAT <- matrix(NA,N,N)
  diag(PIHAT) <- diag
  PIHAT[lower.tri(PIHAT)] <- pi_hat[,1]
  PIHAT[upper.tri(PIHAT)] <- t(PIHAT)[upper.tri(PIHAT)]
  idplink <- function (N)
  {
    id <- function(N, diag=TRUE)
    {
      irep <- function(i) rep(i,i)
      iter <- function(i) 1:i
      i <- unlist(lapply(1:N,irep))
      j <- unlist(lapply(lapply(1:N,iter),cbind))
      ij <- cbind(i,j)
      if (diag) ij
      else subset(ij,ij[,1]!=ij[,2])
    }
    ij <- id(N, diag=FALSE)
    ij[order(ij[,2]),2:1]
  }
  fid1 <- g[1,"FID1"]
  iid1 <- g[1,"IID1"]
  fid <- c(fid1,g[1:(N-1),"FID2"])
  iid <- c(iid1,g[1:(N-1),"IID2"])
  ji <- idplink(N)
  idn <- 1:N
  ii <- cbind(i=idn,j=idn)
  ij <- rbind(ji,ii)
  idx <- ij[order(ij[,2]),2:1]
  invisible(list(pihat=PIHAT[lower.tri(PIHAT,diag=TRUE)],PIHAT=PIHAT,id=cbind(fid,iid),idx=idx))
}

WriteGRMSAS <- function(grmlist,outfile="gwas")
{
  for(p in c("foreign")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
        warning(paste("WriteGRMSAS needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  with(grmlist,
  {
    parm <- 1
    colnames(GRM) <- row <- paste("COL",1:dim(id)[1],sep="")
    names(id) <- c("idn","id")
    write.dta(data.frame(id,parm,row,GRM),paste(outfile,".dta",sep=""))
  })
}

ReadGRMPCA <- function(prefix)
{
  values <- scan(paste(prefix,".eigenval",sep=""))
  vectors <- read.table(paste(prefix,".eigenvec",sep=""))
  N <- dim(vectors)[1]
  id <- vectors[,1:2]
  vectors <- as.matrix(vectors[,-(1:2)])
  s <- (values>0)
  posGRM <- vectors[,s]%*%diag(values[s])%*%t(vectors[,s])
  s <- (values<0)
  negGRM <- vectors[,s]%*%diag(values[s])%*%t(vectors[,s])
  invisible(list(N=N,posGRM=posGRM,negGRM=negGRM,id=id))
}

# GLGC-GIANT code for inverse normal transformation on x with missing data 

invnormal <- function(x)
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))) 

# log(p) for a normal deviate z

logp <- function(z) log(2) + pnorm(-abs(z), lower.tail = TRUE, log.p = TRUE)

# log10(p) for a normal deviate z

log10p <- function(z)
  log(2, base=10)+pnorm(-abs(z), lower.tail=TRUE, log.p=TRUE)/log(10)

# gc.lambda and miamiplot functions hosted at CEU by Daniel R Barnes
# A simplified version is as follows,
# obs <- median(chisq)
# exp <- qchisq(0.5, 1) # 0.4549364
# lambda <- obs/exp
# see also estlambda from GenABEL and qq.chisq from snpStats

gc.lambda <- function(p) {
  p <- p[!is.na(p)]
  n <- length(p)

  obs <- qchisq(p,1,lower.tail=FALSE)
  exp <- qchisq(1:n/n,1,lower.tail=FALSE)

  lambda <- median(obs)/median(exp)
  return(lambda)
}

lambda1000 <- function(lambda, ncases, ncontrols)
  1 + (lambda - 1) * (1 / ncases + 1 / ncontrols)/( 1 / 1000 + 1 / 1000)

miamiplot <- function (x, chr = "CHR", bp = "BP", p = "P", pr = "PR", snp = "SNP", 
    col = c("midnightblue", "chartreuse4"), col2 = c("royalblue1", 
        "seagreen1"), ymax = NULL, highlight = NULL, highlight.add = NULL, 
    pch = 19, cex = 0.75, cex.lab = 1, xlab = "Chromosome", ylab = "-log10(P) [y>0]; log10(P) [y<0]", 
    lcols = c("red", "black"), lwds = c(5, 2), ltys = c(1, 2), 
    main = "", ...) 
{
    P = index = NULL
    PR = index = NULL
    if (!(chr %in% names(x))) 
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x))) 
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x))) 
        stop(paste("Column", p, "not found!"))
    if (!(pr %in% names(x))) 
        stop(paste("Column", pr, "not found!"))
    if (!(snp %in% names(x))) 
        if (!is.numeric(x[[chr]])) 
            stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]])) 
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]])) 
        stop(paste(p, "column should be numeric."))
    if (!is.numeric(x[[pr]])) 
        stop(paste(pr, "column should be numeric."))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
        PR = x[[pr]])
    if (!is.null(x[[snp]])) 
        d = transform(d, SNP = x[[snp]])
    d = subset(d[order(d$CHR, d$BP), ], (P > 0 & P <= 1 & is.numeric(P) & 
        PR > 0 & PR <= 1 & is.numeric(PR)))
    d$logp = -log10(d$P)
    d$logpr = log10(d$PR)
    d$pos = NA
    ymax = ceiling(max(d$logp) + 2)
    ymin = floor(min(d$logpr) - 2)
    d$index = NA
    ind = 0
    for (i in unique(d$CHR)) {
        ind = ind + 1
        d[d$CHR == i, ]$index = ind
    }
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
        d$pos = d$BP
        ticks = floor(length(d$pos))/2 + 1
        xlabel = paste("Chromosome", unique(d$CHR), "position")
        labs = ticks
    }
    else {
        lastbase = 0
        ticks = NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i, ]$pos = d[d$index == i, ]$BP
            }
            else {
                lastbase = lastbase + tail(subset(d, index == 
                  i - 1)$BP, 1)
                d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
                  lastbase
            }
            ticks = c(ticks, d[d$index == i, ]$pos[floor(length(d[d$index == 
                i, ]$pos)/2) + 1])
        }
        xlabel = "Chromosome"
        labs = unique(d$CHR)
    }
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)
    plot(NULL, xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
        xlim = c(xmin, xmax), ylim = c(ymin, ymax), main = main, 
        xlab = xlab, ylab = ylab, las = 1, pch = pch, cex.lab = cex.lab, 
        ...)
    if (nchr == 1) {
        axis(1, ...)
    }
    else {
        axis(1, at = ticks, labels = labs, ...)
    }
    col = rep(col, max(d$CHR))
    col2 = rep(col2, max(d$CHR))
    if (nchr == 1) {
        with(d, points(pos, logpr, cex = cex, pch = pch, ...))
        with(d, points(pos, logp, cex = cex, pch = pch, ...))
    }
    else {
        icol = 1
        for (i in unique(d$index)) {
            with(d[d$index == unique(d$index)[i], ], points(pos, 
                logpr, col = col2[icol], cex = cex, pch = pch, 
                ...))
            with(d[d$index == unique(d$index)[i], ], points(pos, 
                logp, col = col[icol], cex = cex, pch = pch, 
                ...))
            icol = icol + 1
        }
    }
    abline(h = -log10(5e-08), col = lcols[2], lwd = lwds[2], 
        lty = ltys[2])
    abline(h = log10(5e-08), col = lcols[2], lwd = lwds[2], lty = ltys[2])
    abline(h = 0, col = lcols[1], lwd = lwds[1], lty = ltys[1])
    if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP))) 
            warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight = d[which(d$SNP %in% highlight), ]
        with(d.highlight, points(pos, logp, col = "red3", cex = cex, 
            pch = pch, ...))
        with(d.highlight, points(pos, logpr, col = "red3", cex = cex, 
            pch = pch, ...))
    }
    if (!is.null(highlight.add)) {
        print("yessssssssssssssssss")
        if (any(!(highlight.add %in% d$SNP))) 
            warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight.add = d[which(d$SNP %in% highlight.add), 
            ]
        with(d.highlight.add, points(pos, logp, col = "darkgreen", 
            cex = cex, pch = pch, ...))
        with(d.highlight.add, points(pos, logpr, col = "darkgreen", 
            cex = cex, pch = pch, ...))
    }
}

cis.vs.trans.classification <- function(hits, panel, id, radius=1e6)
# cis.vs.trans.classification(hits=jma.cojo, panel=inf1, id="uniprot")
{
  p.start <- p.end <- Chr <- p.chr <- bp <- dist.inds <- same.inds <- NA

# "Thu Nov  8 12:13:07 2018"

# author jp549@cam.ac.uk

# identify cis vs trans hits

# rule: a cis acting variant lies within the region
# from 1MB upstream of the start position to 1MB downstream of the end position 
# of the gene that encodes the protein being tested in the GWAS

# All signals that are outside this window will be defined as trans

# add a prefix 'p.' so we know these cols refer to the protein being GWAS'd

  colnames(panel) <- paste0("p.", colnames(panel))

# map on to the hits file, using UniProtID as the common reference

  hits_panel <- merge(x=hits, y=panel, by.x=id, by.y=paste0('p.',id), all.x=TRUE)

# classify into cis and trans

# set cis as -1MB upstream to +1MB downstream

  N <- nrow(hits_panel)
  hits_panel <- within(hits_panel,
  { 
    cis.start <- p.start - radius
    if (any(cis.start < 0 )) cis.start[which(cis.start<0)] <- 0
    cis.end <- p.end + radius

# any variant on a different chromosome to the gene encoding the target protein is not cis

    dist.inds <<- which(Chr != p.chr)
    cis <- rep(NA, N)
    if (length(dist.inds)>0)  cis[dist.inds] <- FALSE

# for ones on the same chr, we can't be sure without looking at position

    same.inds <<- which(Chr == p.chr)

# see if variant lies in the cis region

    if (length(same.inds)>0) cis[same.inds] <- bp[same.inds] > cis.start[same.inds] & bp[same.inds] < cis.end[same.inds]
    cis.trans <- rep(NA, N)
    cis.trans[cis==TRUE] <- "cis"
    cis.trans[cis==FALSE] <- "trans"
  })

# split by protein

  list.by.prot <- split(hits_panel, f=with(hits_panel,p.gene))

# get the breakdown of cis vs trans per protein
# sapply(list.by.prot, function(x) table(with(x, cis.trans)))

  x <- with(hits_panel,table(p.gene, cis.trans))
  s <- sum(x)
  total <- apply(x,2,sum)
  xx <- rbind(x,total)
  total <- apply(xx,1,sum)
  x <- cbind(xx,total)
  invisible(list(data=hits_panel,table=x,total=s))
}

cnvplot <- function(data)
# cnvplot(cnv)
{
  d <- within(data,{chr<-replace(chr,chr=="X",23); chr<-replace(chr,chr=="Y",24)})
  pos <- vector("numeric")
  n <- length(table(with(data,chr)))
  for (x in 1:n) pos[x] <- with(subset(d,chr==paste(x)),{max(end)})
  CM <- cumsum(pos)
  par(xaxt = "n", yaxt = "n")
  xy <- xy.coords(c(0,CM), seq(1,90,by=90/(1+n)))
  plot(xy$x, xy$y, type = "n", ann = FALSE, axes = FALSE)
  colors <- rep(c("red","blue"),n)
  par(xaxt = "s", yaxt = "s", xpd = TRUE)
  xy <- function(x) if (x<23) x else if (x==23) "X" else if (x==24) "Y";
  for (x in 1:n) with(subset(d,chr==paste(x)), {
      l <- ifelse(x==1,0,CM[x-1])
      segments(l+start,freq,l+end,freq,lwd="3",col=colors[x])
      text(ifelse(x == 1, CM[x]/2, (CM[x] + CM[x-1])/2), 0, pos = 1, offset = 0.5, xy(x), cex=0.4)
  })
  segments(0,0,CM[x],0)
  axis(2,line=-0.5)
  title(xlab="Chromosome",ylab="Frequency",line=2)
}

circos.cnvplot <- function(data)
# circos.cnvplot(cnv)
{
  for(p in c("circlize")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
        warning(paste("circos.cnvplot needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  start <- end <- NA
  cnv <- within(subset(data,start<=end),{chr=paste0("chr",chr)})
  requireNamespace("circlize")
  circlize::circos.par(start.degree = 50, track.height = 0.3, cell.padding = c(0, 0, 0, 0))
  circlize::circos.initializeWithIdeogram(species="hg19", track.height = 0.05, ideogram.height = 0.06)
  circlize::circos.genomicTrackPlotRegion(cnv, ylim = c(0, 50), panel.fun = function(region,value,...) {
                      color <- as.numeric(gsub("chr","",circlize::get.current.chromosome()))
                      with(cbind(region,value),circlize::circos.segments(start,freq,end,freq,col=color,lwd=1))
  })
  circlize::circos.clear()
}

circos.cis.vs.trans.plot <- function(hits, panel, id, radius=1e6)
# circos.cis.vs.trans.plot(hits="INF1.clumped", panel=inf1, id="uniprot")
{
  bp <- NA
  for(p in c("circlize")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
        warning(paste("circos.cis.vs.trans.plot needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  requireNamespace("circlize")
  clumped <- read.table(hits,as.is=TRUE,header=TRUE)
  hits <- merge(clumped[c("CHR","BP","SNP","prot")],panel[c("prot","uniprot")],by="prot")
  names(hits) <- c("prot","Chr","bp","SNP","uniprot")
  cvt <- cis.vs.trans.classification(hits,panel,id,radius)
  with(cvt,summary(data))
  circlize::circos.par(start.degree = 90, track.height = 0.1, cell.padding = c(0, 0, 0, 0))
  circlize::circos.initializeWithIdeogram(species="hg19", track.height = 0.05, ideogram.height = 0.06)
  ann <- panel[c("chr","start","end","gene")]
  ann <- within(ann, {chr=paste0("chr",chr);start=start-radius;end <- end+radius})
  ann[with(ann,start<0),"start"] <- 0
  circlize::circos.genomicLabels(ann,labels.column = 4, side="inside")
  b1 <- with(cvt,data)[c("Chr","bp")]
  b1 <- within(b1,{Chr=paste0("chr",Chr);start=bp-1})
  names(b1) <- c("chr","end","start")
  b2 <- with(cvt,data)[c("p.chr","cis.start","cis.end","p.gene","p.prot")]
  b2 <- within(b2,{p.chr=paste0("chr",p.chr)})
  names(b2) <- c("chr","start","end","gene","prot")
  colors <- rep(NA,nrow(with(cvt,data)))
  colors[with(cvt,data)["cis.trans"]=="cis"] <- 10
  colors[with(cvt,data)["cis.trans"]=="trans"] <- 12
  circlize::circos.genomicLink(b1, b2, col = colors, border = colors, directional=1, lwd = 1.6)
  circlize::circos.clear()
}

circos.mhtplot <- function(data, glist)
# g <- c("IRS1","SPRY2","FTO","GRIK3","SNED1","HTR1A","MARCH3","WISP3","PPP1R3B","RP1L1","FDFT1","SLC39A14","GFRA1","MC4R")
# circos.mhtplot(mhtdata,g)
{
  pos <- gene <- NULL
  for(p in c("circlize")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
        warning(paste("circos.mhtplot needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  requireNamespace("circlize")
  d <- within(data, {
    chr <- paste0("chr",chr)
    start <- pos - 1
    end <- pos
  })[c("chr","start","end","p","gene")]
  hd <-	subset(data, gene %in% glist)[c("chr","start","end","gene")]
  hd <-	within(hd, {chr	<- paste0("chr",chr)})
  ann <- data.frame()
  for (g in glist)
  {
     m <- subset(hd, gene==g)
     n <- round(nrow(m) / 2 + 0.5)
     ann <- rbind(ann,m[n,])
  }
  circlize::circos.par(start.degree = 90, track.height = 0.4, cell.padding = c(0, 0, 0, 0))
  circlize::circos.initializeWithIdeogram(species = "hg18", track.height = 0.05, ideogram.height = 0.06)
  circlize::circos.genomicTrackPlotRegion(d[c("chr","start","end","p")], ylim = c(0, 15), 
         panel.fun = function(region, value, ...) {
           color <- as.numeric(gsub("chr", "", circlize::get.current.chromosome()))
           with(cbind(region, value), circlize::circos.genomicPoints(region, -log10(value), cex=0.3, col = color))
  })
  circlize::circos.genomicLabels(ann, labels.column = 4, side = "inside")
  circlize::circos.clear()
}

cs <- function(tbl, b="Effect", se="StdErr", log_p=NULL, cutoff=0.95)
# credible set based on METAL sumstats
{
  requireNamespace("matrixStats")
  tbl <- within(tbl, {
           if (is.null(log_p)) z <- tbl[[b]]/tbl[[se]]
           else z <- qnorm(tbl[[log_p]], log.p=TRUE)
           z2 <- z * z / 2
           d <- matrixStats::logSumExp(z2)
           log_ppa <- z2 - d
           ppa <- exp(log_ppa)
        })
  ord <- with(tbl, order(ppa,decreasing = TRUE))
  tbl_ord <- within(tbl[ord,], {cppa <- cumsum(ppa)})
  last <- which(with(tbl_ord,cppa) >= cutoff)[1]
  tbl[ord[1:last],]
}

xy <- function(x) if (x<23) x else if (x==23) "X" else if (x==24) "Y"

ixy <- function(x) if (x=="X") 23 else if (x=="Y") 24 else x

hg38 <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,
          114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415)
hg19 <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,
          115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566)
hg18 <- c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,146274826,140273252,135374737,134452384,132349534,
          114142980,106368585,100338915,88827254,78774742,76117153,63811651,62435964,46944323,49691432,154913754,57772954)

grid2d <- function(chrlen=hg19, plot=TRUE, cex=0.6)
{
  CM <- cumsum(chrlen)
  n <- length(chrlen)
  xy <- xy.coords(c(0,CM), c(0,CM))
  if (plot)
  {
    par(xaxt = "n", yaxt = "n")
    plot(xy$x, xy$y, type = "n", ann = FALSE, axes = FALSE)
    par(xaxt = "s", yaxt = "s", xpd = TRUE)
    for (x in 1:n) {
        segments(CM[x],0,CM[x],CM[n],col="black")
        segments(0,CM[x],CM[n],CM[x],col="black")
        text(ifelse(x == 1, CM[x]/2, (CM[x] + CM[x-1])/2), 0, pos = 1, offset = 0.5, xy(x), cex=cex)
        text(0, ifelse(x == 1, CM[x]/2, (CM[x] + CM[x-1])/2), pos = 2, offset = 0.5, xy(x), cex=cex)
    }
    segments(0,0,CM[n],0)
    segments(0,0,0,CM[n])
    title(xlab="pQTL position",ylab="protein position",line=2)
  }
  return(list(n=n, CM=c(0,CM)))
}

mhtplot2d <- function(data, plot=TRUE, cex=0.6)
# mhtplot2d(read.table("INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE))
{
  r <- grid2d(plot=plot)
  n <- with(r, n)
  CM <- with(r, CM)
  d <- data[c("SNP","Chr","bp","p.chr","p.start","p.end","p.target.short","p.gene","log10p","cis")]
  chr1 <- d[["Chr"]]
  chr1[chr1=="X"] <- 23
  chr1[chr1=="Y"] <- 24
  pos1 <- CM[chr1] + d[["bp"]]
  chr2 <- d[["p.chr"]]
  chr2[chr2=="X"] <- 23
  chr2[chr2=="Y"] <- 24
  mid <- (d[["p.start"]] + d[["p.end"]])/2
  pos2 <- CM[chr2] + mid
  if (plot) {
     points(pos1,pos2,cex=cex,col=ifelse(d[["cis"]],"red","blue"),pch=19)
     legend("top",legend=c("cis","trans"),box.lty=0,cex=cex,col=c("red","blue"),horiz=TRUE,inset=c(0,1),xpd=TRUE,pch=19)
  }
  return(data.frame(id=d[["SNP"]],chr1=chr1,pos1=d[["bp"]],chr2=chr2,pos2=mid,x=pos1,y=pos2,
         target=d[["p.target.short"]],gene=d[["p.gene"]],log10p=with(d,log10p),col=ifelse(d[["cis"]],"blue","red")))
}

snptest_sample <- function(data,sample_file="snptest.sample",ID_1="ID_1",ID_2="ID_2",missing="missing",C=NULL,D=NULL,P=NULL)
{
  cat(ID_1,ID_2,missing,C,D,P,file=sample_file)
  cat("\n",file=sample_file,append=TRUE)
  len_C <- length(C)
  len_D <- length(D)
  len_P <- length(P)
  cat("0 0 0",rep("C",len_C),rep("D",len_D),rep("P",len_P),file=sample_file,append=TRUE)
  cat("\n",file=sample_file,append=TRUE)
  write.table(data[c(ID_1,ID_2,missing,C,D,P)],file=sample_file,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
}

chr_pos_a1_a2 <- function(chr,pos,a1,a2,prefix="chr",seps=c(":","_","_"),uppercase=TRUE)
{
  chr <- paste0(prefix,chr)
  chrpos <- paste(chr,pos,sep=seps[1])
  a1a2 <- paste(a1,a2,sep=seps[3])
  a2a1 <- paste(a2,a1,sep=seps[3])
  swap <- (a1 > a2)
  a1a2[swap] <- a2a1[swap]
  a1a2.lower <- tolower(a1a2)
  a1a2.upper <- toupper(a1a2)
  if(uppercase) paste(chrpos,a1a2.upper,sep=seps[2]) else paste(chrpos,a1a2.lower,sep=seps[2])
}

inv_chr_pos_a1_a2 <- function(chr_pos_a1_a2,prefix="chr",seps=c(":","_","_"))
{
  if ((seps[1]==seps[2])&(seps[2]==seps[3]))
  {
    s <- sapply(chr_pos_a1_a2,strsplit,seps[1])
    chr <- lapply(s,"[",1)
    pos <- lapply(s,"[",2)
    a1 <- lapply(s,"[",3)
    a2 <- lapply(s,"[",4)
  } else if ((seps[1]!=seps[2])&(seps[2]==seps[3]))
  {
    s <- sapply(chr_pos_a1_a2,strsplit,seps[2])
    chrpos <- lapply(s,"[",1)
    s1 <- sapply(chrpos,strsplit,seps[1])
    chr <- lapply(s1,"[",1)
    pos <- lapply(s1,"[",2)
    a1 <- lapply(s,"[",2)
    a2 <- lapply(s,"[",3)
  } else if ((seps[1]!=seps[2])&(seps[2]!=seps[3]))
  {
    s <- sapply(chr_pos_a1_a2,strsplit,seps[2])
    chrpos <- lapply(s,"[",1)
    s1 <- sapply(chrpos,strsplit,seps[1])
    chr <- lapply(s1,"[",1)
    pos <- lapply(s1,"[",2)
    s2 <- lapply(s,"[",2)
    s3 <- sapply(s2,strsplit,seps[3])
    a1 <- lapply(s3,"[",1)
    a2 <- lapply(s3,"[",2)
  }
  if (prefix=="") chr <- gsub("chr","",chr)
  s <- data.frame(chr=unlist(chr),pos=unlist(pos),a1=unlist(a1),a2=unlist(a2))
  names(s) <- c("chr","pos","a1","a2")
  return(s)
}

mhtplot3d <- function(xyz="INF1.merge.cis.vs.trans",
                      cols=c("id","chr1","pos1","chr2","pos2","gene","target","log10p","x","y","col"),
                      xy.scale=c(1.3e8,1.3e8),marker.size=4,log10p.max=400,
                      prefix=c("Sentinel","CHR","POS","CHR","POS","Gene","Target","-log10(p)"),
                      postfix="\u003c/br>",
                      json.file="d3.json",pretty=TRUE)
{
  for(p in c("jsonlite", "plotly")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
        warning(paste("mhtplot3d needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  src <- list(
    x = list(
      layout = list(
        margin = list(b = 40, l = 60, t = 25, r = 10),
        scene = list(
          xaxis = list(title = "pQTL position",
                       tickmode = "array", 
                       autotick = FALSE, 
                       tick0 = 1, 
                       dtick = 1, 
                       ticklen = 0, 
                       tickwidth = 0, 
                       tickfont = list(size = 10),
                       tickvals = as.list(1:24),
                       ticktext = as.list(c(1:22,"X","Y"))
          ),
          yaxis = list(title = "Protein position",
                       tickmode = "array",
                       autotick = FALSE,
                       tick0 = 1,
                       dtick = 1,
                       ticklen = 0,
                       tickwidth = 0,
                       tickfont = list (size = 10),
                       tickvals = as.list(1:24),
                       ticktext = as.list(c(1:22,"X","Y"))
          ),
          zaxis = list(title = "-log10(p)", tickfont = list(size = 10)),
          camera = list(eye=list(x = -1.3, y = -1.2, z = 1.1)),
          aspectmode = "manual",
          aspectratio = list(x = 0.9, y = 1, z = 0.6)
        ),
        legend = list(x = 10, y = 0.5),
        xaxis = list(domain=list(0,1)),
        yaxis = list(domain=list(0,1)),
        title = "Scatterplot of sentinels",
        showlegend = TRUE
      ),
      source = "A",
      config = list(
        modeBarButtonsToAdd = list(name = "Collaborate"),
        modeBarButtonsToRemove = list("sendDataToCloud")
      ),
      data <- list(),
      base_url = "https://plot.ly"
    ),
    evals = list(),
    jsHooks = list()
  )
  d <- read.table(xyz,as.is=TRUE,header=TRUE)
  r <- mhtplot2d(d, plot=FALSE)
  cuts <- with(r, abs(log10p) > log10p.max)
  r <- within(r,{x=x/xy.scale[1]; y=y/xy.scale[2]; log10p[!cuts] <- abs(log10p[!cuts]); log10p[cuts] <- log10p.max})
  fixes <- function(col,d) paste(paste(prefix[col],d[,col],sep=":"),postfix)
  cistrans <- function(br,name,color)
  {
    t <- subset(r[cols],col==br)
    cuts <- with(t, abs(log10p) == log10p.max)
    n.cuts <- with(t,length(log10p[cuts]))
    s <- rep('circle',nrow(t))
    s[cuts] <- 'diamond'
    list(x=t$x, y=t$y, z=t$log10p,
         text=as.list(apply(sapply(c(1:8),fixes,t),1,paste,collapse=" ")),
         type="scatter3d",
         mode="markers", 
         name=name,
         marker=list(color=color,symbol=s,size=marker.size))
  }
  cis <- cistrans("blue","cis","rgba(191,56,42,1)")
  trans <- cistrans("red","trans","rgba(12,75,142,1)")
  src$x$data <- list(cis,trans)
  if(!is.null(json.file))
  {
    json <- jsonlite::toJSON(src,auto_unbox=TRUE,pretty=pretty)
    sink(json.file)
    print(json)
    sink()
  }
  p <- plotly::plot_ly()
  p <- with(src$x$data[[1]],plotly::add_trace(p, x=x, y=y, z=z, marker=marker, mode=mode, name=name, text=text, type=type))
  p <- with(src$x$data[[2]],plotly::add_trace(p, x=x, y=y, z=z, marker=marker, mode=mode, name=name, text=text, type=type))
  p <- with(src$x$layout, plotly::layout(p, scene=scene, xaxis=xaxis, yaxis=yaxis, margin=margin, title=title, showlegend=showlegend))
}

# https://plot.ly/r/reference/#scatter3d
# sed -i 's|<\\/br>|\\u003c/br>|g' d3.json
# plotly::toRGB( c('#BF382A', '#0C4B8E')) ==> "rgba(191,56,42,1)" "rgba(12,75,142,1)"

pvalue <- function (z, decimals = 2)
{
    lp <- -log10p(z)
    exponent <- ceiling(lp)
    base <- 10^-(lp - exponent)
    paste0(round(base, decimals), "e", -exponent)
}

log10pvalue <- function(p=NULL,base=NULL,exponent=NULL)
{
  if(!is.null(p))
  {
    p <- format(p,scientific=TRUE)
    p2 <- strsplit(p,"e")
    base <- as.numeric(lapply(p2,"[",1))
    exponent <- as.numeric(lapply(p2,"[",2))
  } else if(is.null(base) | is.null(exponent)) stop("base and exponent should both be specified")
  log10(base)+exponent
}

# Adapted from code by Felix Day 16/9/2015
# with reference from Bowden J, et al. (2015). 44(2):512-525.
weighted.median <- function(x, w)
{
   x.order <- x[order(x)]
   w.order <- w[order(x)]
   ws <- cumsum(w.order)-0.5*w.order
   ws <- ws / sum(w.order)
   below <- max(which(ws < 0.5))
   x.order[below] + (x.order[below+1]-x.order[below]) * (0.5-ws[below]) / (ws[below+1]-ws[below])
}

mr.boot = function(bXG, sebXG, bYG, sebYG, w, n.boot=1000, method="median")
{
   boot <- vector('numeric')
   for(i in 1:n.boot)
   {
       bXG.boot <- rnorm(length(bXG), mean = bXG, sd = sebXG)
       bYG.boot <- rnorm(length(bYG), mean = bYG, sd = sebYG)
       if(method=="Egger")
       {
         bYG.boot = bYG.boot*sign(bXG.boot)
         bXG.boot = abs(bXG.boot)
         boot[i] = coef(summary(lm(bYG.boot~bXG.boot,weights=sebYG^-2)))[2,1]
       } else {
         bIV.boot <- bYG.boot / bXG.boot
         boot[i] <- weighted.median(bIV.boot, w)
      }
   }
   sd(boot)
}

gsmr <- function(data, X, Y, alpha=0.05, other_plots=FALSE)
{
   c <- qnorm(alpha/2,lower.tail=FALSE)
   nxy <- vector("list", length(X)*length(Y))
   k <- 1
   ES_plot <- vector('list')
   for(i in Y)
   {
       dat <- subset(data.frame(SNP = data[["SNP"]],
                                bzx = eval(parse(text = paste0("data$b.", X))), SEzx = eval(parse(text = paste0("data$SE.", X))),
                                bzy = eval(parse(text = paste0("data$b.", i))), SEzy = eval(parse(text = paste0("data$SE.", i)))),
                    !is.na(bzx+SEzx+bzy+SEzy))
       SNP <- dat[["SNP"]]
       bzx <- dat[["bzx"]]; SEzx <- dat[["SEzx"]]
       bzy <- dat[["bzy"]]; SEzy <- dat[["SEzy"]]
       nxy[k] <- paste0(X, ".", i)
       m1 <- lm(bzy~bzx-1, weights=SEzy^-2); summary.m1 <- summary(m1); coef.summary.m1 <- coef(summary.m1)
       m <- lm(bzy~bzx, weights=SEzy^-2); summary.m <- summary(m); coef.summary.m <- coef(summary.m)
       bIVW <- coef.summary.m1[1,1]
       sebIVW <- coef.summary.m1[1,2] / min(with(summary.m1, sigma), 1)
       bEGGER <- coef.summary.m[2,1]
       sebEGGER <- coef.summary.m[2,2] / min(with(summary.m, sigma), 1)
       intEGGER <- coef.summary.m[1,1]
       seintEGGER <- coef.summary.m[1,2]
       bIV <- bzy / bzx
       weights <- (SEzy / bzx)^-2
       bIVW <- sum(bzy*bzx*SEzy^-2) / sum(bzx^2 * SEzy^-2)
       bWM <- weighted.median(bIV, weights)
       sebWM <- mr.boot(bzx, SEzx, bzy, SEzy, weights)
       penalty <- pchisq(weights*(bIV-bIVW)^2, df=1, lower.tail=FALSE)
       penalty.weights <- weights*pmin(1, penalty * 20)
       bPWM <- weighted.median(bIV, penalty.weights)
       sebPWM <- mr.boot(bzx, SEzx, bzy, SEzy, penalty.weights)
       res <- metafor::rma(bIV, 1/sqrt(weights), method="FE")
       CochQ <- res$QE
       CochQp <- res$QEp
       graph_title <- paste("SNP effect on", i)
       x_title <- paste("Effect on", X)
       y_title <- paste("Effect on", i)
       if(other_plots)
       {
         metafor::funnel(res, main=paste(graph_title), xlab=y_title)
         abline(v=0, lty="dashed", col="red")
         metafor::forest(res, main=paste(graph_title), addfit=FALSE, slab=SNP, xlab=y_title)
       }
       bzxLCL <- bzx-c*SEzx; bzxUCL <- bzx+c*SEzx
       bzyLCL <- bzy-c*SEzy; bzyUCL <- bzy+c*SEzy
       ES_plot[[i]] <- ggplot2::ggplot(data.frame(bzx, SEzx, bzy, SEzy, bzxLCL, bzxUCL, bzyLCL, bzyUCL),
                                ggplot2::aes(x=bzx, y=bzy)) +
                                ggplot2::geom_point() +
                                ggplot2::theme_bw() +
                                ggplot2::geom_errorbar(ggplot2::aes(ymin=bzyLCL, ymax=bzyUCL)) +
                                ggplot2::geom_errorbarh(ggplot2::aes(xmin=bzxLCL, xmax=bzxUCL)) +
                                ggplot2::geom_abline(intercept=0, slope=0, size=1) +
                                ggplot2::geom_abline(intercept=0, slope=bIVW, size=1, colour="red", show.legend=TRUE) +
                                ggplot2::geom_abline(intercept=0, slope=bEGGER, size=1, colour="blue", show.legend=TRUE) +
                                ggplot2::geom_abline(intercept=0, slope=bWM, size=1, colour="green", show.legend=TRUE) +
                                ggplot2::geom_abline(intercept=0, slope=bPWM, size=1, colour="orange", show.legend=TRUE) +
                                ggplot2::geom_vline(xintercept=0, size=1) +
                                cowplot::theme_half_open(12) +
                                ggplot2::theme(plot.margin = ggplot2::margin(6, 4, 6, 4)) +
                                ggplot2::ggtitle(graph_title) +
                                ggplot2::xlab(x_title) +
                                ggplot2::ylab(y_title)
       plot(ES_plot[[i]])
       ColOut <- parse(text = paste0(X, ".", i))
       assign(paste(ColOut),c(paste(ColOut),bIVW,sebIVW,CochQ,CochQp,bEGGER,sebEGGER,intEGGER,seintEGGER,bWM,sebWM,bPWM,sebPWM))
       k <- k+1
  }
  r <- c("IV","bIVW","sebIVW","CochQ","CochQp","bEGGER","sebEGGER","intEGGER","seintEGGER","bWM","sebWM","bPWM","sebPWM")
  for (p in nxy) r <- rbind(r, eval(parse(text = paste(p))))
  invisible(list(r=r,plots=ES_plot))
}
