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

lambda1000 <- function(lambda, ncases, ncontrols)
  1 + (lambda - 1) * (1 / ncases + 1 / ncontrols)/( 1 / 1000 + 1 / 1000)

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

sun3d <- function(xyz="INF1.merge.cis.vs.trans",
                      cols=c("id","chr1","pos1","chr2","pos2","gene","target","log10p","x","y","col"),
                      xy.scale=c(1.3e8,1.3e8),marker.size=4,log10p.max=400,
                      prefix=c("Sentinel","CHR","POS","CHR","POS","Gene","Target","-log10(p)"),
                      postfix="\u003c/br>",
                      json.file="d3.json",pretty=TRUE)
{
  for(p in c("jsonlite", "plotly")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
        warning(paste("sun3d needs package `", p, "' to be fully functional; please install", sep=""))
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
  d <- read.csv(xyz,as.is=TRUE)
  r <- pqtl2dplot(d, plot=FALSE)
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

xy <- function(x) if (x<23) x else if (x==23) "X" else if (x==24) "Y"
ixy <- function(x) if (x=="X") 23 else if (x=="Y") 24 else x

grid2d <- function(chrlen, plot=TRUE, cex=0.6)
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
  invisible(list(n=n, CM=c(0,CM)))
}

allDuplicated <- function(x)
{
  front <- duplicated(x)
  back <- duplicated(x, fromLast = TRUE)
  all_dup <- front | back
  return(all_dup)
}

mht.control <- function(type="p", usepos=FALSE, logscale=TRUE, base=10, cutoffs=NULL, colors=NULL,
                        labels=NULL, srt=45, gap=NULL, cex=0.4, yline=3, xline=3)
               list(type=type, usepos=usepos, logscale=logscale, base=base, cutoffs=cutoffs, colors=colors,
                    labels=labels, srt=srt, gap=gap, cex=cex, yline=yline, xline=xline)

hmht.control <- function(data=NULL, colors=NULL, yoffset=0.25, cex=1.5, boxed=FALSE)
                list(data=data,colors=colors,yoffset=yoffset,cex=cex,boxed=boxed)

textbox <- function(label, name=NULL, gp=NULL, vp=NULL)
{
  gt <- grid::gTree(label=label, name=name, gp=gp, vp=vp, cl="textboxtree")
  grid::grid.draw(gt)
}

makeContent.textboxtree <- function(x)
{
  t <- grid::textGrob(x$label, name="text")
  rr <- grid::roundrectGrob(width=1.5*grid::grobWidth(t), height=1.5*grid::grobHeight(t), name="box")
  grid::setChildren(x, grid::gList(t, rr))
}
