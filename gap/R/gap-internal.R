allDuplicated <- function(x)
{
  front <- duplicated(x)
  back <- duplicated(x, fromLast = TRUE)
  all_dup <- front | back
  return(all_dup)
}

#' C version of g2c
#' @noRd

g2a.c <- function (g)
{
    d <- 1 + 8 * (g - 1)
    u <- 1 + ((1 + (sqrt(d) - 1) - 1) / 2)
    u <- ceiling(u)
    l = g - u * (u - 1) / 2
    list (l=l,u=u)
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

#' Allele indices for a given haplotype ID in a multiallelic system
#' @section Usage:
#' revhap(loci,hapid)
#' @noRd
#'
#' @keywords internal

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

#' as a2g
#' @noRd

gcode <- function(a1,a2) {

  i <- ifelse(a1 < a2,a2,a1)
  j <- ifelse(a1 < a2,a1,a2)
  genocode <- i*(i-1)/2 + j
  return(genocode)

}

#' Recovery of alleles from genotype(s)
#' @noRd

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

#' an experimental function for PLEM format
#' @noRd

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

#' PLEM format
#' @noRd

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

#' A function to combine imputation results
#' @noRd

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

#' Variance of a ratio
#' @noRd
#'
#' @keywords internal

VR <- function(v1,vv1,v2,vv2,c12)
{
  nv <- v2^2*vv1 + v1^2*vv2 - 2*v1*v2*c12
  dv <- v2^4
  nv/dv
}

lambda1000 <- function(lambda, ncases, ncontrols)
  1 + (lambda - 1) * (1 / ncases + 1 / ncontrols)/( 1 / 1000 + 1 / 1000)

#' 3D Manhattan plot according to Sun, et al. (2018).
#'
#' @section Usage:
#' sun3d(xyz="INF1.merge.cis.vs.trans",
#'       cols=c("id","chr1","pos1","chr2","pos2","gene","target","log10p","x","y","col"),
#'       xy.scale=c(1.3e8,1.3e8),marker.size=4,log10p.max=400,
#'       prefix=c("Sentinel","CHR","POS","CHR","POS","Gene","Target","-log10(p)"),
#'       postfix="\u003c/br>",
#'       json.file="d3.json",pretty=TRUE)
#' @noRd
#'
#' @keywords internal

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
  r <- qtl2dplot(d, plot=FALSE)
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
