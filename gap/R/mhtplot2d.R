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
  invisible(list(n=n, CM=c(0,CM)))
}

#' 2D Manhattan plot
#'
#' @md
#' @param d Data to be used.
#' @param snp_name variant name.
#' @param snp_pos variant position.
#' @param snp_chr variant chromosome.
#' @param gene_chr gene chromosome.
#' @param gene_start gene start position.
#' @param gene_end gene end position.
#' @param protein protein name.
#' @param gene gene name.
#' @param lp log10(p).
#' @param cis cis variant when TRUE.
#' @param plot to plot when TRUE.
#' @param cex extension factor.
#' @export
#' @return positional information.
#' @examples
#'
#' INF <- Sys.getenv("INF")
#' d <- read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans"),as.is=TRUE)
#' r <- mhtplot2d(d)

mhtplot2d <- function(d, snp_name="SNP", snp_chr="Chr", snp_pos="bp",
                      gene_chr="p.chr", gene_start="p.start", gene_end="p.end",
                      protein="p.target.short", gene="p.gene", lp="log10p",
                      cis="cis",
                      plot=TRUE, cex=0.6)
{
  r <- grid2d(plot=plot)
  n <- with(r, n)
  CM <- with(r, CM)
  chr1 <- d[[snp_chr]]
  chr1[chr1=="X"] <- 23
  chr1[chr1=="Y"] <- 24
  pos1 <- CM[chr1] + d[[snp_pos]]
  chr2 <- d[[gene_chr]]
  chr2[chr2=="X"] <- 23
  chr2[chr2=="Y"] <- 24
  mid <- (d[[gene_start]] + d[[gene_end]])/2
  pos2 <- CM[chr2] + mid
  if (plot) {
     points(pos1, pos2, cex=cex, col=ifelse(d[[cis]],"red","blue"), pch=19)
     legend("top", legend=c("cis","trans"), box.lty=0, cex=cex, col=c("red","blue"),
            horiz=TRUE, inset=c(0,1), xpd=TRUE, pch=19)
  }
  return(list(n=n, CM=CM, data=data.frame(id=d[[snp_name]],
                                          chr1=chr1, pos1=d[[snp_pos]],
                                          chr2=chr2, pos2=mid, x=pos1, y=pos2,
                                          target=d[[protein]], gene=d[[gene]], log10p=d[[lp]],
                                          col=ifelse(d[[cis]],"blue","red")
  )))
}
