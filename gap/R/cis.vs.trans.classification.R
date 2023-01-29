#' A cis/trans classifier
#'
#' @param hits Data to be used, which contains prot, Chr, bp, id and/or other information such as SNPid.
#' @param panel Panel data.
#' @param id Identifier.
#' @param radius The flanking distance for variants.
#'
#' @details
#' The function classifies variants into cis/trans category according to a panel which contains id, chr, start, end, gene variables.
#'
#' @export
#' @return
#' The cis/trans classification.
#' @examples
#' cis.vs.trans.classification(hits=jma.cojo, panel=inf1, id="uniprot")
#' \dontrun{
#' INF <- Sys.getenv("INF")
#' f <- file.path(INF,"work","INF1.merge")
#' clumped <- read.delim(f,as.is=TRUE)
#' hits <- merge(clumped[c("CHR","POS","MarkerName","prot","log10p")],
#'               inf1[c("prot","uniprot")],by="prot")
#' names(hits) <- c("prot","Chr","bp","SNP","log10p","uniprot")
#' cistrans <- cis.vs.trans.classification(hits,inf1,"uniprot")
#' cis.vs.trans <- with(cistrans,data)
#' knitr::kable(with(cistrans,table),caption="Table 1. cis/trans classification")
#' with(cistrans,total)
#' }
#' @author James Peters

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
    cis.trans[cis] <- "cis"
    cis.trans[!cis] <- "trans"
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
