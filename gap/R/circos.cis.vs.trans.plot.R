#' circos plot of cis/trans classification
#'
#' @param hits A text file as input data with varibles named "CHR","BP","SNP","prot".
#' @param panel Protein panel with prot(ein), uniprot (id) and "chr","start","end","gene".
#' @param id Identifier.
#' @param radius The flanking distance as cis.
#'
#' @details
#' The function implements a circos plot at the early stage of SCALLOP-INF meta-analysis.
#'
#' @export
#' @return None.
#' @examples
#' \dontrun{
#'   circos.cis.vs.trans.plot(hits="INF1.merge", panel=inf1, id="uniprot")
#' }

circos.cis.vs.trans.plot <- function(hits, panel, id, radius=1e6)
{
  bp <- NA
  for(p in c("circlize","dplyr")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!requireNamespace(p, quietly = TRUE))
        warning(paste("circos.cis.vs.trans.plot needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  requireNamespace("circlize")
  requireNamespace("dplyr")
  clumped <- read.table(hits,as.is=TRUE,header=TRUE)
  hits <- merge(clumped[c("CHR","BP","SNP","prot")],panel[c("prot","uniprot")],by="prot") %>%
          setNames(c("prot","Chr","bp","SNP","uniprot"))
  cvt <- cis.vs.trans.classification(hits,panel,id,radius)
  with(cvt,summary(data))
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, track.height = 0.1, cell.padding = c(0, 0, 0, 0))
  circlize::circos.initializeWithIdeogram(species="hg19", track.height = 0.05, ideogram.height = 0.06)
  ann <- panel[c("chr","start","end","gene")]
  ann <- within(ann, {chr=paste0("chr",chr);start=start-radius;end <- end+radius})
  ann[with(ann,start<0),"start"] <- 0
  circlize::circos.genomicLabels(ann,labels.column = 4, cex=0.8, font=3, side="inside")
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
}
