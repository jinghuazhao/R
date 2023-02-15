#' A QTL cis/trans classifier
#'
#' @param geneSNP data.frame with columns on gene, SNP and biomarker (e.g., expression, protein).
#' @param SNPPos data.frame containing SNP, chromosome and position.
#' @param genePos data.frame containing gene, chromosome, start and end positions.
#' @param radius flanking distance.
#'
#' @details
#' The function obtains QTL (simply called SNP here) cis/trans classification based on gene positions.
#'
#' @export
#' @return
#' It returns a geneSNP-prefixed data.frame with the following columns:
#' - geneChrom gene chromosome.
#' - geneStart gene start.
#' - geneEnd gene end.
#' - SNPChrom pQTL chromosome.
#' - SNPPos pQTL position.
#' - Type cis/trans labels.
#'
#' @examples
#' \dontrun{
#'   merged <- read.delim("INF1.merge",as.is=TRUE)
#'   hits <- merge(merged[c("CHR","POS","MarkerName","prot","log10p")],
#'                 inf1[c("prot","uniprot")],by="prot")
#'   names(hits) <- c("prot","Chr","bp","SNP","log10p","uniprot")
#'
#'   options(width=200)
#'   geneSNP <- merge(hits[c("prot","SNP","log10p")],
#'                    inf1[c("prot","gene")],by="prot")[c("gene","SNP","prot","log10p")]
#'   SNPPos <- hits[c("SNP","Chr","bp")]
#'   genePos <- inf1[c("gene","chr","start","end")]
#'   cvt <- qtlClassifier(geneSNP,SNPPos,genePos,1e6)
#'   cvt
#'   cistrans <- cis.vs.trans.classification(hits,inf1,"uniprot")
#'   cis.vs.trans <- with(cistrans,data)
#'   cistrans.check <- merge(cvt[c("gene","SNP","Type")],cis.vs.trans[c("p.gene","SNP","cis.trans")],
#'                           by.x=c("gene","SNP"),by.y=c("p.gene","SNP"))
#'   with(cistrans.check,table(Type,cis.trans))
#' }
#'
#' @note This is adapted from iBMQ/eqtlClassifier as an xQTL (x=e, p, me, ...) classifier.
#' @seealso [`cis.vs.trans.classification`]

qtlClassifier <- function (geneSNP, SNPPos, genePos, radius)
{
    if (dim(geneSNP)[2] < 3) {
        stop("The geneSNP data.frame need at least 3 columns.")
    }
    if (dim(SNPPos)[2] != 3) {
        stop("The SNP position data.frame need 3 columns.")
    }
    if (dim(genePos)[2] != 4) {
        stop("The gene position data.frame need 4 columns.")
    }
    res1 <- 0
    res2 <- 0
    res3 <- 0
    res4 <- 0
    res5 <- 0
    res6 <- ""
    for (i in 1:nrow(geneSNP)) {
        qtl_gene <- as.character(geneSNP[i, 1])
        bool <- toupper(as.character(genePos[, 1])) %in% toupper(qtl_gene)
        if (any(bool)) {
            genechr <- genePos[bool, 2][1]
            genestart <- genePos[bool, 3][1]
            genestop <- genePos[bool, 4][1]
        }
        else {
            cat("gene:", qtl_gene, "\t: position is missing\n")
            genechr <- "NA"
            genestart <- "NA"
            genestop <- "NA"
        }
        qtl_snp <- as.character(geneSNP[i, 2])
        bool <- toupper(as.character(SNPPos[, 1])) %in% toupper(qtl_snp)
        if (any(bool)) {
            chrsnp <- SNPPos[bool, 2][1]
            pos <- SNPPos[bool, 3][1]
        }
        else {
            cat("SNP", qtl_snp, "\t: position is missing\n")
            chrsnp <- "NA"
            pos <- "NA"
        }
        if (genechr != "NA" & chrsnp != "NA") {
            if (genechr == chrsnp) {
                if (abs(as.numeric(genestart) - as.numeric(pos)) <= 
                  radius) {
                  type = "cis"
                }
                else {
                  type = "trans"
                }
            }
            else {
                type = "trans"
            }
        }
        else {
            type = "NA"
        }
        res1 <- c(res1, genechr)
        res2 <- c(res2, genestart)
        res3 <- c(res3, genestop)
        res4 <- c(res4, chrsnp)
        res5 <- c(res5, pos)
        res6 <- c(res6, type)
    }
    res <- data.frame(res1, as.numeric(res2), as.numeric(res3), res4, as.numeric(res5), res6)
    colnames(res) <- c("geneChrom", "geneStart", "geneEnd", "SNPChrom", "SNPPos", "Type")
    return(cbind(geneSNP, res[-1, ]))
}
