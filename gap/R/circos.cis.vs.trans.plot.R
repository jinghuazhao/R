#' Circos plot of cis vs trans protein associations
#'
#' Visualizes cis- and trans-protein associations along chromosomes using a
#' circular genome plot (Circos plot).
#'
#' The function classifies protein-QTL (pQTL) hits into cis and trans
#' associations and displays genomic links between protein loci and target
#' genes.
#'
#' This version is robust to differences in column naming conventions
#' (e.g. dot vs underscore formats).
#'
#' @param f Character string. Path to a tab-delimited file containing
#'   pQTL hits. Must include columns: `CHR`, `BP`, `SNP`, and `prot`.
#'
#' @param panel Data frame. Annotation panel containing at least:
#'   `prot`, `uniprot`, `chr`, `start`, `end`, and `gene`.
#'
#' @param id Character identifier passed to
#'   `cis.vs.trans.classification()`.
#'
#' @param radius Numeric distance threshold defining the cis region
#'   (default: `1e6` base pairs).
#'
#' @return
#' A Circos plot is drawn to the current graphics device.
#' The function returns invisibly.
#'
#' @details
#' This function relies on **circlize** for visualization and assumes that
#' `cis.vs.trans.classification()` returns a list containing a data frame
#' `data` with cis/trans annotations.
#'
#' Column names may vary depending on package version. The implementation
#' attempts flexible matching of:
#'
#' - `p.chr` or `p_chr`
#' - `cis.start` or `cis_start`
#' - `cis.end` or `cis_end`
#' - `p.gene` or `p_gene`
#' - `p.prot` or `p_prot` or `prot`
#'
#' @export
#'
#' @examples
#' \dontrun{
#' f <- file.path(find.package("pQTLtools"),"tests","INF1.merge")
#' circos.cis.vs.trans.plot(f, inf1, id = "uniprot")
#' }
#'
circos.cis.vs.trans.plot <- function (f, panel, id, radius = 1e6)
{
    ## ----------------------------
    ## dependencies
    ## ----------------------------
    for (p in c("circlize", "dplyr")) {
        if (!requireNamespace(p, quietly = TRUE)) {
            warning(sprintf("Missing package: %s", p))
        }
    }

    requireNamespace("circlize")
    requireNamespace("dplyr")

    ## ----------------------------
    ## load data
    ## ----------------------------
    clumped <- read.table(f, as.is = TRUE, header = TRUE)

    hits2 <- merge(
        clumped[c("CHR", "BP", "SNP", "prot")],
        panel[c("prot", "uniprot")],
        by = "prot"
    )

    names(hits2) <- c("prot", "Chr", "bp", "SNP", "uniprot")

    ## ----------------------------
    ## classification
    ## ----------------------------
    cvt <- cis.vs.trans.classification(hits2, panel, id, radius)
    cd <- cvt$data

    if (!is.data.frame(cd)) {
        stop("cvt$data is not a data.frame")
    }

    ## ----------------------------
    ## helper: flexible column getter
    ## ----------------------------
    getcol <- function(data, options) {
        hit <- options[options %in% names(data)][1]
        if (is.na(hit)) return(NULL)
        data[[hit]]
    }

    chr   <- getcol(cd, c("p.chr", "p_chr"))
    start <- getcol(cd, c("cis.start", "cis_start"))
    end   <- getcol(cd, c("cis.end", "cis_end"))
    gene  <- getcol(cd, c("p.gene", "p_gene"))
    prot <- getcol(cd, c("p.prot", "p_prot", "prot"))

    if (any(sapply(list(chr, start, end, gene, prot), is.null))) {
        stop("Missing required columns in cvt$data. Check structure.")
    }

    ## ----------------------------
    ## circos setup
    ## ----------------------------
    circlize::circos.clear()
    circlize::circos.par(
        start.degree = 90,
        track.height = 0.1,
        cell.padding = c(0, 0, 0, 0)
    )

    circlize::circos.initializeWithIdeogram(
        species = "hg19",
        track.height = 0.05,
        ideogram.height = 0.06
    )

    ## ----------------------------
    ## annotation
    ## ----------------------------
    ann <- panel[c("chr", "start", "end", "gene")]

    ann$chr <- paste0("chr", ann$chr)
    ann$start <- pmax(0, ann$start - radius)
    ann$end <- ann$end + radius

    circlize::circos.genomicLabels(
        ann,
        labels.column = 4,
        cex = 0.8,
        font = 3,
        side = "inside"
    )

    ## ----------------------------
    ## cis anchors
    ## ----------------------------
    b1 <- data.frame(
        chr   = paste0("chr", hits2$Chr),
        start = hits2$bp - 1,
        end   = hits2$bp
    )

    ## ----------------------------
    ## cis/trans links
    ## ----------------------------
    b2 <- data.frame(
        chr   = paste0("chr", chr),
        start = start,
        end   = end,
        gene  = gene,
        prot  = prot
    )

    ## ----------------------------
    ## colors
    ## ----------------------------
    if (!"cis.trans" %in% names(cd)) {
        stop("cis.trans column missing in cvt$data")
    }

    colors <- rep(NA, nrow(cd))
    colors[cd$cis.trans == "cis"] <- 10
    colors[cd$cis.trans == "trans"] <- 12

    ## ----------------------------
    ## plot
    ## ----------------------------
    circlize::circos.genomicLink(
        b1,
        b2,
        col = colors,
        border = colors,
        directional = 1,
        lwd = 1.6
    )

    invisible(NULL)
}
