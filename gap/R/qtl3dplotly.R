#' 3D QTL plot (cis/trans genomic visualization)
#'
#' Interactive 3D visualization of QTL–gene relationships using a scatter3D plot.
#'
#' This function extends `gap::qtl2dplot()` into a 3D representation by adding
#' a z-axis (typically \-log10(p) or similar association strength).
#'
#' Features:
#'
#' - caps extreme values at `zmax`
#' - colors cis/trans associations
#' - provides interactive hover labels
#'
#' @param d Data frame in `qtl2dplot()`-compatible format.
#' @param chrlen Chromosome lengths (e.g. `gap::hg19`, `hg38`).
#' @param zmax Maximum value for z-axis truncation.
#' @param qtl.id Prefix for SNP identifier.
#' @param qtl.prefix Prefix for QTL coordinate label.
#' @param qtl.gene Prefix for gene coordinate label.
#' @param target.type Label for target entity type.
#' @param TSS Logical; use transcription start site if TRUE.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param ... Additional arguments passed to `qtl2dplot()`.
#'
#' @return A `plotly` 3D scatter object.
#' @examples
#' \dontrun{
#' suppressMessages(library(dplyr))
#' INF <- Sys.getenv("INF")
#' d <- read.csv(
#'   file.path(INF, "work", "INF1.merge.cis.vs.trans"),
#'   as.is = TRUE
#' ) %>%
#'   dplyr::mutate(log10p = -log10p)
#' r <- qtl3dplotly(d, zmax = 300)
#' htmlwidgets::saveWidget(r, file = "INF1.qtl3dplotly.html")
#' r
#' }
#'
#' @export
#'
qtl3dplotly <- function(
    d,
    chrlen = gap::hg19,
    zmax = 300,
    qtl.id = "SNPid:",
    qtl.prefix = "QTL:",
    qtl.gene = "Gene:",
    target.type = "Protein",
    TSS = FALSE,
    xlab = "QTL position",
    ylab = "Gene position",
    ...
) {

    # base projection
    t2d <- qtl2dplot(
        d,
        chrlen,
        TSS = TSS,
        plot = FALSE,
        ...
    )

    names(t2d$data)[names(t2d$data) == "id"] <- "t2d_id"
    t2d_pos <- t2d$data

    # labels
    t2d_pos$snpid <- paste(qtl.id, t2d_pos$t2d_id)

    t2d_pos$pos_qtl <- paste0(
        qtl.prefix, t2d_pos$chr1, ":", t2d_pos$pos1
    )

    t2d_pos$pos_gene <- paste0(
        qtl.gene, t2d_pos$chr2, ":", t2d_pos$pos2
    )

    t2d_pos$target_gene <- paste0(
        target.type,
        " (gene): ",
        t2d_pos$target,
        " (",
        t2d_pos$gene,
        ")"
    )

    t2d_pos$lp <- paste("value:", t2d_pos$value)

    t2d_pos$text <- paste(
        t2d_pos$snpid,
        t2d_pos$pos_qtl,
        t2d_pos$pos_gene,
        t2d_pos$target_gene,
        t2d_pos$lp,
        sep = "\n"
    )

    # cap z values
    t2d_pos$z <- ifelse(t2d_pos$value <= zmax, t2d_pos$value, zmax)

    # axis scaffolding (reuse CM)
    n <- t2d$n
    CM <- t2d$CM

    tkvals <- tktxts <- vector()

    for (x in seq_len(n)) {
        tkvals[x] <- ifelse(
            x == 1,
            CM[x] / 2,
            (CM[x] + CM[x - 1]) / 2
        )
        tktxts[x] <- xy(x)
    }

    # plotly axes
    xaxis <- list(
        title = xlab,
        tickmode = "array",
        tickvals = tkvals,
        ticktext = tktxts
    )

    yaxis <- list(
        title = ylab,
        tickmode = "array",
        tickvals = tkvals,
        ticktext = tktxts
    )

    zaxis <- list(
        title = "-log10(p)"
    )

    # colors
    cols <- c("cis" = "#BF382A", "trans" = "#0C4B8E")

    fig <- plotly::plot_ly(
        data = t2d_pos,
        x = ~x,
        y = ~y,
        z = ~z,
        color = ~cistrans,
        colors = cols,
        type = "scatter3d",
        mode = "markers",
        text = ~text,
        hoverinfo = "text"
    ) %>%
        plotly::layout(
            scene = list(
                xaxis = xaxis,
                yaxis = yaxis,
                zaxis = zaxis,
                aspectratio = list(x = 0.9, y = 1, z = 0.6)
            ),
            title = "3D QTL cis/trans plot",
            showlegend = TRUE
        )

    fig
}
