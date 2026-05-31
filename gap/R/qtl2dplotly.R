#' 2D QTL cis/trans visualization (Plotly)
#'
#' Interactive visualization of QTL–gene relationships in 2D genomic space,
#' separating cis and trans associations into stable Plotly traces.
#'
#' This function wraps `gap::qtl2dplot()` and converts the output into an
#' interactive Plotly scatter plot.
#'
#' Unlike the original implementation, this version:
#'
#' - avoids `color = ~cistrans` (prevents row dropping in Plotly)
#' - uses explicit traces for cis/trans groups
#' - ensures stable HTML widget rendering
#' - eliminates "Ignoring observations" warnings
#'
#' @param d Data frame containing QTL data. Must include SNP, chromosome,
#' position, gene mapping, and cis/trans classification.
#' @param chrlen Chromosome length reference (default: `gap::hg19`).
#' @param qtl.id Prefix for SNP identifier in tooltip.
#' @param qtl.prefix Prefix for QTL coordinate label.
#' @param qtl.gene Prefix for gene coordinate label.
#' @param target.type Label for target entity type.
#' @param TSS Logical. If TRUE, uses transcription start site.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param ... Additional arguments passed to `gap::qtl2dplot()`.
#'
#' @return A `plotly` object.
#' @examples
#' \dontrun{
#' INF <- Sys.getenv("INF")
#' d <- read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans"),as.is=TRUE)
#' r <- qtl2dplotly(d)
#' htmlwidgets::saveWidget(r,file="INF1.qtl2dplotly.html")
#' r
#' }
#'
#' @export
#'
qtl2dplotly <- function (
    d,
    chrlen = gap::hg19,
    qtl.id = "SNPid:",
    qtl.prefix = "QTL:",
    qtl.gene = "Gene:",
    target.type = "Protein",
    TSS = FALSE,
    xlab = "QTL position",
    ylab = "Gene position",
    ...
) {

    t2d <- qtl2dplot(
        d,
        chrlen,
        TSS = TSS,
        plot = FALSE,
        ...
    )

    names(t2d$data)[names(t2d$data) == "id"] <- "t2d_id"
    t2d_pos <- t2d$data

    t2d_pos$snpid <- paste(qtl.id, t2d_pos$t2d_id)

    t2d_pos$pos_qtl <- paste0(qtl.prefix, t2d_pos$chr1, ":", t2d_pos$pos1)
    t2d_pos$pos_gene <- paste0(qtl.gene, t2d_pos$chr2, ":", t2d_pos$pos2)

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

    t2d_pos$cistrans <- as.factor(t2d_pos$cistrans)

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

    axes <- list(
        tickmode = "array",
        tickvals = tkvals,
        ticktext = tktxts
    )

    xaxis <- c(title = xlab, axes)
    yaxis <- c(title = ylab, axes)

    cols <- c("cis" = "#BF382A", "trans" = "#0C4B8E")

    fig <- plotly::plot_ly()

    for (lvl in levels(t2d_pos$cistrans)) {

        df <- t2d_pos[t2d_pos$cistrans == lvl, ]

        fig <- fig %>%
            plotly::add_trace(
                data = df,
                x = ~x,
                y = ~y,
                type = "scatter",
                mode = "markers",
                name = lvl,
                text = ~text,
                hoverinfo = "text",
                marker = list(
                    size = 11,
                    color = cols[[lvl]]
                )
            )
    }

    fig %>%
        plotly::layout(
            xaxis = xaxis,
            yaxis = yaxis,
            showlegend = TRUE
        )
}
