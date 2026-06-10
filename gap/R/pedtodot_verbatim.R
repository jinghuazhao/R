#' Pedigree-drawing with graphviz
#'
#' @param f A data.frame containing pedigrees, each with pedigree id, individual id, father id, mother id, sex and affection status.
#'
#' @details
#' Read a GAS or LINKAGE pedigree and write a Graphviz DOT representation for each pedigree.
#'
#' This is a verbatim translation of the original `pedtodot` implemented in Bash/awk.
#' To check independently, try `xclip < pedtodot_verbatim.R` and paste into an R session.
#'
#' @export
#' @return
#' For each pedigree, a Graphviz DOT file is written to the working directory.
#'
#' @examples
#' \dontrun{
#' # pedigree p3 in pedtodot
#'   pedtodot_verbatim(p3)
#' }
#' @note Adapted from Bash/awk script by David Duffy

pedtodot_verbatim <- function(f)
{
  shape <- shade <- character()
  shape["f"] <- "box,regular=1"
  shape["1"] <- "box,regular=1"
  shape["m"] <- "circle"
  shape["2"] <- "circle"
  shape["u"] <- "diamond"
  shape["0"] <- "diamond"
  shade["y"] <- "grey"
  shade["2"] <- "grey"
  shade["n"] <- "white"
  shade["1"] <- "white"
  shade["x"] <- "white"
  shade["0"] <- "white"
  eol <- "\n"
  pid <- f[[1]]
  uid <- unique(pid)
  for (pidx in seq_along(uid))
  {
    p <- uid[pidx]
    cat("Pedigree ",p,"\n")
    ped <-f[f[[1]]==p, drop=FALSE]
    if (nrow(ped) == 0) next
    sink(paste0(p, ".dot"))
    sex <- aff <- character()
    marriage <- list()
    child <- list()
    for (i in seq_len(nrow(ped)))
    {
       f2 <- ped[i,2]; f3 <- ped[i,3]; f4 <- ped[i,4]; f5 <- ped[i,5]; f6 <- ped[i,6]
       f2 <- as.character(f2)
       f3 <- as.character(f3)
       f4 <- as.character(f4)
       sex[f2] <- f5
       aff[f2] <- "x"
       if (grepl("^[012nyx]$", f6)) aff[f2] <- f6
       if (!is.na(f3) && !is.na(f4) && !f3 %in% c("x","0",".","")) {
          parents <- paste0(f3,"-",f4)
          n <- marriage[[parents]]
          if (is.null(n)) n <- 0
          n <- n + 1
          marriage[[parents]] <- n
          child[paste0(parents,"-",n)] <- f2
       }
    }
    cat (paste0("digraph Ped_", p, " {", eol))
    cat ("node [shape=diamond] ;",eol)
    cat ("ratio =\"auto\" ;",eol)
    cat ("mincross = 2.0 ;",eol)
    cat (paste0("label=\"Pedigree ", p, "\" ;",eol))
    cat ("rotate=0 ;",eol)
    for (id in names(sex)) {
        sx <- sex[id]
        if (is.null(sx) || length(sx)==0 || is.na(sx)) sx <- "u"
        af <- aff[id]
        if (is.null(af) || length(af)==0 || is.na(af)) af <- "x"
        sh <- shape[sx]; if (is.na(sh)) sh <- "ellipse"
        sd <- shade[af]; if (is.na(sd)) sd <- "white"
        cat(paste0("\"", id, "\" [shape=", sh, ", style=filled,fillcolor=", sd, "] ;", eol))
    }
    for (pair in sort(names(marriage))) {
        par <- strsplit(pair, "-", fixed=TRUE)[[1]]
        mating_t <- paste0("\"t_", par[1], "x", par[2], "\"")
        mating_b <- paste0("\"b_", par[1], "x", par[2], "\"")
        cat(paste0(mating_t, " [shape=diamond,style=filled,label=\"\",height=.1,width=.1] ;", eol))
        cat(paste0(mating_b, " [shape=diamond,style=filled,label=\"\",height=.1,width=.1] ;", eol))
        cat(paste0("\"", par[1], "\" -> ", mating_t, " [dir=none, weight=1, penwidth=3.0] ;", eol))
        cat(paste0("\"", par[2], "\" -> ", mating_t, " [dir=none, weight=1, penwidth=3.0] ;", eol))
        cat(paste0(mating_t, " -> ", mating_b, " [dir=none, weight=1, penwidth=3.0] ;", eol))
        n_child <- marriage[pair]
        if (is.na(n_child)) next
        for (k in seq_len(n_child)) {
            cat(paste0(mating_b, " -> \"", child[paste0(pair, "-", k)], "\" [dir=none, weight=2] ;", eol))
        }
    }
    cat ("}",eol)
    sink()
  }
}
