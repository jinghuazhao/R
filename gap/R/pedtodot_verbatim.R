#' Pedigree-drawing with graphviz
#'
#' Read a GAS or LINKAGE format pedigree, return a digraph in the dot language and optionally
#' call dot/neato to make pedigree drawing.
#'
#' @param f A data.frame containing pedigrees, each with pedigree id, individual id, father id, mother id, sex and affection status.
#' @param run A flag to run dot/neato on the generated .dot file(s).
#' @param toDOT A flag to generate script for DOT::dot().
#' @param ... Other flag(s) for DOT::dot().
#'
#' @details
#' This is a verbatim translation of the original pedtodot implemneted in Bash/awk in contrast to `pedtodot` which was largely a mirror.
#' To check independently, try `xsel -i < <(cat pedtodot_verbatim.R)` or `cat pedtodot_verbatim.R | xsel -i` and paste into an R session.
#'
#' @export
#' @return
#' No value is returned but outputs in .dot, .pdf, and .svg.
#'
#' @examples
#' \dontrun{
#' # pedigree p3 in pedtodot
#'   pedtodot_verbatim(p3,run=TRUE,toDOT=TRUE,return="verbatim")
#' }
#' @note Adapted from Bash/awk script by David Duffy

pedtodot_verbatim <- function(f,run=FALSE,toDOT=FALSE,...)
{
  shape <- shade <- array()
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
  eol <- ifelse(toDOT,"","\n")

  uid <- unique(f[,1])
  for (pidx in 1:length(uid))
  {
    p <- uid[pidx]
    cat("Pedigree ",p,"\n")
    sink(paste0(p,ifelse(toDOT,".R",".dot")))
    ped <- subset(f,f[,1]==p)
    sex <- aff <- array()
    marriage <- child <- array(0)
    for (i in 1:nrow(ped))
    {
       f2 <- ped[i,2]; f3 <- ped[i,3]; f4 <- ped[i,4]; f5 <- ped[i,5]; f6 <- ped[i,6]
       sex[paste(f2)] <- paste(f5)
       aff[paste(f2)] <- "x"
       if (grepl(f6,"012nyx")) aff[f2] <- f6
       if (f3 != "x" & f3 != "0") {
          parents <- paste0(f3,"-",f4)
          if(is.na(marriage[parents])) marriage[parents] <- 1
          else marriage[parents] <- marriage[parents] + 1
          child[paste0(parents,"-",marriage[parents])] <- f2
       }
    }
    cat (paste0(ifelse(toDOT,"DOT::dot('",""),"digraph Ped_", p, " {", eol))
    cat ("node [shape=diamond] ;",eol)
    cat ("ratio =\"auto\" ;",eol)
    cat ("mincross = 2.0 ;",eol)
    cat (paste0("label=\"Pedigree ", p, "\" ;",eol))
    cat ("rotate=0 ;",eol)
    for (s in 1:(length(sex)-1)) {
        cat (paste0("\"", names(sex[s+1]), "\" [shape=", shape[sex[s+1]], ",", " style=filled,fillcolor=", shade[aff[s+1]], "] ;",eol))
    }
    for (m in 1:(length(marriage)-1)) {
        par <- unlist(strsplit(names(marriage[m+1]),"-"))
        mating_t <- paste0("\"t_", par[1], "x", par[2], "\"")
        mating_b <- paste0("\"b_", par[1], "x", par[2], "\"")
        cat (paste0(mating_t, " [shape=diamond,style=filled,label=\"\",height=.1,width=.1] ;",eol))
        cat (paste0(mating_b, " [shape=diamond,style=filled,label=\"\",height=.1,width=.1] ;",eol))
        cat (paste0("\"", par[1], "\" -> ", mating_t, " [dir=none, weight=1, penwidth=3.0] ;",eol))
        cat (paste0("\"", par[2], "\" -> ", mating_t, " [dir=none, weight=1, penwidth=3.0] ;",eol))
        cat (paste0(mating_t, " -> ", mating_b, " [dir=none, weight=1, penwidth=3.0] ;",eol))
        pair <- paste0(par[1],"-",par[2])
        for (k in 1:marriage[pair]) {
            cat (paste0(mating_b, " -> \"", child[paste0(pair,"-",k)], "\"", " [dir=none, weight=2] ;",eol))
        }
    }
    cat (paste0(ifelse(toDOT,"}')","}"),eol))
    sink()
    if (run) {
       if(!toDOT) {
         cat(sprintf("running dot/neato on %s\n",p))
         system(sprintf("dot -Tpdf %s.dot -o %s_dot.pdf",p,p))
         system(sprintf("neato -Tsvg %s.dot -o %s_neato.svg",p,p))
       } else {
         cat(sprintf("running DOT::dot on %s\n",p))
         requireNamespace("DOT")
         source(paste0(p,".R"))
       }
    }
  }
}
