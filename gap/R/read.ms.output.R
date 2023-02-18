#' A utility function to read ms output
#'
#' @param msout an ms output.
#' @param is.file a flag indicating ms output as a system file or an R object.
#' @param xpose a flag to obtain the tranposed format as it is (when TRUE).
#' @param verbose when TRUE, display on screen every 1000 for large nsam.
#' @param outfile to save the haplotypes in a tab-delimited ASCII file.
#' @param outfileonly to reset gametes to NA when nsam/nreps is very large and is useful with outfile.
#'
#' @details
#' This function reads in the output of the program ms, a program to generate
#' samples under a variety of neutral models.
#'
#' The argument indicates either a file name or a vector of character strings,
#' one string for each line of the output of ms. As with the second case, it
#' is appropriate with system(,intern=TRUE), see example below.
#'
#' @export
#' @return
#' The returned value is a list storing the results:
#' - call system call to ms.
#' - seed random number seed to ms.
#' - nsam number of copies of the locus in each sample.
#' - nreps the number of independent samples to generate.
#' - segsites a vector of the numbers of segregating sites.
#' - times vectors of time to most recent ancester (TMRCA) and total tree lengths.
#' - positions positions of polymorphic sites on a scale of (0,1).
#' - gametes a list of haplotype arrays.
#' - probs the probability of the specified number of segregating sites 
#' given the genealogical history of the sample and the value to -t option.
#'
#' @references
#' \insertRef{hudson02}{gap}
#'
#' Press WH, SA Teukolsky, WT Vetterling, BP Flannery (1992). Numerical Recipes in C. Cambridge University Press, Cambridge.
#'
#' @examples
#' \dontrun{
#' # Assuming ms is on the path
#'
#' system("ms 5 4 -s 5 > ms.out")
#' msout1 <- read.ms.output("ms.out")
#'
#' system("ms 50 4 -s 5 > ms.out")
#' msout2 <- read.ms.output("ms.out",outfile="out",outfileonly=TRUE)
#'
#' msout <- system("ms 5 4 -s 5 -L", intern=TRUE)
#' msout3 <- read.ms.output(msout,FALSE)
#' }
#'
#' @author D Davison, RR Hudson, JH Zhao
#' @keywords utilities

read.ms.output <- function(msout, is.file=TRUE, xpose=TRUE, verbose=TRUE, outfile=NULL, outfileonly=FALSE)
{
    if (is.file) msout <- scan(file=msout, what=character(0), sep="\n", quiet=TRUE)
    if (is.na(msout[1])) stop("Usage: read.ms.output(msout), or read.ms.output(filename)")
    nsam <- as.integer(strsplit(msout[1], split=" ")[[1]][2])
    ndraws <- as.integer(strsplit(msout[1], split=" ")[[1]][3])
    result <- gamlist <- positions <- list()
    marker <- grep("prob",msout)
    probs <- sapply(strsplit(msout[marker], split=":"), function(vec) as.numeric(vec[2]))
    marker <- grep("time",msout)
    times <- sapply(strsplit(msout[marker], split="\t"), function(vec){ as.numeric(vec[2:3])})
    marker <- grep("segsites", msout)
    stopifnot(length(marker) == ndraws)
    segsites <- sapply(strsplit(msout[marker], split=" "), function(vec) as.integer(vec[2]))
    if (!is.null(outfile)) of <- file(outfile,"w")
    for (draw in seq(along=marker))
    {
        if (verbose) if (!(draw %% 1000)) cat(draw, " ")
        if (segsites[draw] > 0)
        {
            tpos <- strsplit(msout[marker[draw]+1], split=" ")
            positions[[draw]] <- as.numeric(tpos[[1]][2:(segsites[draw]+1)]) 
            haplotypes <- msout[(marker[draw] + 2):(marker[draw] + 2 + nsam - 1)]
            if (!is.null(outfile)) cat(paste(draw,1:nsam,haplotypes,sep="\t"),file=of,sep="\n")
            haplotypes <- strsplit(haplotypes, split="")
            h <- sapply(haplotypes, function(el) c(as.integer(el)))
            if(segsites[draw] == 1) h <- t(as.matrix(h))
        }
        else {
            h <- matrix(nrow=0, ncol=nsam)
            positions[[draw]]<- NA
        }
        if (xpose)
        {
           s <- h
           colnames(s) <- 1:nsam
        }
        else {
           s <- t(h)
           row.names(s) <- 1:nsam
        }
        if (outfileonly) gamlist[[draw]] <- NA else gamlist[[draw]] <- s
        stopifnot(all(dim(h) == c(segsites[draw],nsam))) 
    }
    if (verbose) cat("\n")
    if (!is.null(outfile)) close(of)
    z <- list(call=msout[1], seed=as.numeric(strsplit(msout[2]," ")[[1]]), nsam=nsam, nreps=ndraws,
              segsites=segsites, probs=probs, times=t(times), positions=positions, gametes=gamlist)
    invisible(z)
}
