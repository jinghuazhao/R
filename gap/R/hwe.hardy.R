#' Hardy-Weinberg equilibrium test using MCMC
#'
#' Hardy-Weinberg equilibrium test by MCMC
#'
#' @param a an array containing the genotype counts, as integer.
#' @param alleles number of allele at the locus, greater than or equal to 3, as integer.
#' @param seed pseudo-random number seed, as integer.
#' @param sample optional, parameters for MCMC containing number of chunks,
#'               size of a chunk and burn-in steps, as integer.
#'
#' @note Codes are commented for taking x a genotype object, as genotype to prepare
#' \code{a} and \code{alleles} on the fly.
#'
#' @export
#' @return
#' The returned value is a list containing:
#' \describe{
#'  \item{method}{Hardy-Weinberg equilibrium test using MCMC}
#'  \item{data.name}{name of used data if \code{x} is given}
#'  \item{p.value}{Monte Carlo p value}
#'  \item{p.value.se}{standard error of Monte Carlo p value}
#'  \item{switches}{percentage of switches (partial, full and altogether)}
#' }
#' @author Sun-Wei Guo, Jing Hua Zhao, Gregor Gorjanc
#'
#' @source https://sites.stat.washington.edu/thompson/Genepi/pangaea.shtml
#'
#' @references
#' Guo, S.-W. and E. A. Thompson (1992) Performing the exact test of
#' Hardy-Weinberg proportion for multiple alleles. Biometrics. 48:361--372.
#'
#' @seealso \code{\link[gap]{hwe}}, \code{\link[genetics]{HWE.test}}, \code{\link[genetics]{genotype}}
#'
#' @examples
#' \dontrun{
#'  # example 2 from hwe.doc:
#'    a<-c(
#'    3,
#'    4, 2,
#'    2, 2, 2,
#'    3, 3, 2, 1,
#'    0, 1, 0, 0, 0,
#'    0, 0, 0, 0, 0, 1,
#'    0, 0, 1, 0, 0, 0, 0,
#'    0, 0, 0, 2, 1, 0, 0, 0)
#'    ex2 <- hwe.hardy(a=a,alleles=8)
#'
#'    # example using HLA
#'    data(hla)
#'    x <- hla[,3:4]
#'    y <- pgc(x,handle.miss=0,with.id=1)
#'    n.alleles <- max(x,na.rm=TRUE)
#'    z <- vector("numeric",n.alleles*(n.alleles+1)/2)
#'    z[y$idsave] <- y$wt
#'    hwe.hardy(a=z,alleles=n.alleles)
#'
#'    # with use of class 'genotype'
#'    # this is to be fixed
#'    library(genetics)
#'    hlagen <- genotype(a1=x$DQR.a1, a2=x$DQR.a2,
#'                       alleles=sort(unique(c(x$DQR.a1, x$DQR.a2))))
#'    hwe.hardy(hlagen)
#'
#'    # comparison with hwe
#'    hwe(z,data.type="count")
#'
#'    # to create input file for HARDY
#'    print.tri<-function (xx,n) {
#'        cat(n,"\n")
#'        for(i in 1:n) {
#'            for(j in 1:i) {
#'                cat(xx[i,j]," ")
#'            }
#'        cat("\n")
#'        }
#'        cat("100 170 1000\n")
#'    }
#'    xx<-matrix(0,n.alleles,n.alleles)
#'    xxx<-lower.tri(xx,diag=TRUE)
#'    xx[xxx]<-z
#'    sink("z.dat")
#'    print.tri(xx,n.alleles)
#'    sink()
#'    # now call as: hwe z.dat z.out
#' }
#'
#' @note Adapted from HARDY, testable with -Dexecutable as standalone program.
#'
#' keywords htest

hwe.hardy<-function(a, alleles=3, seed=3000, sample=c(1000, 1000, 5000)) {
#   require(genetics)
#    if (!missing(x)) {
#        if (!is.genotype(x)) {
#            stop("'x' must be of class 'genotype' or 'haplotype'")
#        } else {
#            # Get genotype counts
#            tab <- table(factor(allele(x, 1), levels=allele.names(x)),
#                         factor(allele(x, 2), levels=allele.names(x)))
#            a <- as.integer(t(tab)[lower.tri(t(tab), diag=T)])
#            a <- a[order(a)]
#            # Get number of alleles
#            alleles <- length(allele.names(x))
#        }        
#    }
    if (alleles<3) stop("number of alleles should be at least 3")
    p <- 1.0
    se <- 0.0
    swp <- rep(0,3)
    z <- .C("hwe_hardy", a=as.integer(a), alleles=as.integer(alleles),
            seed=as.integer(seed), gss=as.integer(sample),
            p=as.double(p), se=as.double(se), swp=as.double(swp),
            PACKAGE="gap")
    # Printout (partly taken from htest)
    z$method <- "Hardy-Weinberg equilibrium test using MCMC"
    z$data.name <- deparse(substitute(x))
    cat("\n")
    writeLines(strwrap(z$method, prefix = "\t"))
    cat("\n")
    cat("data: ", z$data.name, "\n")
    out <- character()
    fp <- format.pval(z$p, digits=4)
    out <- c(out, paste("p-value",
                  if(substr(fp,1,1) == "<") fp else paste("=",fp)))    
    out <- c(out, paste("p-value.se", "=", format(round(z$se, 4))))
    writeLines(strwrap(paste(out, collapse = ", ")))    
    cat("percentage of switches:\n")    
    cat(paste(" - partial", "=", format(round(z$swp[1], 2))), "\n")
    cat(paste(" - full", "=", format(round(z$swp[2], 2))), "\n")
    cat(paste(" - altogether", "=", format(round(z$swp[3], 2))), "\n\n")
    RVAL <- list(method=z$method, data.name=z$data.name, p.value=z$p, 
                 p.value.se=z$se, switches=z$swp)
    names(RVAL$switches) <- c("partial", "full", "altogether")
    return(invisible(RVAL))
}
