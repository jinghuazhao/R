#' Multiple imputation analysis for hap
#'
#' This command reads outputs from hap session that uses multiple imputations, i.e. -mi# option. To
#' simplify matters it assumes -ss option is specified together with -mi option there.
#'
#' This is a very naive version of MIANALYZE, but can produce results for PROC MIANALYZE of SAS.
#'
#' It simply extracts outputs from hap.
#'
#' @param hapfile hap haplotype output file name.
#' @param assfile hap assignment output file name.
#' @param miafile mia output file name.
#' @param so to generate results according to subject order.
#' @param ns do not sort in subject order.
#' @param mi number of multiple imputations used in hap.
#' @param allsnps all loci are SNPs.
#' @param sas produce SAS data step program.
#'
#' @export
#' @return
#' The returned value is a list.
#'
#' @references
#' Zhao JH and W Qian (2003) Association analysis of unrelated individuals
#' using polymorphic genetic markers. RSS 2003, Hassalt, Belgium
#'
#' Clayton DG (2001) SNPHAP. https://github.com/chr1swallace/snphap.
#'
#' @seealso \code{\link[gap]{hap}}
#'
#' @examples
#' \dontrun{
#' # 4 SNP example, to generate hap.out and assign.out alone
#' data(fsnps)
#' hap(id=fsnps[,1],data=fsnps[,3:10],nloci=4)
#'
#' # to generate results of imputations
#' control <- hap.control(ss=1,mi=5)
#' hap(id=fsnps[,1],data=fsnps[,3:10],nloci=4,control=control)
#'
#' # to extract information from the second run above
#' mia(so=1,ns=1,mi=5)
#' file.show("mia.out")
#'
#' ## commands to check out where the output files are as follows:
#' ## Windows
#' # system("command.com")
#' ## Unix
#' # system("csh")
#' }
#'
#' @note adapted from hap, in fact cline.c and cline.h are not used.
#' keywords utilities

mia<-function(hapfile="hap.out",assfile="assign.out",miafile="mia.out",so=0,ns=0,mi=0,allsnps=0,sas=0)
{
  # to call up mi.inference here

  z<-.C("mia_c",hapfile=as.character(hapfile),assfile=as.character(assfile),
        miafile=as.character(miafile),so=as.integer(so),ns=as.integer(ns),mi=as.integer(mi),
        allsnps=as.integer(allsnps),sas=as.integer(sas),PACKAGE="gap")
}
