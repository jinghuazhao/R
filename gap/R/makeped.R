#' A function to prepare pedigrees in post-MAKEPED format
#'
#' @param pifile input filename.
#' @param pofile output filename.
#' @param auto.select no loops in pedigrees and probands are selected automatically? 0=no, 1=yes.
#' @param with.loop input data with loops? 0=no, 1=yes.
#' @param loop.file filename containing pedigree id and an individual id for each loop, set if with.loop=1.
#' @param auto.proband probands are selected automatically? 0=no, 1=yes.
#' @param proband.file filename containing pedigree id and proband id, set if auto.proband=0 (not implemented).
#'
#' @details
#' Many computer programs for genetic data analysis requires pedigree data to be in the so-called
#' ``post-MAKEPED'' format. This function performs this translation and allows for some
#' inconsistences to be detected.
#'
#' The first four columns of the input file contains the following information:
#'
#' pedigree ID, individual ID, father's ID, mother's ID, sex
#'
#' Either father's or mother's id is set to 0 for founders, i.e. individuals with no parents.
#' Numeric coding for sex is 0=unknown, 1=male, 2=female. These can be followed by satellite
#' information such as disease phenotype and marker information.
#'
#' The output file has extra information extracted from data above.
#'
#' Before invoking makeped, input file, loop file and proband file have to be prepared.
#'
#' By default, auto.select=1, so translation proceeds without considering loops and proband statuses.
#' If there are loops in the pedigrees, then set auto.select=0, with.loop=1, loop.file="filespec".
#'
#' There may be several versions of makeped available, but their differences with this port should
#' be minor.
#'
#' @export
#' @examples
#' \dontrun{
#' cwd <- getwd()
#' cs.dir <- file.path(find.package("gap.examples"),"tests","kinship")
#' setwd(cs.dir)
#' dir()
#' makeped("ped7.pre","ped7.ped",0,1,"ped7.lop")
#' setwd(cwd)
#' # https://lab.rockefeller.edu/ott/
#' }
#'
#' @note adapted from makeped.c by W Li and others.
#' keywords datagen

makeped<-function(pifile="pedfile.pre",pofile="pedfile.ped",auto.select=1,
                  with.loop=0,loop.file=NA,auto.proband=1,proband.file=NA)
{
  z<-.C("makeped_c",pifile=as.character(pifile),pofile=as.character(pofile),autoselect=as.integer(auto.select),
        withloop=as.integer(with.loop),loopfile=as.character(loop.file),
        autoproband=as.integer(auto.proband),probandfile=as.character(proband.file),PACKAGE="gap")
}
#
# 18-10-03 start to implement
# 19-10-03 worked but could not figure out loops
# 20-10-03 polished and tested with Abbas Parsian's homozygosity mapping pedigrees, add documentation
