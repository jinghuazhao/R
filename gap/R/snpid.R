#' SNP id by chr:pos+a1/a2
#'
#' @param chr Chromosome.
#' @param pos Position.
#' @param a1 Allele 1.
#' @param a2 Allele 2.
#' @param prefix Prefix of the identifier.
#' @param seps Delimiters.
#' @param uppercase A flag to return in upper case.
#'
#' @details
#' This function generates unique identifiers for variants
#'
#' @export
#' @return Identifier.
#' @examples
#' # rs12075
#' chr_pos_a1_a2(1,159175354,"A","G",prefix="chr",seps=c(":","_","_"),uppercase=TRUE)

chr_pos_a1_a2 <- function(chr,pos,a1,a2,prefix="chr",seps=c(":","_","_"),uppercase=TRUE)
{
  chr <- paste0(prefix,chr)
  chrpos <- paste(chr,pos,sep=seps[1])
  a1a2 <- paste(a1,a2,sep=seps[3])
  a2a1 <- paste(a2,a1,sep=seps[3])
  swap <- (a1 > a2)
  a1a2[swap] <- a2a1[swap]
  a1a2.lower <- tolower(a1a2)
  a1a2.upper <- toupper(a1a2)
  if(uppercase) paste(chrpos,a1a2.upper,sep=seps[2]) else paste(chrpos,a1a2.lower,sep=seps[2])
}

#' Retrieval of chr:pos+a1/a2 according to SNP id
#'
#' This function obtains information embedded in  unique identifiers.
#'
#' @param chr_pos_a1_a2 SNP id.
#' @param prefix Prefix of the identifier.
#' @param seps Delimiters of fields.
#'
#' @export
#' @return
#' A data.frame with the following variables:
#' - chr Chromosome.
#' - pos Position.
#' - a1 Allele 1.
#' - a2 Allele 2.
#'
#' @examples
#' # rs12075
#' inv_chr_pos_a1_a2("chr1:159175354_A_G",prefix="chr",seps=c(":","_","_"))

inv_chr_pos_a1_a2 <- function(chr_pos_a1_a2,prefix="chr",seps=c(":","_","_"))
{
  if ((seps[1]==seps[2])&(seps[2]==seps[3]))
  {
    s <- sapply(chr_pos_a1_a2,strsplit,seps[1])
    chr <- lapply(s,"[",1)
    pos <- lapply(s,"[",2)
    a1 <- lapply(s,"[",3)
    a2 <- lapply(s,"[",4)
  } else if ((seps[1]!=seps[2])&(seps[2]==seps[3]))
  {
    s <- sapply(chr_pos_a1_a2,strsplit,seps[2])
    chrpos <- lapply(s,"[",1)
    s1 <- sapply(chrpos,strsplit,seps[1])
    chr <- lapply(s1,"[",1)
    pos <- lapply(s1,"[",2)
    a1 <- lapply(s,"[",2)
    a2 <- lapply(s,"[",3)
  } else if ((seps[1]!=seps[2])&(seps[2]!=seps[3]))
  {
    s <- sapply(chr_pos_a1_a2,strsplit,seps[2])
    chrpos <- lapply(s,"[",1)
    s1 <- sapply(chrpos,strsplit,seps[1])
    chr <- lapply(s1,"[",1)
    pos <- lapply(s1,"[",2)
    s2 <- lapply(s,"[",2)
    s3 <- sapply(s2,strsplit,seps[3])
    a1 <- lapply(s3,"[",1)
    a2 <- lapply(s3,"[",2)
  }
  if (prefix=="") chr <- gsub("chr","",chr)
  s <- data.frame(chr=unlist(chr),pos=unlist(pos),a1=unlist(a1),a2=unlist(a2))
  names(s) <- c("chr","pos","a1","a2")
  return(s)
}
