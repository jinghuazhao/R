#' Preparing weight for GENECOUNTING
#'
#' This function is a R port of the GENECOUNTING/PREPARE program which takes
#' an array of genotyep data and collapses individuals with the same multilocus
#' genotype. This function can also be used to prepare for the genotype table in testing
#' Hardy-Weinberg equilibrium.
#'
#' @param data the multilocus genotype data for a set of individuals.
#' @param handle.miss a flag to indicate if missing data is kept, 0 = no, 1 = yes.
#' @param is.genotype a flag to indicate if the data is already in the form of genotype identifiers.
#' @param with.id a flag to indicate if the unique multilocus genotype identifier is generated.
#'
#' @export
#' @return
#' The returned value is a list containing:
#' \describe{
#' \item{cdata}{the collapsed genotype data}
#' \item{wt}{the frequency weight}
#' \item{obscom}{the observed number of combinations or genotypes}
#' \item{idsave}{optional, available only if with.id = 1}
#' }
#'
#' @references
#' Zhao JH, Sham PC (2003). Generic number system and haplotype analysis. Comp Prog Meth Biomed 70:1-9
#'
#' @seealso \code{\link[gap]{genecounting}},\code{\link[gap]{hwe.hardy}}
#'
#' @examples
#' \dontrun{
#' require(gap.datasets)
#' data(hla)
#' x <- hla[,3:8]
#'
#' # do not handle missing data
#' y<-pgc(x,handle.miss=0,with.id=1)
#' hla.gc<-genecounting(y$cdata,y$wt)
#'
#' # handle missing but with multilocus genotype identifier
#' pgc(x,handle.miss=1,with.id=1)
#'
#' # handle missing data with no identifier
#' pgc(x,handle.miss=1,with.id=0)
#' }
#'
#' @author Jing Hua Zhao
#' @note Built on pgc.c.
#' @keywords utilities

pgc <- function (data,handle.miss=1,is.genotype=0,with.id=0)
{
    nobs <- dim(data)[1]
    nloci2 <- dim(data)[2]
    if (is.genotype)
    {
       nloci <- nloci2
       data<-cbind(data,data)
       a1 <- a2 <- gid <- 0
       for (i in 1:nobs)
       {
           row.i <- data[i,]
           for (j in 1:nloci)
           {
               .C("g2a_",s=as.integer(row.i[j]),a1=as.integer(a1),a2=as.integer(a2),gid=as.integer(gid),PACKAGE="gap")
               data[i,2*j-1] <- a1
               data[i,2*j] <- a2
           }
       }
    }
    else nloci <- nloci2/2
    data <- as.matrix(data)
    stack <- rbind(data[,(2*1:nloci)-1],data[,(2*1:nloci)])
    alleles <- apply(stack,2,max)
    idsave <- wt <- array(0,nobs)
    obscom <- nobs
    data <- t(data)
    gret <- matrix(array(0,nobs*nloci2),nrow=nobs)
    z <- .C("pgc_c",gdata=as.integer(data),handlemiss=as.integer(handle.miss),nobs=as.integer(nobs),nloci=as.integer(nloci),
            alleles=as.integer(alleles), wt=as.integer(wt),gret=as.integer(gret),
            withid=as.integer(with.id),idsave=as.double(idsave),obscom=as.integer(obscom),PACKAGE="gap")
    subset <- 1:(z$obscom)
    gret <- matrix(z$gret,nrow=nloci2)[,subset]
    if (with.id) list(cdata=t(gret),obscom=z$obscom,idsave=z$idsave[subset],wt=z$wt[subset])
    else list(cdata=t(gret),obscom=z$obscom,wt=z$wt[subset])
}

# R port of GENECOUNTING/PREPARE
# 29-1-2004 start implementing
# 30-1-2004 in shape
# 31-1-2004 working
