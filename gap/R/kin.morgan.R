#' kinship matrix for simple pedigree
#'
#' kinship matrix according to Morgan v2.1.
#'
#' @param ped individual's id, father's id and mother's id.
#' @param verbose an option to print out the original pedigree.
#'
#' @export
#' @return The returned value is a list containing:
#' \describe{
#' \item{kin}{the kinship matrix in vector form}
#' \item{kin.matrix}{the kinship matrix}
#' }
#'
#' @references
#' Morgan V2.1 \url{https://sites.stat.washington.edu/thompson/Genepi/MORGAN/Morgan.shtml}
#'
#' @seealso \code{\link[gap]{gif}}
#'
#' @examples
#' \dontrun{
#' # Werner syndrome pedigree
#' werner<-c(
#'  1, 0,  0,  1,
#'  2, 0,  0,  2,
#'  3, 0,  0,  2,
#'  4, 1,  2,  1,
#'  5, 0,  0,  1,
#'  6, 1,  2,  2,
#'  7, 1,  2,  2,
#'  8, 0,  0,  1,
#'  9, 4,  3,  2,
#' 10, 5,  6,  1,
#' 11, 5,  6,  2,
#' 12, 8,  7,  1,
#' 13,10,  9,  2,
#' 14,12, 11,  1,
#' 15,14, 13,  1)
#' werner<-t(matrix(werner,nrow=4))
#' kin.morgan(werner[,1:3])
#' }
#'
#' @author Morgan development team, Jing Hua Zhao
#' @note The input data is required to be sorted so that parents preceed their children.
#' @keywords datagen

kin.morgan<-function(ped,verbose=FALSE)
{
   v2k <- function (pars) 
   # 5/6/2004
   # this function is in spirit similar to g2a and make.del(mvnmle)
   {
       k <- floor((-1 + sqrt(1 + 8 * length(pars)))/2)
       mymatrix <- diag(1:k)
       if (k > 1) {
           for (i in 1:k) {
               mymatrix[i, 1:i] <- mymatrix[1:i,i] <- pars[1:i]
               pars <- pars[-(1:i)]
           }
       }
       mymatrix
   }
   pedsize<-dim(ped)[1]
   kin<-rep(0,pedsize*(pedsize+1)/2)
   id <- ped[,1]
   father <- ped[,2]
   mother <- ped[,3]
   tid <- c(id,father,mother)
   uid <- unique(tid[tid!=0])
   iid <- match(id,uid)
   fid <- match(father,uid,nomatch=0)
   mid <- match(mother,uid,nomatch=0)
   peddata <- rbind(id,father,mother)
   pedindex <- rbind(iid,fid,mid)
   if (verbose)
   {
      cat("The original pedigree IDs and their indices:\n")
      print(t(rbind(peddata,pedindex)))
   }
   z<-.C("kin_morgan",data=as.integer(peddata),pedsize=as.integer(pedsize),
         pedinex=as.integer(pedindex),kin=as.double(array(kin)),PACKAGE="gap")
   kin.matrix=v2k(z$kin)
   colnames(kin.matrix) <- rownames(kin.matrix) <- id
   list(kin=z$kin,kin.matrix=kin.matrix)
}
