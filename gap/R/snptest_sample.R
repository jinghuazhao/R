#' A utility to generate SNPTEST sample file
#'
#' @param data Data to be used.
#' @param sample_file Output filename.
#' @param ID_1 ID_1 as in the sample file.
#' @param ID_2 ID_2 as in the sample file.
#' @param missing Missing data column.
#' @param C Continuous variables.
#' @param D Discrete variables.
#' @param P Phenotypic variables.
#' @export
#' @return
#' Output file in SNPTEST's sample format.
#' @examples
#' \dontrun{
#' d <- data.frame(ID_1=1,ID_2=1,missing=0,PC1=1,PC2=2,D1=1,P1=10)
#' snptest_sample(d,C=paste0("PC",1:2),D=paste0("D",1:1),P=paste0("P",1:1))
#' }

snptest_sample <- function(data,sample_file="snptest.sample",ID_1="ID_1",ID_2="ID_2",missing="missing",C=NULL,D=NULL,P=NULL)
{
  cat(ID_1,ID_2,missing,C,D,P,file=sample_file)
  cat("\n",file=sample_file,append=TRUE)
  len_C <- length(C)
  len_D <- length(D)
  len_P <- length(P)
  cat("0 0 0",rep("C",len_C),rep("D",len_D),rep("P",len_P),file=sample_file,append=TRUE)
  cat("\n",file=sample_file,append=TRUE)
  write.table(data[c(ID_1,ID_2,missing,C,D,P)],file=sample_file,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
}
