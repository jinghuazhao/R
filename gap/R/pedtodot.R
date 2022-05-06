#' Converting pedigree(s) to dot file(s)
#'
#' This function converts GAS or LINKAGE formatted pedigree(s) into .dot file
#' for each pedigree to be used by dot in graphviz, which is a flexible package
#' for graphics freely available.
#'
#' Note that a single PostScript (PDF) file can be obtainaed by dot, fdp, 
#' or neato.
#'
#' dot -Tps <dot file> -o <ps file>  
#'
#' or
#'
#' fdp -Tps <dot file> -o <ps file>  
#'
#' or
#'
#' neato -Tps <dot file> -o <ps file>
#'
#' See relevant documentations for other formats.
#'
#' To preserve the original order of pedigree(s) in the data, you can examine the
#' examples at the end of this document.
#'
#' Under Cygwin/Linux/Unix, the PostScript file can be converted to Portable
#' Document Format (PDF) default to Acrobat.
#'
#' ps2pdf <ps file>
#'
#' Use ps2pdf12, ps2pdf13, or ps2pdf14 for appropriate versions of Acrobat
#' according to information given on the headline of <ps file>.
#'
#' Under Linux, you can also visualize the .dot file directly via command,
#'
#' dotty <dot file> &
#'
#' @param pedfile a pedigree file in GAS or LINKAGE format, note if 
#' individual's ID is character then it is necessary to specify as.is=T
#' in the read.table command.
#' @param makeped a logical variable indicating if the pedigree file is post-makeped.
#' @param sink a logical variable indicating if .dot file(s) are created.
#' @param page a string indicating the page size, e.g, A4, A5, B5, Legal, Letter, 
#' Executive, "x,y", where x, y is the customized page size.
#' @param url Unified Resource Locator (URL) associated with the diagram(s).
#' @param height the height of node(s).
#' @param width the width of node(s).
#' @param rotate if set to 90, the diagram is in landscape.
#' @param dir direction of edges, i.e., "none", "forward","back","both". This will be useful
#' if the diagram is viewed by lneato.
#'
#' @details
#' We can extract the code below (or within pedtodot.Rd) to pedtodot and then
#' use command: 
#'
#' sh pedtodot <pedigree file>
#'
#' @export
#' @return For each pedigree, the function generates a .dot file to be used by dot. The
#' collection of all pedigrees (*.dot) can also be put together.
#' 
#' @seealso package sem in CRAN and Rgraphviz in BioConductor \url{https://www.bioconductor.org/}.
#'
#' @examples
#' \dontrun{
#' # example as in R News and Bioinformatics (see also plot.pedigree in package kinship)
#' # it works from screen paste only
#' p1 <- scan(nlines=16,what=list(0,0,0,0,0,"",""))
#'  1   2   3  2  2  7/7  7/10 
#'  2   0   0  1  1  -/-  -/-  
#'  3   0   0  2  2  7/9  3/10 
#'  4   2   3  2  2  7/9  3/7  
#'  5   2   3  2  1  7/7  7/10 
#'  6   2   3  1  1  7/7  7/10 
#'  7   2   3  2  1  7/7  7/10 
#'  8   0   0  1  1  -/-  -/-  
#'  9   8   4  1  1  7/9  3/10 
#' 10   0   0  2  1  -/-  -/- 
#' 11   2  10  2  1  7/7  7/7 
#' 12   2  10  2  2  6/7  7/7 
#' 13   0   0  1  1  -/-  -/- 
#' 14  13  11  1  1  7/8  7/8 
#' 15   0   0  1  1  -/-  -/- 
#' 16  15  12  2  1  6/6  7/7 
#'
#' p2 <- as.data.frame(p1)
#' names(p2) <-c("id","fid","mid","sex","aff","GABRB1","D4S1645")
#' p3 <- data.frame(pid=10081,p2)
#' attach(p3)
#' pedtodot(p3)
#' #
#' # Three examples of pedigree-drawing
#' # assuming pre-MakePed LINKAGE file in which IDs are characters
#' pre<-read.table("pheno.pre",as.is=TRUE)[,1:6]
#' pedtodot(pre)
#' dir()      
#' # for post-MakePed LINKAGE file in which IDs are integers
#' ped <-read.table("pheno.ped")[,1:10]
#' pedtodot(ped,makeped=TRUE)
#' dir()
#' # for a single file with a list of pedigrees ordered data
#' sink("gaw14.dot")
#' pedtodot(ped,sink=FALSE)
#' sink()
#' file.show("gaw14.dot")
#' # more details
#' pedtodot(ped,sink=FALSE,page="B5",url="https://jinghuazhao.github.io/")
#'
#' # An example from Richard Mott and in the demo
#' filespec <- system.file("tests/ped.1.3.pre")
#' pre <- read.table(filespec,as.is=TRUE)
#' pre
#' pedtodot(pre,dir="forward")
#' }
#'
#' @author David Duffy, Jing Hua Zhao
#' @note This is based on the gawk script program pedtodot by David Duffy with minor changes.
#' @keywords dplot

pedtodot <- function(pedfile,makeped=FALSE,sink=TRUE,page="B5",
            url="https://jinghuazhao.github.io/",height=0.5,width=0.75,rotate=0,dir="none")
{
  if (makeped) ped <- pedfile[,-c(5,6,7,9)]
  else ped <- pedfile
  pedigree <- ped[,1]
  member <- ped[,2]
  father <- ped[,3]
  mother <- ped[,4]
  sex <- ped[,5]
  aff <- ped[,6]
  page.int <- charmatch(page,c("A4","A5","B5","Legal","Letter","Executive"))
  pagesize <- c("8.2677165,11.692913",
                "5.83,8.27",
                "7.17,10.12",
                "8.5,14",
                "8.5,11",
                "7.25,10.5") 
  ashape <- matrix(c(
             "m","box,regular=1",
             "1","box,regular=1",
             "f","circle",
             "2","circle"),ncol=2,byrow=T)
  ashade <- matrix(c(
             "y","style=filled,color=grey",
             "2","style=filled,color=grey",
             "n","style=\"setlinewidth(2)\"",
             "1","style=\"setlinewidth(2)\"",
             "x","green",
             "0","green"),ncol=2,byrow=T)
  ssize <- dim(ped)[1]
  shape <- shade <- rep('1',ssize)
  for (s in 1:ssize) {
      for (t in 1:4) if (sex[s]==ashape[t,1]) shape[s] <- ashape[t,2]
      for (t in 1:6) if (aff[s]==ashade[t,1]) shade[s] <- ashade[t,2]
  }
  uid <- unique(pedigree)
  for (j in 1:length(uid))
  {
    if(sink) cat(paste("[",uid[j],"]",sep=""))
    if(sink) sink(paste(uid[j],".dot",sep=""))
    cat(paste("digraph ped_",uid[j],sep=""),"{\n")
    if (page!="") {
       if (is.na(page.int)) cat(paste("page=\"", page, "\"",sep="")," ;\n")
       else if (page.int>0) cat(paste("page=\"", pagesize[page.int], "\"",sep="")," ;\n")
    }
    cat("ratio=\"auto\" ;\n")
    cat("mincross = 2.0 ;\n")
    cat("label=\"pedigree",uid[j],"\" ;\n")
    cat(paste("rotate=",rotate,sep="")," ;\n")
    if(url!="") cat(paste("URL=\"",url,"\"",sep="")," ;\n")
    selected <- pedigree==uid[j]
    id.j <- member[selected]
    dad.j <- father[selected]
    mom.j <- mother[selected]
    sex.j <- sex[selected]
    aff.j <- aff[selected]
  # Pedigree diagrams via kinship:
  # ped.j <- pedigree(id=id.j, dadid=dad.j, momid=mom.j, sex=sex.j, affected=aff.j)
  # plot(ped.j)
    shape.j <- shape[selected]
    shade.j <- shade[selected]
    n <- length(id.j)
    for (s in 1:n) cat(paste("\"", id.j[s], "\" [shape=", sep=""), shape.j[s], 
                      ",height=", height, ",width=",width, shade.j[s], "] ;\n")
    fid <- match(dad.j,id.j)
    mid <- match(mom.j,id.j)
    fid <- fid[!is.na(fid)]
    mid <- mid[!is.na(mid)]
    marriage <- matrix(rep(0,3*n*(n+1)/2),ncol=3)
    child <- array(rep('0',n*n*(n+1)/2+2),dim=c(n*(n+1)/2,n+2))
    k <- 1
    for (s in 1:n) {
       s1 <- fid[k]
       s2 <- mid[k]
       l <- min(s1,s2)
       u <- max(s1,s2)
       if (dad.j[s]!="x" && dad.j[s]!="0") {
          loc <- u*(u-1)/2 + l
          marriage[loc,1] <- s1
          marriage[loc,2] <- s2
          marriage[loc,3] <- marriage[loc,3] + 1
        # child[loc,1] <- s1
        # child[loc,2] <- s2
          child[loc,marriage[loc,3]+2] <- id.j[s]
          k <- k + 1
       }
    }
    marriage <- as.data.frame(marriage)
    child <- as.data.frame(child)
    married <- marriage[marriage[,3]>0,]
    n <- dim(married)[1]
    for (m in 1:n) {
        s1 <- married[m,1]
        s2 <- married[m,2]
        l <- min(s1,s2)
        u <- max(s1,s2)
        loc <- u*(u-1)/2 + l
        s1 <- id.j[s1]
        s2 <- id.j[s2]
        mating <- paste("\"", s1, "x", s2, "\"",sep="")
        cat(mating, "[shape=diamond,style=filled,label=\"\",height=.1,width=.1] ;\n")
        cat(paste("\"", s1, "\"",sep="")," -> ", mating, paste(" [dir=",dir, ",weight=1]",sep="")," ;\n")
        cat(paste("\"", s2, "\"",sep="")," -> ", mating, paste(" [dir=",dir, ",weight=1]",sep="")," ;\n")
        for (k in 1:married[m,3]) {
            cat(mating, " -> ",paste("\"", child[loc,k+2], "\"",sep=""), paste("[dir=",dir, ",weight=2]",sep="")," ;\n")
        }
    }
    cat("}\n")
    if(sink) sink()
  }
  cat("\n")
}

# History
# 02/01/2005 start experiment
# 03/01/2005 keep pedtodot in .Rd file with further work
# 04/01/2005 success with the use of 2-d array in no need of sort, put to gap
# 10/03/2008 use of setlinewidth(2) for unaffected
#
# The R/S program for GAW14 was originally written for Xiaoyan Liu's data from Stata
# the data with the following variables as required in genassoc by David Clayton
# pedigree, member, father, mother, sex, affected
# pedigrees 27, 106 have loops which can be broken via individuals 4, 5
# i.e., pedigrees 10051, 10065, individuals 10000529, 10001161
#
