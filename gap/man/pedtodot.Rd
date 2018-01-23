\name{pedtodot}
\alias{pedtodot}
\title{Converting pedigree(s) to dot file(s)}
\description{
This function converts GAS or LINKAGE formatted pedigree(s) into .dot file
for each pedigree to be used by dot in graphviz, which is a flexible package
for graphics freely available.

Note that a single PostScript (PDF) file can be obtainaed by dot, fdp, 
or neato.

dot -Tps <dot file> -o <ps file>  

or

fdp -Tps <dot file> -o <ps file>  

or

neato -Tps <dot file> -o <ps file>

See relevant documentations for other formats.

To preserve the original order of pedigree(s) in the data, you can examine the
examples at the end of this document.

Under Cygwin/Linux/Unix, the PostScript file can be converted to Portable
Document Format (PDF) default to Acrobat.

ps2pdf <ps file>

Use ps2pdf12, ps2pdf13, or ps2pdf14 for appropriate versions of Acrobat
according to information given on the headline of <ps file>.

Under Linux, you can also visualize the .dot file directly via command,

dotty <dot file> &

}

\usage{
pedtodot(pedfile,makeped=FALSE,sink=TRUE,page="B5",
         url="http://www.mrc-epid.cam.ac.uk",
         height=0.5,width=0.75,rotate=0,dir="none")
}

\arguments{
\item{pedfile}{a pedigree file in GAS or LINKAGE format, note if 
individual's ID is character then it is necessary to specify as.is=T
in the read.table command}
\item{makeped}{a logical variable indicating if the pedigree file is post-makeped}
\item{sink}{a logical variable indicating if .dot file(s) are created}
\item{page}{a string indicating the page size, e.g, A4, A5, B5, Legal, Letter, 
Executive, "x,y", where x, y is the customized page size}
\item{url}{Unified Resource Locator (URL) associated with the diagram(s)}
\item{height}{the height of node(s)}
\item{width}{the width of node(s)}
\item{rotate}{if set to 90, the diagram is in landscape}
\item{dir}{direction of edges, i.e., "none", "forward","back","both". This will be useful
if the diagram is viewed by lneato}
}

\details{
We can extract the code below (or within pedtodot.Rd) to pedtodot and then
use command: 

sh pedtodot <pedigree file>

}

\value{
For each pedigree, the function generates a .dot file to be used by dot. The
collection of all pedigrees (*.dot) can also be put together.
}

\seealso{package sem in CRAN and Rgraphviz in BioConductor \url{http://www.bioconductor.org}}

\examples{\dontrun{
# example as in R News and Bioinformatics (see also plot.pedigree in package kinship)
# it works from screen paste only
p1 <- scan(nlines=16,what=list(0,0,0,0,0,"",""))
 1   2   3  2  2  7/7  7/10 
 2   0   0  1  1  -/-  -/-  
 3   0   0  2  2  7/9  3/10 
 4   2   3  2  2  7/9  3/7  
 5   2   3  2  1  7/7  7/10 
 6   2   3  1  1  7/7  7/10 
 7   2   3  2  1  7/7  7/10 
 8   0   0  1  1  -/-  -/-  
 9   8   4  1  1  7/9  3/10 
10   0   0  2  1  -/-  -/- 
11   2  10  2  1  7/7  7/7 
12   2  10  2  2  6/7  7/7 
13   0   0  1  1  -/-  -/- 
14  13  11  1  1  7/8  7/8 
15   0   0  1  1  -/-  -/- 
16  15  12  2  1  6/6  7/7 

p2 <- as.data.frame(p1)
names(p2) <-c("id","fid","mid","sex","aff","GABRB1","D4S1645")
p3 <- data.frame(pid=10081,p2)
attach(p3)
pedtodot(p3)
#
# Three examples of pedigree-drawing
# assuming pre-MakePed LINKAGE file in which IDs are characters
pre<-read.table("pheno.pre",as.is=TRUE)[,1:6]
pedtodot(pre)
dir()      
# for post-MakePed LINKAGE file in which IDs are integers
ped <-read.table("pheno.ped")[,1:10]
pedtodot(ped,makeped=TRUE)
dir()
# for a single file with a list of pedigrees ordered data
sink("gaw14.dot")
pedtodot(ped,sink=FALSE)
sink()
file.show("gaw14.dot")
# more details
pedtodot(ped,sink=FALSE,page="B5",url="http://www.mrc-epid.cam.ac.uk/")

# An example from Richard Mott and in the demo
filespec <- system.file("tests/ped.1.3.pre")
pre <- read.table(filespec,as.is=TRUE)
pre
pedtodot(pre,dir="forward")
}
}
\author{David Duffy, Jing Hua Zhao}
\note{This is based on the gawk script program pedtodot by David Duffy with minor changes}
\keyword{dplot}
