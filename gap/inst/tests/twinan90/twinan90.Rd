\name{twinan90}
\alias{twinan90}
\title{Classic twin models}
\usage{twinan90(mzdat,dzdat,vname='mzdz',xlamb=1,const=0,vmiss=-9,path=1,
           ped=0,nvar=1,form='((1x,a1,5x,F6.2))')}
\description{
Classic twin models

The function also allows for test for normality and Box-Cox transformation.
It further allows for output data in FISHER format.
}

\arguments{
\item{mzdat}{two columns of data for MZ twins}
\item{dzdat}{two columns of data for DZ twins}
\item{vname}{variable name}
\item{xlamb}{the Box-Cox transformation parameter}
\item{const}{constant to be added to each data value}
\item{vmiss}{missing value indicator}
\item{path}{path analysis to be conducted}
\item{ped}{if 1, to output data in FISHER format}
\item{nvar}{number of variables in the data file}
\item{form}{the Fortran format for the data file}
}

\value{
The returned values is in a list with the following components, while the MLEs are in 
two system files (.log and .out):

\item{h2}{The heritability estimate based on 2(rMZ-rDZ), where rMZ and rDZ are the intraclass
correlation coefficients for MZ and DZ twin pairs}
\item{seh2}{The standard error for the above statistic}
\item{nMZ}{The number of MZ pairs used}
\item{nDZ}{The number of DZ pairs used}
\item{rMZ}{The intraclass correlation for MZ}
\item{rDZ}{The intraclass correlation for DZ}
\item{covMZ}{The variance-covariance matrix for MZ}
\item{covDZ}{The variance-covariance matrix for DZ}
\item{vname.log}{the diagnostic file}
\item{vname.out}{the output file}
}

\section{References}{ 
Williams CJ, Christian JC, Norton JA Jr. (1992) TWINAN90: A FORTRAN program 
for conducting ANOVA-based and likelihood-based analyses of twin data. Computer Methods and 
Programs in Biomedicine 38:(2-3):167-176 

Williams CJ (1993). On the covariance between parameter estimates in models of twin 
data. Biometrics. 49(2):557-68
}

\examples{
\dontrun{
filespec <- file.path(.path.package("gap"),"tests/twinan90/mzdz.dat")
mzdz <- matrix(scan(filespec,skip=1),ncol=2,byrow=T)
mzdat <- mzdz[1:131,]
dzdat <- mzdz[132:206,]
twinan90(mzdat,dzdat,xlamb=2)

# the normality test can be conducted as follows,
ks.test(mzdat,"pnorm")
ks.test(dzdat,"pnorm")
qqnorm(mzdat)
qqnorm(dzdat)
# contact e-mail for Prof Williams
# http://www.webpages.uidaho.edu/~chrisw/
}
}
\author{Chris Williams, Jing Hua Zhao}
\note{built on newtw5.for}
\keyword{models}
