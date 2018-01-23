\name{whscore}
\alias{whscore}
\title{Whittemore-Halpern scores for allele-sharing}
\usage{whscore(allele,type)}
\description{
Allele sharing score statistics
}

\arguments{
  \item{allele}{a matrix of alleles of affected pedigree members}
  \item{type}{0 = pairs, 1 = all}
}

\value{
The returned value is the value of score statistic

}

\section{References}{
Kruglyak L, Daly MJ, Reeve-Daly MP, Lander ES (1996) Parametric and Nonparametric 
linkage analysis: a unified multipoint approach. Am. J. Hum. Genet. 58:1347-1363

Whittemore AS, Halpern J (1994) A class of tests for linkage using affected 
pedigree members. Biometrics 50:118-127

Whittemore AS, Halpern J (1994) Probability of gene identity by descent: 
computation and applications. Biometrics 50:109-117
}

\examples{
\dontrun{
c<-matrix(c(1,1,1,2,2,2),ncol=2)
whscore(c,type=1)
whscore(c,type=2)
}
}
\author{Leonid Kruglyak, Jing Hua Zhao}
\note{adapted from GENEHUNTER}
\keyword{utilities}
