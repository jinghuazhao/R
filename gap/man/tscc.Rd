\name{tscc}
\alias{tscc}
\title{Power calculation for two-stage case-control design}
\usage{tscc(model,GRR,p1,n1,n2,M,alpha.genome,pi.samples,pi.markers,K)}
\description{
This function gives power estimates for two-stage case-control design for genetic association.

The false positive rates are calculated as follows,

\deqn{P(|z1|>C1)P(|z2|>C2,sign(z1)=sign(z2))} and \deqn{P(|z1|>C1)P(|zj|>Cj||z1|>C1)}
for replication-based and joint analyses, respectively; where C1, C2, and Cj
are threshoulds at stages 1, 2 replication and joint analysis, 

\deqn{ z1 = z(p1,p2,n1,n2,pi.samples)}
\deqn{z2 = z(p1,p2,n1,n2,1-pi.samples)}
\deqn{zj = sqrt(pi.samples)*z1+sqrt(1-pi.samples)*z2}
}

\arguments{
  \item{model}{any in c("multiplicative","additive","dominant","recessive")}
  \item{GRR}{genotype relative risk}
  \item{p1}{the estimated risk allele frequency in cases}
  \item{n1}{total number of cases}
  \item{n2}{total number of controls}
  \item{M}{total number of markers}
  \item{alpha.genome}{false positive rate at genome level}
  \item{pi.samples}{sample\% to be genotyped at stage 1}
  \item{pi.markers}{markers\% to be selected (also used as the false positive rate at stage 1)}
  \item{K}{the population prevalence}
}

\value{
The returned value is a list containing a copy of the input plus output as follows,
  \item{model}{any in c("multiplicative","additive","dominant","recessive")}
  \item{GRR}{genotype relative risk}
  \item{p1}{the estimated risk allele frequency in cases}
  \item{pprime}{expected risk allele frequency in cases}
  \item{p}{expected risk allele frequency in controls}
  \item{n1}{total number of cases}
  \item{n2}{total number of controls}
  \item{M}{total number of markers}
  \item{alpha.genome}{false positive rate at genome level}
  \item{pi.samples}{sample\% to be genotyped at stage 1}
  \item{pi.markers}{markers\% to be selected (also used as the false positive rate at stage 1)}
  \item{K}{the population prevalence}
  \item{C}{threshoulds for no stage, stage 1, stage 2, joint analysis}
  \item{power}{power corresponding to C}
}

\section{References}{
Skol AD, Scott LJ, Abecasis GR, Boehkne M (2006).
Joint analysis in more efficient than replication-based aalysis for two-stage
genome-wide association studies. Nature Genetics 38:209-213
}

\examples{
K <- 0.1
p1 <- 0.4
n1 <- 1000
n2 <- 1000 
M <- 300000
alpha.genome <- 0.05
GRR <- 1.4
p1 <- 0.4
pi.samples <- 0.2
pi.markers <- 0.1

options(echo=FALSE)
cat("sample\%,marker\%,GRR,(thresholds x 4)(power estimates x 4)\n")
for(GRR in c(1.3,1.35,1.40)) {
   cat("\n")
   for(pi.samples in c(1.0,0.5,0.4,0.3,0.2)) {
      if(pi.samples==1.0) s <- 1.0
      else s <- c(0.1,0.05,0.01)
      for(pi.markers in s)
      {
        x <- tscc("multiplicative",GRR,p1,n1,n2,M,alpha.genome,
                  pi.samples,pi.markers,K)
        l <- c(pi.samples,pi.markers,GRR,x$C,x$power)
        l <- sprintf("\%.2f \%.2f \%.2f, \%.2f \%.2f \%.2f \%.2f, \%.2f \%.2f \%.2f \%.2f",
                     l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],l[10],l[11])
        cat(l,"\n")
      }
      cat("\n")
   }
}
options(echo=TRUE)
}
\author{Jing Hua Zhao}
\note{solve.skol is adapted from CaTS}
\keyword{misc}
