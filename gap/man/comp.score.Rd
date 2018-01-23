\name{comp.score}
\alias{comp.score}
\title{score statistics for testing genetic linkage of quantitative trait}
\usage{comp.score(ibddata="ibd_dist.out", phenotype="pheno.dat", mean=0,
                  var=1, h2=0.3)}
\description{
The function empirically estimate the variance of the score functions.
The variance-covariance matrix consists of two parts: the additive
part and the part for the individual-specific environmental effect.
Other reasonable decompositions are possible. 

This program has the following improvement over "score.r":

   1. It works with selected nuclear families

   2. Trait data on parents (one parent or two parents), if available,
      are utilized.

   3. Besides a statistic assuming no locus-specific dominance effect,
      it also computes a statistic that allows for such effect.
      It computes two statistics instead of one.  

Function "merge" is used to merge the IBD data for a pair with the
transformed trait data (i.e., \eqn{w_kw_l}).
}

\arguments{
\item{ibddata}{The output file from GENEHUNTER using command "dump ibd".
The default file name is \eqn{ibd_dist.out}.}

\item{phenotype}{The file of pedigree structure and trait value.  
The default file name is "pheno.dat". Columns (no headings) are: 
family ID, person ID, father ID, mother ID, gender, trait value,
where Family ID and person ID must be numbers, not characters.
Use character "NA" for missing phenotypes.}

\item{mean}{(population) mean of the trait, with a default value of 0}

\item{var}{(population) variance of the trait, with a default value of 1}

\item{h2}{heritability of the trait, with a default value of 0.3}
}     

\value{
a matrix with each row containing the location and the statistics and their p-values.
}

\section{References}{
Kruglyak L, Daly MJ, Reeve-Daly MP, Lander ES (1996) Parametric and Nonparametric 
linkage analysis: a unified multipoint approach. Am J Hum Genet 58:1347-1363

Kruglyak L, Lander ES (1998) Faster multipoint linkage analysis using Fourier transforms
J Comp Bio 1998 5:1-7

Wang K (2005) A likelihood approach for quantitative-trait-locus mapping with
selected pedigrees. Biometrics 61:465-473
}

\examples{
\dontrun{
# An example based on GENEHUNTER version 2.1, with quantitative trait data in file
# "pheno.dat" generated from the  standard normal distribution. The following
# exmaple shows that it is possible to automatically call GENEHUNTER using R
# function "system".

cwd <- getwd()
cs.dir <- file.path(path.package("gap"),"tests/comp.score")
setwd(cs.dir)
dir()
# system("gh < gh.inp")
cs.default <- comp.score()
setwd(cwd)
}
}
\author{Yingwei Peng, Kai Wang}
\note{Adapt from score2.r}
\keyword{htest}
