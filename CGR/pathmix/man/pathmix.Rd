\name{pathmix}
\alias{pathmix}
\title{
Path analysis of pairs of relatives
}

\description{
This package implements the path models for quantitative phenotypes such as blood
pressure as described in Morton et al. (1983). The raw data are condensed into 
estimates of familial correlations and then analysed according to models parametrised 
for nuclear families and fitted by the maximum likelihood estimates (MLEs).

Let \eqn{r_i}, \eqn{i=1,\ldots,r_m}{i=1,...,m} denote \eqn{m} observed correlations 
with corresponding sample sizes \eqn{n_i}. Let \eqn{z_i} denote Fisher's
z-transformation of \eqn{r_i=1/2\ln[(1+r_i)/(1-r_i)]}. If the covariance matrix 
between the z values is denoted by \eqn{\Sigma}, the log-likelihood of the observed
data can be written, assuming multivariate normality for all the z valuess, as
\deqn{\ln L = - \chi^2/2 + constant}{ln L = - X^2/2 + constant} where
\deqn{\chi^2 = (z - \bar z)'\Sigma^{-1}(z - \bar z)}{X^2 = (z - zbar)' Sigma^(-1)(z-zbar)}
ignoring the bias term which has a negligibale effect on the results. Here z is a
column vector of the observed z-transforms, and \eqn{\bar z}{zbar} a column vector
of z-transforms of the expected correlations \eqn{\rho_i} which are derived as functions
of path coefficients under a given model. In the likelihood function, a power of the
determinant \eqn{|\Sigma|}{|Sigma|} was incorporated into the `constant', emipirically
it is found that this quantity hardly changes for different models. Likelihood ratio
tests (LRTs) can be carried out between nested models containing different number of
parameters according to the standard asymptotic theory.

Specification of the covariance matrix \eqn{\Sigma} depends on how the data are
obtained. If \eqn{z_i} are indepedent (estimated from different samples), covariance
between \eqn{z_i} and \eqn{z_j} is zero (\eqn{i \neq j}{i<>j}). Then it is known that
asymptotic properties hold good even for small to moderate sample sizes and the 
\eqn{\chi^2}{X^2} expression simplies to 
\deqn{\chi^2 = \sum_{i=1}^mn_i(z_i-\bar z_i)^2}{X^2 = sum_i (z_i - zbar)^2} by taking
the variance of \eqn{z_i} as \eqn{1/n_i} ignoring refinements.

It is possible to perform heterogeneity testing among multiple studies under models
described above.

The path diagram (.pdf) together with the graphviz diagram (.dot) are given in doc 
directory of the installed package and inst/doc of the source package.

The path model assumes that family environment (C), or the transmissible environment,
acts additively with genotype (G) to produce a phenotype (P), all interactions being
negligible. Environmental indicies (I) are created and treated as correlated phenotypes.
Furthermore, the model allows for assortative mating through a copath (p) and social
homogamy (H) is treated in terms of simple paths (m, u) lying in the interval [-1,1].
Specific maternal effects are included by distinguishing \eqn{f_F} and \eqn{f_M} the
effects of parents' environments on that of a child reared by them, with F and M denoting
father and mother respectively. Intergenerational differences in genetic and
environmental effects are retained such that the genetic and environmental (cultural)
heritabilities in children are \eqn{h^2} and \eqn{c^2} respectively, and \eqn{h^2z^2}
and \eqn{c^2y^2} in adults. A non-transmitted sibship environment (B) is incorporated
such that the correlation between sibling environments due to non-transmitted 
environment is \eqn{b^2} in young sibs, and \eqn{b^2x^2} in adult sibs. Similarly, the
path coefficients associated with environmental indices are \eqn{i} and \eqn{iv}
respectively for young children and adults, those of genotype to indices \eqn{j} and
\eqn{jw} for young children and adults. The path coefficients from parents' genotypes
to children are assumed 1/2 ignoring age effect. The correlation between an individual's
adulthood and childhood genotype is represented by \eqn{r_G}.

Further details concerning the expected correlations, nonlinear constraints (which 
differentiate the need of optimisation routine ALMINI rather than GEMINI) can be found
in Morton et al. (1983).
}

\usage{
pathmix(iop=3, datfile="data", jobfile="job", profile="prolix", 
        terfile="summary")
}

\arguments{
\item{iop}{iop=1,2,3,4 for ALMINI test, GEMINI test, path3a, path3b}
\item{datfile}{The data file}
\item{jobfile}{The job file. The original skeleton is given as follows,

\preformatted{
PX(1)(ID=1, PO=2, PH=4, IN=3)
FM(27)(18X, F3.0, 2X, A1, 2X, 2F8.4)
TR(9)(4)(4)(0.,1.)
CC comments
EOF
XY(ALL)
PA(R11=R12=0.938)
IT(R11=R12)(h=0.01,T=0.001)
EOF
PA(H=0.85,...)
IT(H,C)(H=0.001,T=0.0001,INLC=0)
IT(...)
PA(...)
IT(...)
EOF
DATA=TONG
SUMMARY=WW1
PROLIX=WW2
EOF
inidividual correlations
EOF
Pooled observed correlations
EOF
EOT
}

The six files (marked by EOF) were required by PATHMIX for the following information:

\enumerate{
\item The first set of parentheses on the PX-control specifies estimation of 
correlations from nuclear family data and fitting path models to observed correlations.
The second set of parentheses on the PX-control specifies input data fields, 

\tabular{ll}{
ID \tab family identification numbber (F field) \cr
PO \tab position of an individual within a family \cr
   \tab (1 = father, 2 = model, 3 = young child, 4 = adult child \cr
PH \tab phenotype \cr
IN \tab index
}

The first set of FM-control specifies record length in words, and the second identifies
the input data fields. TR-control is for transformation.

This ends with a CC-control.

\item It begins with an XY-control, specifying types of correlations to be estimated.
Consider N individuals with observations X and Y, which could correspond to that of 
parents and of children, eight unknowns could be involved through multivariate normality,
which are denoted as U1, U2, V1, V2, R1, R2, R11 and R12, respectively.
Several sets of XY, PA, and IT can follow one another.

\item PA and IT controls specify the path models to be fitted. INLC is the constraint
status (0 = unconstrainted, 1 = inequality constraints, 2 = equality constraints).

\item This specifies the data file, summary output file and prolix output file.

\item A complete lisst of all individual correlations, either estimated using the XY
control or entered interatively.

\item A summary of of the observed correlations
}}
\item{profile}{The prolix file}
\item{terfile}{The terse or summary file}
}

\value{
At the moment, no value is returned and the outputs is read through the file.show function.
}

\section{References}{
Morton N.E., Rao D.C & Lalouel J-M (1983). Methods in Genetic Epidemiology. 
Karger PO Box, CH-4009 Basel (Switzerland).
}

\examples{\dontrun{
# test for ALMINI, check files ALMTEST.OUT and ALMTEST.PLX
pathmix(1)
# test for GEMINI, check files GEMTEST.OUT and GEMTEST.PLX
pathmix(2)
# path3a
pathmix(3, datfile="TEST.DAT", jobfile="test.jf")
# path3b
pathmix(4, datfile="TEST.DAT", jobfile="test.jf")
}}

\note{
The software documentation is given in the inst/doc (source package) or doc 
(the installed package) directory.

}}
\keyword{models}
