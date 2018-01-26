\name{pointer}
\alias{pointer}
\title{
Complex segregation analysis using pointers
}

\description{
This package implements mixed models of Morton \& MacLean (1974), Lalouel and Morton (1981), Morton et al. 
(1983).

We assume that a continuous variable, x, results from the independent contributions of a major locus, a
multifactorial transmissible component and random environment. Random mating is assumed. The major locus,
either autosomal or sex-linked, has two alleles yielding three genotypes, with effects and prior
probabilities given in Lalouel \& Morton (1981). Mendelian transmission is assumed, but this has been
extended. Let \eqn{c\sim N(0,C)}{c~N(0,C)}, \eqn{e\sim N(0,E)}{e~N(0,E)}, \eqn{x\sim N(u,V)}{x~N(u,V)},
with V = G + C + E, where u and G are the mean effect and variance due to the major locus, as we assume
cov(c,e)=cov(g,c)=cov(g,e)=0. When generation-specific heritabilities are assumed, we have 
\eqn{X_A = g + c_A + e_A}, \eqn{X_K = g + c_K + e_K}, \eqn{V_A = G + C_A + E_A = V_K = G + C_K + E_K = V}, 
where A and K subscripts designate adults and children respectively. Childhood heritability is denoted
\eqn{H = C_K/V} while adult heritability is \eqn{HZ = C_A/V} with \eqn{Z = C_A/C_K}.

The iterable parameters of the models are:

\tabular{ll}{
V \tab = total variance of x. Set to 1 when only affection status is defined. \cr
u \tab = mean of x. Set of 0 when only affection status is defined. \cr
d \tab = dominance \cr
q \tab = gene frequency of the allele leading to affection or elevation of x. \cr 
H \tab = childhood heritability, \eqn{H = C_K/V}. \cr
Z \tab = ratio of adulthood to childhood heritability, \eqn{Z = C_A/C_K}. \cr
x \tab = proportion of sporadic cases due to new mutants \cr
W \tab = additional random envrionmental effects relating quantitative trait x \cr
  \tab  and libability to affection y, when so defined. If affection is defined \cr
  \tab  on x, set W = 0. \cr
\eqn{\tau_1} \tab = p(AA transmits A). \cr
\eqn{\tau_2} \tab = p(Aa transmits A). \cr
\eqn{\tau_3} \tab = p(aa transmits A). \cr
R \tab = parent-offspring correlation of multifactorial transmissible components.
}

The model requires specification of other parameters estimated prior to segregation analysis: 
(1). thresholds \eqn{Z_{xi}, i=1, \ldots, n_z} or affection rates \eqn{I_i, i=1, \ldots, n_I} when
phenotype specification concerns only a quantitative trait of affection status but not both in a given
set of data; otherwise, affection rates \eqn{I_i, i=1,\ldots,n_I} but only one threshold \eqn{Z_X} for
selection on a quantitative trait can be defined, and (2) ascertainment probabilities \eqn{\pi_i, i=1,
\ldots,n_Z}
}

\usage{
pointer(datfile="poidat",jobfile="poijob",profile="poipro",
        terfile="poiter",control=control.pointer())
}

\arguments{
\item{datfile}{The data file. It contains individual records indicating family identification number,
position within family, and other relevant information. Records for one or both parents need not be given.
Only one pointer of each  type may be specified. Data must be sorted by the GR field, and by PO within ID.}
\item{jobfile}{The job file. The major control is
\preformatted{
PT(ID=f1, PO=f2, AF=f3, QU=f4, LI=f5, PT=f6, PR=f7, 
   GN=f8, GR=f9, SX=f10)
(J or C)
}
The first set of parentheses specifies input data fields. The second set determines joint (J) or 
conditional (C) probabilities, with the default being C, i.e., the probability of the phenotypes of
children conditional on phenotypes of parents and pointers. Joint probabilities of parents and children, 
conditional on phenotype of the pointer, are not applicable to selection through parents of children.

The input variables are in part optional:
\itemize{
\item ID (required) -- family identification as A or F field.
\item PO (required) -- position of an individual within a family. 1 = father, 2 = mother, 3 = child,
     4 = father's pointer, 5 = mother's pointer, and 6 = children's pointer.
\item QU (optional if AF given) -- value of a quantitative trait, which should be transformed to
      control skewness. Blank is interpreted as unknown, and is not rejected.
\item AF (optional if QU given) -- affection status, 0 = normal, 1 = affected, and 9 = unknown. \cr
Both afffection status and quantitative trait will be considered in computing probabilities 
when both are specified, unless W = 0.
\item LI (required if the affection field is specified) -- the value i in this field specifies the
individual's liability class (1,...,9)
\item PT (required only if pointers are included in data file) -- for pointer records, specifies the
degree of relationship to the pointer. Must be A2 field (table  XIII and XIV)
\item PR (required for incomplete selection) -- proband field. Selection for a given sibship is 
complete if the PR field is not defined or all children are coded 0 in this field. Selection is 
multiple incomplete if at least one child is coded \eqn{>} 0.
\item GN (optional) -- generation field. 1 = child as adult, and 2 = child as juvenile. 

If not specified, all sibships are taken as juvenile. Pointers are assigned to the pointee
generation.
\item GR (optional) -- group identification field. The group may be pedigree, index phenotype, or other
classifier used by the RA control. Data must be presorted by this field, if defined.
\item SX (required for sex-linkage) -- sex of individual. 1 = male, 2 = female.

FM, SI, SD, and TR read controls are accepted, together with special PI, LI, LZ, and ZY 
controls. The PI control has the form PI(\eqn{\pi_1},...,\eqn{\pi_{10}}), where \eqn{\pi_i}
is the value of pi applied to sibships with proband code i. The LI (required if an LI field was specified
in the major control) has the form LI(\eqn{I_1}, ..., \eqn{L_{10}}), where \eqn{I_i} is the cumulative 
incidence for liability indicator i. For sex-linkage there are two sets of parentheses, giving male and 
female incidences, respectively.

The LZ control has the form LZ(\eqn{Z_1}, ..., \eqn{Z_{10}}) for autosomes, LZ(\eqn{Z_{1M}}, ..., 
\eqn{Z_{10M}})(\eqn{Z_{1F}},...,\eqn{Z_{10F}}) for
sex-linkage, where \eqn{z_i} is the threshold corresponding to liablity indicator i.

The ZY control is of the form ZY(Y).

For all subsequent controls in the job file, individuals with values of the quantitative trait greater
than Y will be aken as affected and of unknown quantitative value. This is not reflected in the mating 
type distribution on the summary file.

Controls for analysis are PA, IT, RA, and RP. 

The PA control is of the form. 

PA(V-v, U=u, D=d, T=t, Q=q, H=h, Z=z, M=m, W=w, 
   Ti=t1, T2=t2, T3=t3, R=r) (N=n) 

where genetic parameters are given trial values in the first parentheses and the optional second set
specifies the number of integration classes from 1-11 with default 5. The default for aparameter not
specified on the PA control is 0 for U, D, T, Q, W, H, X and M, but 1 for V and Z; it is 1, 1/2, 0 
for T1, T2, T3 respectively, and 1/2 for R. The variance component W is available only for the autosomal 
case. When W \eqn{>} 0, both x and affection status are considered if both are recorded for some individuals,
assuming ascertainment through affection status only.

A PA control provides trial values for the following IT control, which in turn provides trial values
for the next IT, and so forth.

The IT control is of the form IT(\eqn{p_1},...) (H=h, T=t), where the \eqn{p_i} are genetic parameters 
to be estimated and the second set of paretheses contains parameters for numerical analysis by GEMINI, 
with defaults H=0.001, T=0.001.

The RA and RP controls giving likelihood ratios between 2 hypotheses are of the form 

(RA or RP)(V=v, U=u, T=t, Q=q, H=h, Z=z, X=x, M=m, W=w)(T) 

The parameters for the first hypothesis are taken from the preceding PA or IT control. The RA control
gives a likelihood ratio for each family, and (if a GR field was specified) for each group. An RP control
gives only the total for each group. Values of 2LN(L1/L2) outside +/- T are listed in the summary file.
The job file terminates with a CC control.
}}
\item{profile}{The prolix file}
\item{terfile}{The terse or summary file}
\item{control}{see \code{\link[pointer]{control.pointer}}}
}

\value{
At the moment, no value is returned and the outputs is read through the file.show function.
}

\seealso{
\code{\link[pointer]{control.pointer}}
}

\section{References}{
Morton NE, MacLean CJ (1974). Analysis of family resemblance: 3. complex segregation of quantitative
traits. Am J Hum Genet 26: 489-503

Lalouel JM, Morton NE (1981). Complex segregation analysis using pointers. Hum Hered 31, 312-21
}

\examples{\dontrun{
# the documentation example
pointer("poidat","poijob","poipro","poiter")
file.show("poipro")
file.show("poiter")
}}

\note{
The original documentation is given as follows.
\preformatted{

Pointer - autosomal version
Written in FORTRAN, this version for a SUN workstation.


This program is in three parts.
1. emx - a module used in splitting the data file
	 for heterogeneity testing.
2. nucfama - process a file of family data.
3. rpointa - (exec = pointr), estimate segregation parameters.


Files in this directory.

source.emx
source.nuc
sorce.rpo  - the set of source code. May be split to subroutines
	     with the SUN FORTRAN 'fsplit' command or compiled 
	     with f77 to give three exec. files (emx, nucfama,
	     pointr)

make.emx
make.nuc
make.rpo  -  individual makefiles. After 'fsplit' copy each to 
	     a file 'Makefile' and give command 'make' for 
	     compilation & linking.

pointer  -   a shell script to set up a 'cp' file that gives
	     file assignments.
poidat
poijob   -   data & job files, an example
poiter
poipro   -   output from pointer from the example
cppoi    -   an example 'cp' file.

The documentation for this program is given in 
Morton N.E., Rao D.C & Lalouel J-M (1983).
Methods in Genetic Epidemiology. Karger 
PO Box, CH-4009 Basel (Switzerland).
ISBN 3-8055-3668-2
}}
\keyword{models}
