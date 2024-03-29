% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/METAL_forestplot.R
\name{METAL_forestplot}
\alias{METAL_forestplot}
\title{forest plot as R/meta's forest for METAL outputs}
\usage{
METAL_forestplot(
  tbl,
  all,
  rsid,
  flag = "",
  package = "meta",
  method = "REML",
  split = FALSE,
  ...
)
}
\arguments{
\item{tbl}{Meta-anslysis summary statistics.}

\item{all}{statistics from all contributing studies.}

\item{rsid}{SNPID-rsid mapping file.}

\item{flag}{a variable in tbl such as cis/trans type.}

\item{package}{"meta" or "metafor" package.}

\item{method}{an explcit flag for fixed/random effects model.}

\item{split}{when TRUE, individual prot-MarkerName.pdf will be generated.}

\item{...}{Additional arguments to \code{meta::forest} or \code{metafor::forest}.}
}
\value{
It will generate a forest plot specified by pdf for direction-adjusted effect sizes.
}
\description{
forest plot as R/meta's forest for METAL outputs
}
\details{
This functions takes a meta-data from METAL (tbl) and data from contributing studies (all)
for forest plot. It also takes a SNPID-rsid mapping (rsid) as contributing studies often
involve discrepancies in rsid so it is appropriate to use SNPID, i.e., chr:pos_A1_A2 (A1<=A2).

The study-specific and total sample sizes (\code{N}) can be customised from METAL commands. By default, the input triplets each contain
a \code{MarkerName} variable which is the unique SNP identifier (e.g., chr:pos:a1:a2) and the \code{tbl} argument has variables
\code{A1} and \code{A2} as produced by METAL while the \code{all} argument has \code{EFFECT_ALLELE} and \code{REFERENCE_ALLELE} as with a \code{study} variable
indicating study name. Another variable common the \code{tbl} and \code{all} is \code{prot} variable as the function was developed in a protein
based meta-analysis. As noted above, the documentation example also has variable \code{N}.
From these all information is in place for generation of a list of forest plots through a batch run.

CUSTOMVARIABLE N\cr
LABEL N as N\cr
WEIGHTLABEL N
}
\examples{
\dontrun{
 data(OPG, package="gap.datasets")
 meta::settings.meta(method.tau="DL")
 METAL_forestplot(OPGtbl,OPGall,OPGrsid,width=8.75,height=5,digits.TE=2,digits.se=2,
                  col.diamond="black",col.inside="black",col.square="black")
 METAL_forestplot(OPGtbl,OPGall,OPGrsid,package="metafor",method="FE",xlab="Effect",
                  showweights=TRUE)
}

}
\references{
\insertRef{schwarzer07}{gap}

\insertRef{willer10}{gap}
}
\author{
Jing Hua Zhao
}
\keyword{distribution}
\keyword{hplot}
