\name{control.pointer}
\alias{control.pointer}
\title{
Control routine for pointer
}

\description{
Set various parameters for function pointer
}

\usage{control.pointer(sexlink=FALSE, split=FALSE, 
                       mating.type='00,01,02,10,11,12,20,21,22',
                       pointer.selection='0,4,5,6', 
                       pointer.degree='1,2,3,4,5', 
                       ascertainment='C')}

\arguments{
\item{sexlink}{TRUE if X-linked data are used}
\item{split}{To allow for any combination of mating type, pointer, pointer degree, and
ascertainment. The intersection of te criteria is put in one file, the complement in another.
They are analysed separately, giving a heterogeneity test.}
\item{mating.type}{set(s) of mating types listed above}
\item{pointer.selection}{set(s) of pointers to be selected}
\item{pointer.degree}{pointer degree(s) to be selected}
\item{ascertainment}{'C' for compelete ascertainment, 'M' for multiple ascertainment}
}

\value{
a list with values specified. 
}

\seealso{
\code{\link[pointer]{pointer}}
}

\examples{\dontrun{
# Sibships with one or more unclassified parent and no pointer in one file,
# all other sibships in the other
control.pointer(mating.type="09,90,99",pointer.selection="0")
}}

\keyword{models}
