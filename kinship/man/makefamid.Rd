\name{makefamid}
\alias{makefamid}
\title{
  Identify family groups
}
\description{
Given a set of parentage relationships, this subdivides a set of
subjects into families. 
}
\usage{
makefamid(id, father.id, mother.id)
}
\arguments{
\item{id}{
a vector of unique subject identifiers
}
\item{father.id }{
for each subject, the identifier of their biolgical father
}
\item{mother.id }{
for each subject, the identifier of thier biological mother
}
}
\value{
a vector of family identifiers.
Individuals who are not blood relatives of anyone else in the data
set as assigned a family id of 0.
}
\details{
This function may be useful to create a family identifier if none
exists in the data (rare),
to check for anomalies in a given family identifier
(see the \code{familycheck} function),
or to create a more space and time efficient kinship matrix by
separating out marry-ins without children as 'unrelated'.
}
\seealso{
\code{\link{makefamid}, \link{kinship}, \link{makekinship}}
}
\examples{
\dontrun{
> newid <- makefamid(cdata$gid, cdata$dadid, cdata$momid)
> table(newid==0)
 FALSE TRUE 
 17859 8191
# So nearly 1/3 of the individuals are not blood relatives.

> kin1 <- makekinship(cdata$famid, cdata$gid, cdata$dadid, cdata$momid)
> kin2 <- makekinship(newid, cdata$gid, cdata$dadid, cdata$momid, unique=0)
> dim(kin2)
[1] 26050 26050
> dim(kin1)
[1] 26050 26050

> length(kin2@blocks)/length(kin1@blocks)
[1] 0.542462
# Basing kin1 on newid rather than cdata$famid (where marry-ins were each
#   labeled as members of one of the 426 families) reduced its size by just
#   less than half.
}
}
\keyword{manip}
% docclass is function
% Converted by Sd2Rd version 1.21.
