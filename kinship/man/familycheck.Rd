\name{familycheck}
\alias{familycheck}
\title{
  Error check for a family classification
}
\description{
Given a family id vector, also compute the familial grouping from
first principles using the parenting data, and compare the results.
}
\usage{
familycheck(famid, id, father.id, mother.id, newfam)
}
\arguments{
\item{famid}{
a vector of family identifiers
}
\item{id}{
a vector of unique subject identifiers
}
\item{father.id }{
vector containing the id of the biological father
}
\item{mother.id }{
vector containing the id of the biological mother
}
\item{newfam}{
the result of a call to \code{makefamid}.
If this has allready been computed by the user,
adding it as an argument shortens the running time somewhat.
}
}
\value{
a data frame with one row for each unique family id in the
\code{famid} argument.
Components of the output are


\item{famid}{
the family id, as entered into the data set
}
\item{n}{
number of subjects in the family
}
\item{unrelated}{
number of them that appear to be unrelated to anyone else 
in the entire pedigree set.  
This is usually marry-ins with no 
children (in the pedigree), and if so are not a problem.
}
\item{split}{
number of unique "new" family ids.
If this is 0, it means that no one in this "family" is related to
anyone else (not good); 
1 = everythings is fine;
2+= the family appears to be a set of disjoint trees.  Are you
                missing some of the people?
}
\item{join}{
number of other families that had a unique famid, but are actually
joined to this one.  0 is the hope.
If there are any joins, then an attribute "join" is attached.  It will be
a matrix with famid as row labels, new-family-id as the columns, and
the number of subjects as entries.  
}
}
\details{
The \code{makefamid} function is used to create a
de novo family id from the parentage data, and this is compared to the
family id given in the data.
}
\seealso{
\code{\link{makefamid}, \link{makekinship}}
}
\examples{
\dontrun{
#
# This is from a pedigree that had some identifier errors
#
> checkit<- familycheck(ids2$famid, ids2$gid, ids2$fatherid, ids2$motherid)
> table(checkit$split)    # should be all 1's
   0   1 2 
 112 424 4
# Shows 112 of the "families" were actually isolated individuals,
#  and that 4 of the families actually split into 2.  
#  In one case, a mistyped father id caused one child, along with his spouse
#  and children, to be "set adrift" from the connected pedigree.

> table(checkit$join)
   0 1 2 
 531 6 3
#
# There are 6 families with 1 other joined to them (3 pairs), and 3 with
#  2 others added to them (one triplet).
# For instance, a single mistyped father id of someone in family 319,
#  which was by bad luck the id of someone else in family 339,
#  was sufficient to join two groups.
> attr(checkit, 'join')
    [,1] [,2] [,3] [,4] [,5] [,6] [,7] 
 31   78    0    0    0    0    0    0
 32    3   15    0    0    0    0    0
 33    6    0   12    0    0    0    0
 63    0    0    0   63    0    0    0
 65    0    0    0   17   16    0    0
122    0    0    0    0    0   16    0
127    0    0    0    0    0   30    0
319    0    0    0    0    0    0   20
339    0    0    0    0    0    0   37
}
}
\keyword{manip}
% docclass is function
% Converted by Sd2Rd version 1.21.
