# $Id: check.hint.s,v 1.1 2007/12/12 23:13:24 therneau Exp $
# Remove any inconsistencies in the spousal hints.  These arise in autohint
#  with complex pedigrees.  One can have ABA (subject A is on both the
#  left and the right of B), cycles, etc.  We fix this by making all of the
#  hints be consistent pairs: if A is the spouse on the left of B, then B
#  is the spouse on the right of A.  
#  If there are triple, quad, etc marriages, I may have listed a male
#    as the spouse hint of a male or female as female (when two males are
#     side by side with a line between them, or two females).
# This is called both by autohint and align.pedigree.  Why?  Because users
#   can introduce these problems too, when they modify the hints.
check.hint <- function(hints, sex) {
    n <- nrow(hints)

    for (i in (1:n)[hints[,2] !=0]) {
        spouse <- abs(hints[i,2])
        if (sex[i] == sex[spouse]) hints[i,2] <- 0   
        else if (hints[i,2] > 0) {
            # Find and break any spouse-to-the-right loops
            j <- spouse
            while (hints[j,2] > 0) {
                if (hints[j,2] ==i) hints[j,2] <- 0
                else j <- hints[j,2]
                }
            # "Spouse" is to the right of subject i.  Allow no one
            #  else to be joined on the left of spouse.  (Spouse might have
            #  a hint showing someone else on thier right though, in a
            #  second marriage for instance).
            if (hints[spouse,2] <= 0) hints[spouse,2] <- (-i)
            }
        else {
            # Break any spouse-to-the-left loops
            j <- spouse
            while (hints[j,2] < 0) {
                if (hints[j,2] == -i) hints[j,2] <- 0
                else j <- -hints[j,2]
                }
            if (hints[spouse,2] >= 0) hints[spouse,2] <- i
            }
        }
    hints
    }
