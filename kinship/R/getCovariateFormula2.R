# $Id: getCovariateFormula2.s,v 1.1 2003/08/09 21:48:42 Therneau Exp $
#
# Updated code for getCovariateFormula
#  The default S code assumes that | would be the top of the list.  That
#  ain't necessarily so, if there were parenthesis.  Look for the bar.
# This code is identical to getGroupsFormula, except for the final 3 lines.

getCovariateFormula2 <- function(object) {
    form <- formula(object)
    if(!inherits(form, "formula")) {
        stop("\"Form\" argument must be a formula")
        }
    temp <- form[[length(form)]]

    tfun <- function(x) {
        # This function walks the parse tree looking for a vertical bar
        # When it finds one, it splits the formula there, into two peices
        if (mode(x) == '(') Recall(x[[2]])  #step over parentheses
        if (length(x) ==3) {
            if (x[[1]]== as.name('|')) list(x[[2]], x[[3]])
            else {
                temp <- Recall(x[[2]])
                if (length(temp)==2) {
                    x[[2]] <- temp[[2]]
                    list(temp[[1]], x)
                    }
                else {
                    temp <- Recall(x[[3]]) #no | to the left
                    if (length(temp)==2) {
                        x[[3]] <- temp[[1]]
                        list(x, temp[[2]])
                        }
                    else 0
                    }
                }
            }
        }
    
    temp <- tfun(temp)
    if (length(temp)==2)  form[[2]] <- temp[[1]]  # found a | character
    form
    }
