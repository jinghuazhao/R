# $Id: getCrossedTerms.s,v 1.1 2003/08/09 21:48:42 Therneau Exp $
# Routine to random effects formula down into separate crossed
#    terms.  
# The algorithm is pretty simple minded, and may not
#    work for everything. 
#   1. If there is only 1 vertical |, there is nothting to be done
#   2. If there are two vertical bars, find the first + sign after
#       the first bar; this is the break point.
#   3. Squawk if you get lost.
#
# This routine walks up and down the parse tree produced by Splus for a
#  formula.  By using their parse tree, I don't have to understand
#  parentheses.  However, the job would be much much easier if | had a
#  higher precedence than +.  
# splfun() does the real work.  It started simple, and grew case by case as
#  I tried harder formulas.  I'm not sure even I understand it any more.
#  S is not the language of choice for manipulating parse trees.
# 
getCrossedTerms <- function(x) {
    if (class(x)=='formula') x2 <- x[[2]]  #get the heart of the formula
        
    nbar <- function(x) {
        # count the number of bars
        if (length(x)==3)
            1*(x[[1]]==as.name('|')) + Recall(x[[2]]) + Recall(x[[3]])
        else {
            if (mode(x) == "(") Recall(x[[2]])
            else 0  #don't look inside of functions for |
            }
        }
    
    if (nbar(x2) <=1) return (list(x))  #Nothing to be done
    
    # The splitting function: find a | + combination, with optional
    #  / between bar and +, and split the formula there
    # But try never to split apart a parenthesised expression!  (Not
    #  sure I'm always successful at this).
    # lstate = what is to my left, rstate = what is to my right
    #     where 0 = unknown or other, 1=bar, 2= +
    # Success: a list of length 2 is returned if it has found a split
    #    point, or else the (farthest left, farthest right) state of the
    #    final formula.
    # There is no need to recur inside a function like log(), we only walk
    #   through infix operators and parenthesis.
    splfun <- function(x, lstate, rstate) {
        if (mode(x) == '(') {
            return(Recall(x[[2]], 0, 0))  # call it as a disjoint phrase
            }
        if (x[[1]] == as.name('|')) {
            temp1 <- Recall(x[[2]], lstate, 1) 
            if (is.list(temp1)){  #found something to the left
                x[[2]] <- temp1[[2]]
                return(list(temp1[[1]], x))
                }
            temp2 <- Recall(x[[3]], 1, rstate)
            if (is.list(temp2)) { # found something to the right
                x[[3]] <- temp2[[1]]
                return(list(x, temp2[[2]]))
                }
            return(c(temp1[1], temp2[2]))
            }
        
        if (x[[1]] == as.name('+')) {
            temp1 <- Recall(x[[2]], lstate, 2)
            if (is.list(temp1)){  #found something to the left
                x[[2]] <- temp1[[2]]
                return(list(temp1[[1]], x))
                }
            if (temp1[2]==1) # the left phrase as a whole is |
                return(list(x[[2]], x[[3]]))
            
            temp2 <- Recall(x[[3]], 2, rstate)
            if (is.list(temp2)) { # found something to the right
                x[[3]] <- temp2[[1]]
                return(list(x, temp2[[2]]))
                }
            return(c(temp1[1], temp2[2]))
            }
        
        if (x[[1]] == as.name('/') || x[[1]] == as.name('%in%')) {
            temp1 <- Recall(x[[2]], lstate, rstate)
            if (is.list(temp1)) { # found something to the left
                x[[2]] <- temp1[[2]]
                return(list(temp1[[1]], x))
                }
            temp2 <- Recall(x[[3]], lstate, rstate)
            if (is.list(temp2)){  #found something to the right
                x[[3]] <- temp2[[1]]
                return(list(x, temp2[[2]]))
                }
            return(c(temp1[1], temp2[2]))
            }

        if (length(x)==1) return(c(rstate, lstate))
        else return(c(0,0))
        }
    
    final <- splfun(x2, 0, 0)
    if (!is.list(final)) stop("Can't parse formula for random effects")

    while(1) { # Iterate, each pass breaks off one crossed term
        temp <- unlist(lapply(final, nbar))
        if (any(temp==0)) stop("Can't parse random formula for random effects")
        if (all(temp ==1)) break  # no elements left with two | chars
        
        n <- length(final)
        who <- max((1:n)[temp>1])  #choose the rightmost with two |s in it
        
        temp <- splfun(final[[who]], 0, 0)
        if (!is.list(temp)) stop("Can't parse formula for random effects")
        if (who==1) final <- c(temp[1], temp[2], final[-1])
        else if (who==n) final <- c(final[-n], temp[1], temp[2])
        else final <- c(final[1:(who-1)], temp[1], temp[2],
                        final[(who+1):n])
        }
    
    # Turn back into formulas, dropping redundant parenthesis
    lapply(final, function(x) {
        temp <- formula(~ 1)
        if (mode(x) == '(') temp[[2]] <- x[[2]]
        else                temp[[2]] <- x
        temp
        })
    }
