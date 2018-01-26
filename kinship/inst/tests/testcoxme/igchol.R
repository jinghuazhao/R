# Undo a gchol: given gchol(x), returns x
#  Mostly for debugging purposes
#
igchol <- function(x) {
    dd <- diag(x)
    ll <- as.matrix(x)
    ll %*% diag(dd) %*% t(ll)
    }
