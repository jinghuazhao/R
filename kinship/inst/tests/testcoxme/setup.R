options(na.action='na.exclude', contrasts=c('contr.treatment', 'contr.poly'))
# attach("..")
library(kinship)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
date()

