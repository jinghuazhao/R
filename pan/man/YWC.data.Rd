\name{YWC.data}
\docType{data}
\alias{YWC.data}
\title{A data frame for multivariate imputation}
\description{
The data serves as an worked example for PAN as follows,
\preformatted{
head(YWC.data)
Yi.Matrix <- as.matrix(YWC.data[,4:39])
Subj.vector <- as.vector(as.matrix(YWC.data[,2]))
YWC.data$intercept <- 1. 
Pred.Matrix <- as.matrix(YWC.data[,c(40,3)])
xcol <- 1:2
zcol <- 1
a.parameter = ncol(Yi.Matrix)
Binv.parameter = diag(a.parameter*a.parameter)
c.parameter = (zcol*ncol(Yi.Matrix))
Dinv.parameter = diag(c.parameter*c.parameter)
Prior.Distribution <- list(a = a.parameter, 
                           Binv = Binv.parameter, 
                           c = c.parameter, 
                           Dinv = Dinv.parameter)
Imputation.Result <- pan(Yi.Matrix, Subj.vector, Pred.Matrix, xcol, zcol, Prior.Distribution, seed=44896, iter=1000)
}
}
\usage{data(YWC.data)}
\format{A data frame}
\source{Steven J. Pierce}
  
\keyword{datasets}
