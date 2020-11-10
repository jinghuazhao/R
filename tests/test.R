dyn.load("pan.so")
source("pan.R")

# test : A dataset.
# Response variables are : Y1, Y2, Y3, and Y4.
# Type.t and time are predictors.
# X : design matrix of predictors

test <- read.csv('test.csv')
N <- nrow(test)
P <- ncol(test)
X <- with(test, cbind(Type.t == 2, Type.t == 3, Type.t == 4, time))

M <- 10 # number of imputations
prior <- list(a=1, Binv=1, c=1, Dinv=1) # Specified priors
PAN.y1 <- matrix(0, N, M)

# The place where errors occur #
options(echo=FALSE,width=200)
for (m in 1:M){
  cat("m =",m,"\n")
  result <- pan(test$Y1, test$ID, X, 1:4, 4, prior, seed=m, iter=100)
# print(result)
  PAN.y1[,m] <- result$y
}
colnames(PAN.y1) <- paste0("y",1:10)
options(digits=3)
PAN.y1
write.table(format(PAN.y1,digits=3),file="PAN.txt",quote=FALSE,row.names=FALSE)
save.image("test.rda")
