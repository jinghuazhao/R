# correction factor C for heritability under the liability threshold model as a function of K 

C <- function(K) {
  x <- qnorm(1-K)
  z <- dnorm(x)
  return (K*(1-K)/z^2)
}
plot(C,xlab="K",ylab=expression(C))

# correction factor Ccc can also be visualised dynamically

cc <- function(K,f) {
  x <- qnorm(1-K)
  z <- dnorm(x)
  return ((K*(1-K)/z)^2/(f*(1-f)))
}
s1 <- seq(0.01,1,by=0.01)
s2 <- seq(0.01,1,by=0.1)
library(rgl)
plot3d(s1,s2,cc(s1,s2))

