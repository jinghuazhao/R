# R^2 under the liability threshold model

t <- -qnorm(K)
z <- dnorm(t)
m <- z/K
l <- lm(y ~ g)
R2O <- with(l, var(fitted.values)/(P * (1-P))
theta <- m * (P - K)/(1 - K) * (m * (P - K)/(1 - K) - t)
C <- K * (1 - K)/z^2 * K * (1 - K)/(P * (1 - P))
R2 <- R2O * C/(1 + R2O * theta * C)

