R2 <- function (N0, N1, K, P, y, g)
{
    N <- N0 + N1
    t <- -qnorm(K)
    z <- dnorm(t)
    m <- z/K
    m2 <- -m * K/(1 - K)

  # observed scale R2
    l <- lm(y ~ g)
    R2O <- with(l, var(fitted.values)/(N1/N * N0/N))

  # liability scale transformation R2
    theta <- m * (P - K)/(1 - K) * (m * (P - K)/(1 - K) - t)
    C <- K * (1 - K)/z^2 * K * (1 - K)/(P * (1 - P))
    R2 <- R2O * C/(1 + R2O * theta * C)

  # Cox & Snell R2
    l0 <- logLik(glm(y ~ 1, family = binomial(logit)))
    l1 <- logLik(glm(y ~ g, family = binomial(logit)))
    R2CS <- 1 - exp((l0 - l1) * (2/N))

  # Nagelkerke R2
    require(rms)
    l <- lrm(y ~ g)
    R2N <- with(l, stats["R2"])

  # liability scale AUC R2
    require(pROC)
    l <- glm(y ~ g, family = binomial(logit))
    v <- with(l, auc(y, linear.predictors))
    Q <- qnorm(v)
    R2AUC <- 2 * Q^2/((m2 - m)^2 + Q^2 * m * (m - t) + m2 * (m2 - t))

  # probit liability scale R2
    p <- glm(y ~ g, family = binomial(probit))
    R2P <- with(p, var(linear.predictors)/(var(linear.predictors) + 1))

  # logistic liability scale R2
    R2L <- with(l, var(linear.predictors)/(var(linear.predictors) + pi^2/3))

  # weighted probit model R2
    w <- (1 - P) * K/(P * (1 - K))
    wt <- y + 1
    wt[wt == 2] <- w
    p <- glm(y ~ g, weights = wt, family = binomial(probit))
    r <- runif(N, 0, 1)
    sel <- with(p, linear.predictors[y == 0 | r < w])
    R2WP <- var(sel)/(var(sel) + 1)

  # weighted logistic model R2
    l <- glm(y ~ g, weights = wt, family = binomial(logit))
    r <- runif(N, 0, 1)
    sel <- with(l, linear.predictors[y == 0 | r < w])
    R2WL <- var(sel)/(var(sel) + pi^2/3)

    invisible(
      list(R2=R2, R2CS=R2CS, R2N=R2N, R2AUC=R2AUC,
           R2P=R2P, R2L=R2L, R2WP=R2WP, R2WL=R2WL)
    )
}

## Lee SH, Goddard ME, Wray NR, Visscher PM (2012).
## A Better Coefficient of Determination for Genetic Profile Analysis.
## Genetic Epidemiology 36: 214–224
