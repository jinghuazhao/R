#
# Similar to test1.s, but with the rats data.  This has enough groups to
#  force sparse matrix computations.
#
contr.none <- function(n,contrasts=T) {
        if(is.numeric(n) && length(n) == 1.)
                levs <- 1.:n
        else {
                levs <- n
                n <- length(n)
        }
        contr <- array(0., c(n, n), list(levs, levs))
        contr[seq(1., n^2., n + 1.)] <- 1.
	contr
	}
options(contrasts=c('contr.none', 'contr.poly'))

theta <- pi/2
fit0 <- coxme(Surv(time, status) ~ rx, data=rats,
	      random= ~1|litter, variance=theta, iter=0)
tfit <- coxph(Surv(time, status) ~ factor(litter) + rx,
              data=rats, x=T, iter=0)
dt0 <- coxph.detail(tfit)

aeq(apply(dt0$score,2,sum), fit0$u)
h0 <- apply(dt0$imat,1:2,sum) + diag(c(rep(1/theta, 50),0))
h0[1:50,1:50] <- diag(diag(h0)[1:50])
aeq(as.matrix(gchol(h0)), as.matrix(fit0$hmat))
aeq(diag(gchol(h0)), diag(fit0$hmat))
aeq(diag(fit0$var), diag(solve(h0)))

# Now iteration 1
fit1 <- coxme(Surv(time, status) ~ rx, data=rats,
              random= ~1|litter, variance=theta, iter=1)
update0 <- solve(fit0$hmat, fit0$u)
update0[1:50] <- update0[1:50] - mean(update0[1:50])
aeq(update0, c(fit1$frail, coef(fit1)$fixed))
tfit <- coxph(Surv(time, status) ~ factor(litter) + rx,
              data=rats, x=T, iter=0,
              init=c(fit1$frail, coef(fit1)$fixed))
dt1 <- coxph.detail(tfit)

aeq(apply(dt1$score,2,sum)- c(fit1$frail, 0)/theta, fit1$u)
h1 <- apply(dt1$imat,1:2,sum) + diag(c(rep(1/theta, 50),0))
h1[1:50,1:50] <- diag(diag(h1)[1:50])
aeq(as.matrix(gchol(h1)), as.matrix(fit1$hmat))
aeq(diag(gchol(h1)), diag(fit1$hmat))
aeq(diag(fit1$var), diag(solve(h1)))


# And iteration 2
fit2 <- coxme(Surv(time, status) ~ rx, data=rats,
              random= ~1|litter, variance=theta, iter=2)

update1 <- solve(fit1$hmat, fit1$u)
update1[1:50] <- update1[1:50] - mean(update1[1:50])
aeq(update1, c(fit2$frail, coef(fit2)$fixed) -c(fit1$frail, coef(fit1)$fixed))

rm(update0, update1, h0, h1, dt0, dt1, fit0, fit1, fit2, tfit)
options(contrasts=c('contr.treatment', 'contr.poly'))

