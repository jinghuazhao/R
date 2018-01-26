#
# Now do iterations on tdata0, checking it out with coxph.detail
#   Note that the coxph calls give some "singular matrix" warnings
#   which can be ignored
#
theta <- .53
fit0 <- coxme(Surv(time, status) ~ x1 + x2, data=tdata0,
              random= ~1|grp, variance=theta, weight=wt, iter=0)

tfit <- coxph(Surv(time, status) ~ I(grp==1) + I(grp==2) + x1 + x2,
	      data=tdata0, x=T, weight=wt, iter=0)
dt0 <- coxph.detail(tfit)

aeq(apply(dt0$score,2,sum), fit0$u)
h0 <- apply(dt0$imat,1:2,sum) + diag(c(1/theta, 1/theta,0,0))
aeq(as.matrix(gchol(h0)), as.matrix(fit0$hmat))
aeq(solve(fit0$var, full=F), h0)

# Now iteration 1
fit1 <- coxme(Surv(time, status) ~ x1 + x2, data=tdata0,
              random= ~1|grp, variance=theta, weight=wt, iter=1)
aeq(fit0$u %*% fit0$var, c(fit1$frail, coef(fit1)$fixed))
tfit <- coxph(Surv(time, status) ~ I(grp==1) + I(grp==2) + x1 + x2,
	      data=tdata0, x=T, weight=wt, iter=0,
	      init=c(fit1$frail, coef(fit1)$fixed))
dt1 <- coxph.detail(tfit)

aeq(apply(dt1$score,2,sum)- c(fit1$frail, 0,0)/theta, fit1$u)
h1 <- apply(dt1$imat,1:2,sum) + diag(c(1/theta, 1/theta,0,0))
aeq(as.matrix(gchol(h1)), as.matrix(fit1$hmat))
aeq(solve(fit1$var, full=F), h1)

# And iteration 2
fit2 <- coxme(Surv(time, status) ~ x1 + x2, data=tdata0,
              random= ~1|grp, variance=theta, weight=wt, iter=2)
aeq(fit1$u %*% fit1$var, c(fit2$frail, coef(fit2)$fixed) - c(fit1$frail, coef(fit1)$fixed))


#
# Repeat the same process with the (start,stop] data set
#
fit0 <- coxme(Surv(time1, time2, status) ~ x1 + x2, data=tdata0b,
              random= ~1|grp, variance=theta, weight=wt, iter=0)

tfit <- coxph(Surv(time1, time2, status) ~ I(grp==1) + I(grp==2) + x1 + x2,
	      data=tdata0b, x=T, weight=wt, iter=0)
dt0 <- coxph.detail(tfit)

aeq(apply(dt0$score,2,sum), fit0$u)
h0 <- apply(dt0$imat,1:2,sum) + diag(c(1/theta, 1/theta,0,0))
aeq(as.matrix(gchol(h0)), as.matrix(fit0$hmat))
aeq(solve(fit0$var, full=F), h0)

# Now iteration 1
fit1 <- coxme(Surv(time1, time2, status) ~ x1 + x2, data=tdata0b,
              random= ~1|grp, variance=theta, weight=wt, iter=1)
aeq(fit0$u %*% fit0$var, c(fit1$frail, coef(fit1)$fixed))
tfit <- coxph(Surv(time1, time2, status) ~ I(grp==1) + I(grp==2) + x1 + x2,
	      data=tdata0b, x=T, weight=wt, iter=0,
	      init=c(fit1$frail, coef(fit1)$fixed))
dt1 <- coxph.detail(tfit)

aeq(apply(dt1$score,2,sum)- c(fit1$frail, 0,0)/theta, fit1$u)
h1 <- apply(dt1$imat,1:2,sum) + diag(c(1/theta, 1/theta,0,0))
aeq(as.matrix(gchol(h1)), as.matrix(fit1$hmat))
aeq(solve(fit1$var, full=F), h1)

# And iteration 2
fit2 <- coxme(Surv(time1, time2, status) ~ x1 + x2, data=tdata0b,
              random= ~1|grp, variance=theta, weight=wt, iter=2)
aeq(fit1$u %*% fit1$var, c(fit2$frail, coef(fit2)$fixed) - 
                         c(fit1$frail, coef(fit1)$fixed))

rm(fit0, fit1, dt0, dt1, tfit, h0, h1, theta)
rm(fit2)
