#
# trivial test of strata
#
tfit <- coxme(Surv(time, status) ~ x1 + x2, data=tdata0,
              random= ~1|grp, variance=.5, weight=wt, iter=0,
              ties='breslow', init=c(pi,0), sparse.calc=0)

sdata0 <- rbind(tdata0, tdata0)
sdata0$strat <- rep(1:2, each=nrow(tdata0))

sfit0 <- coxme(Surv(time,status) ~ x1 + x2 + strata(strat), data=sdata0,
	       random= ~1|grp, ties='breslow', iter=0, sparse.calc=0,
	       variance=0.5, weight=wt, init=c(pi,0))

sfit1 <- coxme(Surv(time,status) ~ x1 + x2 + strata(strat), data=sdata0,
	       random= ~1|grp, ties='breslow', iter=0, sparse.calc=1,
	       variance=0.5, weight=wt, init=c(pi,0))

aeq(sfit0$u, sfit1$u)
aeq(sfit0$hmat, sfit1$hmat)
aeq(sfit0$u/2, tfit$u)
tpen <- diag(c(2,2,0,0))
aeq(igchol(sfit0$hmat)-tpen, 2*(igchol(tfit$hmat)-tpen))

#
# Repeat for start/stop data
#   The subset makes the test more rigorous (someone exits a risk
#   set without entering it immediatly again).
tfit <- coxme(Surv(time1, time2, status) ~ x1 + x2, data=tdata0b,
              random= ~1|grp, variance=0.6, weight=wt, iter=0,
              ties='efron', init=c(pi,0), subset=(-1))

sdata0b <- rbind(tdata0b[-1,], tdata0b[-1,])
sdata0b$strat <- rep(1:2, each=nrow(tdata0b)-1)

sfit0 <- coxme(Surv(time1,time2,status) ~ x1 + x2 + strata(strat),
	       data=sdata0b, random=~1|grp, iter=0, sparse.calc=0,
	       variance=0.6, weight=wt, init=c(pi,0))
sfit1 <- coxme(Surv(time1,time2,status) ~ x1 + x2 + strata(strat),
	       data=sdata0b, random=~1|grp, iter=0, sparse.calc=1,
	       variance=0.6, weight=wt, init=c(pi,0))
aeq(tfit$u, sfit0$u/2)
aeq(sfit0$u, sfit1$u)
tpen <- diag(c(1/.6, 1/.6, 0,0))
aeq(sfit0$imat, sfit1$imat)
aeq(igchol(sfit0$hmat)-tpen, 2*(igchol(tfit$hmat)-tpen))
