#
# This test does not work out right -- close, but not quite right
#  The issue is REML(coxph) vs ML (coxme) estimates.
#
rfit1 <- coxme(Surv(time, status) ~ rx, random= ~1|litter, data=rats,
              variance=.3)
rfit2 <- coxph(Surv(time, status) ~ rx + 
	         frailty(litter, dist='gauss', theta=.3), rats, eps=1e-10)
rfit3 <- coxph(Surv(time, status) ~ rx + 
	         frailty(litter, dist='gauss', theta=.3, sparse=F), rats)


