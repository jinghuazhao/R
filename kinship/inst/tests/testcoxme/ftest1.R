#
# Try a simple frailty fit to a few families
#   Choose some families with many events
# This fit has no covariates, to make missing value issues go away. ALso, 
#   put the data set in id order.  This makes various coefficients line
#   up in the outputs.  We also re-label the id's for ease of reading.
# The next file, ftest2, does more families and a covariate
#
# data.restore('cdata.sdump')
cdata <- dget("cdata.dput")
keep <- (cdata$famid==8)  #use this one for intitial testing
tdata <- cdata[keep,c("gid", "momid", "dadid", "parity", "famid",
			      "startage", "endage", "cancer", 'sex')]

# Make the kinship matrix of the females
#  for comparison to other frailty models, we want 2*kinship (a correlation
#  matrix)
# The frailty.kin function doesn't deal quite properly with data deleted
#  due to missing -- the final vector of frailty coefs will have one
#  element for each obs present BEFORE missings were excluded. Thus "who2".
#  But it's not worth fixing up what is just a function for testing/validation
tkmat <- kinship(tdata$gid, tdata$dadid, tdata$momid)
who  <- (tdata$sex=='F')
who2 <- (who & !is.na(tdata$startage + tdata$endage + tdata$cancer))
tkmat <- 2* tkmat[who2,who2]

#
# Because the data set is small enough, we can use the frailty.kin
#   function.  It uses the full kmat (with zeros), without the block
#   diagonal computational speedups.
# Get the 0 iteration, 1 iteration, and full solution
tkfit0 <- coxph(Surv(startage, endage, cancer) ~
	             frailty.kin(gid, theta=.8, kmat=tkmat), tdata[who2,], 
		iter=0, x=T)
tkfit1 <- coxph(Surv(startage, endage, cancer) ~
	      frailty.kin(gid, theta=.8, kmat=tkmat), tdata[who2,], iter=1)
tkfit <- coxph(Surv(startage, endage, cancer) ~ 
	      frailty.kin(gid, theta=.8, kmat=tkmat), tdata[who2,])


# Second fit, use the coxme function, which uses the block-diagonal
#   version of the kinship, and an "external" iteration (in its first
#   incarnation)
# The coxme function can deal with deleted rows, so I use "who" (females
#   only, not "who2" (females + no missing values)
#
tkmat2 <- makekinship(tdata$famid, tdata$gid, tdata$momid, tdata$dadid)
aeq(tkmat, 2*as.matrix(tkmat2[who2, who2]))  # just a check

kfit0  <- coxme(Surv(startage, endage, cancer) ~ 1, data=tdata[who,],
                random= ~1|gid, varlist=tkmat2, variance=.8, iter=0)
kfit1  <- coxme(Surv(startage, endage, cancer) ~ 1, data=tdata[who,],
                random= ~1|gid, varlist=tkmat2, variance=.8, iter=1)
kfit   <- coxme(Surv(startage, endage, cancer) ~ 1, data=tdata[who,],
                random= ~1|gid, varlist=tkmat2, variance=.8)


# compare -- they won't be exactly the same due to different iteration
#  paths; the coxme call used some off diagonal elements of imat that
#  tkfit did not 
all(abs(range(tkfit$coef - kfit$frail)) < 1e-4)

#
# Compute the desired update steps in full detail
#   For family 8 there are 19 females.
#      #1 (proband) and #3 have cancer=NA, so get dropped
#      there are 2 events at 28.5 and 38.5, rows 6 and 1 of the 17, respective
#      row 2 of the 17 overlaps no events, O = E = U = zero
#
nf  <- length(tkfit$coef)
dummy <- coxph(Surv(startage, endage, cancer) ~ diag(nf), tdata[who2,], iter=0)
dtdum <- coxph.detail(dummy)
imat.full <- apply(dtdum$imat, 1:2, sum)
u.full    <- apply(dtdum$score, 2,  sum)
pmat <- solve(tkmat)/ .8
update.full<- coxph.wtest(imat.full + pmat, u.full)$solve

# compare to the update step the full matrix gave
aeq(update.full, tkfit1$coef)
aeq(update.full, kfit1$frail)

#
# Now, check out the second iteration too
#
kfit2  <- coxme(Surv(startage, endage, cancer) ~ 1, data=tdata[who,],
                random= ~1|gid, varlist=tkmat2, variance=.8, iter=2)

upen <- c(kfit1$frail %*% pmat)
dummy2 <- coxph(Surv(startage, endage, cancer) ~ diag(nf), tdata[who2,], 
		    iter=0, init=kfit1$frail)
dt2 <- coxph.detail(dummy2)
imat2 <- apply(dt2$imat, 1:2, sum)
u2    <- apply(dt2$score, 2, sum)
update2 <- coxph.wtest(imat2 + pmat, u2-upen)$solve 
aeq(u2-upen, kfit1$u)
aeq(update2, kfit2$frail - kfit1$frail)

rm (u2, upen, update2, dummy2, dt2, imat2, 
    nf, dummy, dtdum, u.full, update.full, imat.full,
    pmat, tkmat, tkmat2, who, who2)
rm(keep, kfit0, kfit1, kfit2, kfit)
rm(tkfit1, tkfit0, tkfit)
rm(tdata)
