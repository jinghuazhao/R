#
# Very simple data set, to test out crossed factors
#
simple <- data.frame(time=c(9:3,1,1), status=c(rep(0,7),1,1),
                     f1=rep(1:3,3), f2=c(rep(1,6), rep(2,3)), x=1:9)
sfit <- coxme(Surv(time, status) ~ x, data=simple, ties='breslow',
              random= ~1|f1/f2, variance=c(1,2), iter=0)

ta <- 2/9
tb <- -1/9
tc <- 4/27
td <- 14/81 
te <- 8/81 
tf <- -2/27
tg <- 2/27
th <- -1/27
ti <- -2/81
tj <- -4/81
tk <- -1/81
tx <- c(-3,0,3,-5,2,-3,3,-1,4,60)/9
itrue <- bdsmatrix(c(ta, tb, tb, tc, tg, tf, th, tf, th,
                         ta, tb, tf, th, tc, tg, tf, th,
                             ta, tf, th, tf, th, tc, tg,
                                 td, ti, tj, ti, tj, ti,
                                     te, ti, tk, ti, tk,
                                         td, ti, tj, ti,
                                             te, ti, tk,
                                                 td, ti, te),
                   blocksize=9, rmat=as.matrix(tx))
ibreslow <- 2*itrue
diag(ibreslow) <- diag(ibreslow) + rep(c(1,1/2,0), c(3,6,1))

aeq(as.matrix(igchol(sfit$hmat)), as.matrix(ibreslow))


sfit2 <- coxme(Surv(time, status) ~ x, data=simple, ties='breslow', 
              random= ~1|f1/f2, variance=c(1,2), iter=0, sparse.calc=1)
aeq(as.matrix(igchol(sfit2$hmat)), as.matrix(ibreslow))


# Now for the Efron approx
sfit <- coxme(Surv(time, status) ~ x, data=simple, ties='efron',
              random= ~1|f1/f2, variance=c(1,2), iter=0)

# the matrix for the first death, where each of the last 2 obs has weight
#  1/2
tx <- c(-54, -13, 67, -132, 78, -68, 55, -4, 71, 1471)
i2 <- bdsmatrix(c( 60, -30, -30,  40,  20, -24,  -6, -24,  -6, 
                        55, -25, -20, -10,  44,  11, -20,  -5, 
                             55, -20, -10, -20,  -5,  44,  11, 
                                  48,  -8, -16,  -4, -16,  -4, 
                                       28,  -8,  -2,  -8,  -2, 
                                            48,  -4, -16,  -4, 
                                                 15,  -4,  -1, 
                                                      48,  -4,     
		                                           15),
                   blocksize=9, rmat=as.matrix(tx))
iefron <- itrue + i2/256
diag(iefron) <- diag(iefron) +  rep(c(1,1/2,0), c(3,6,1)) #add penalty
aeq(as.matrix(igchol(sfit$hmat)), as.matrix(iefron))

sfit2 <- coxme(Surv(time, status) ~ x, data=simple, ties='efron', 
              random= ~1|f1/f2, variance=c(1,2), iter=0, sparse.calc=1)
aeq(as.matrix(igchol(sfit2$hmat)), as.matrix(iefron))

rm(ta, tb, tc, td,te, tf, tg,th, ti,tk,tx)
rm(sfit, itrue)
rm(sfit2, i2, iefron)
rm(simple)
