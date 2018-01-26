# 
# This small test file is found in the "Survival Kit" software
#  of Ducrocq and Solkner
# Actually not useful -- by the end of the example I saw that they don't
#  print results from the random-effects model.
#
ducrocq <- read.table('ducrocq.data', row.names=NULL,
                    col.names=c('id', 'start', 'stop', 'status', 'sex', 'rx'))

dfit1 <- coxph(Surv(start, stop, status) ~ sex, ducrocq,
               method='breslow')
dfit2 <- coxph(Surv(start, stop, status) ~ sex + rx, ducrocq,
               method='breslow')

2* c(dfit1$loglik, dfit2$loglik[2])  #This matches printout on page 60
