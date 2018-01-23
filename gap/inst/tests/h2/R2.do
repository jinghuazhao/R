/* 7-1-2015 MRC-Epid JHZ */

use gcta
// Cox model
set seed 1234567
foreach v in sex age bmi pa_index mdsbcn qge1302sw {
   di "`v'"
   str2ph stcox `v' lat_centre long_centre, bootreps(1000)
   str2d stcox `v' lat_centre long_centre
}
str2ph stcox sex age bmi pa_index mdsbcn qge1302sw, bootreps(1000)
str2d stcox sex age bmi pa_index mdsbcn qge1302sw

set seed 1234567
// cloglog
foreach v in sex age bmi pa_index mdsbcn qge1302sw {
  bootstrap (1-exp(-e(chi2)/e(N))) (1-exp(-e(chi2)/e(N_s))), reps(1000): cloglog diabetes_status `v'
}
bootstrap (1-exp(-e(chi2)/e(N))) (1-exp(-e(chi2)/e(N_s))), reps(1000):/*
*/ cloglog diabetes_status sex age bmi pa_index mdsbcn qge1302sw

set seed 1234567
// probit
foreach v in sex age bmi pa_index mdsbcn qge1302sw {
  quietly cloglog diabetes_status `v'
  local e=e(N_s)
  bootstrap (1-exp(-e(chi2)/e(N))) (1-exp(-e(chi2)/`e')), reps(1000): probit diabetes_status `v'
}
quietly cloglog diabetes_status sex age bmi pa_index mdsbcn qge1302sw
local e=e(N_s)
bootstrap (1-exp(-e(chi2)/e(N))) (1-exp(-e(chi2)/`e')), reps(1000):/*
*/ probit diabetes_status sex age bmi pa_index mdsbcn qge1302sw

set seed 1234567
// logit
foreach v in sex age bmi pa_index mdsbcn qge1302sw {
  quietly cloglog diabetes_status `v'
  local e=e(N_s)
  bootstrap (1-exp(-e(chi2)/e(N))) (1-exp(-e(chi2)/`e')), reps(1000): logit diabetes_status `v'
}
quietly cloglog diabetes_status sex age bmi pa_index mdsbcn qge1302sw
local e=e(N_s)
bootstrap (1-exp(-e(chi2)/e(N))) (1-exp(-e(chi2)/`e')), reps(1000):/*
*/ logit diabetes_status sex age bmi pa_index mdsbcn qge1302sw
