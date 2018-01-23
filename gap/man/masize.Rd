\name{masize}
\alias{masize}
\title{Sample size calculation for mediation analysis}
\usage{masize(model,opts,alpha=0.025,gamma=0.2)}

\description{
The function computes sample size for regression problems where the goal is to assess mediation of the effects of a 
primary predictor by an intermediate variable or mediator.

Mediation has been thought of in terms of the proportion of effect explained, or the relative attenuation of b1, the 
coefficient for the primary predictor X1, when the mediator, X2, is added to the model. The goal is to show that b1*, 
the coefficient for X1 in the reduced model (i.e., the model with only X1, differs from b1, its coefficient in the full 
model (i.e., the model with both X1 and the mediator X2. If X1 and X2 are correlated, then showing that b2, the 
coefficient for X2, differs from zero is equivalent to showing b1* differs from b1. Thus the problem reduces to 
detecting an effect of X2, controlling for X1. In short, it amounts to the more familiar problem of inflating sample 
size to account for loss of precision due to adjustment for X1.

The approach here is to approximate the expected information matrix from the regression model including both X1 and X2, 
to obtain the expected standard error of the estimate of b2, evaluated at the MLE. The sample size follows from 
comparing the Wald test statistic (i.e., the ratio of the estimate of b2 to its SE) to the standard normal distribution, 
with the expected value of the numerator and denominator of the statistic computed under the alternative hypothesis. 
This reflects the Wald test for the statistical significance of a coefficient implemented in most regression packages.

The function provides methods to calculate sample sizes for the mediation problem for linear, logistic, Poisson, and Cox 
regression models in four cases for each model:

\tabular{ll}{
    CpCm \tab continuous primary predictor, continuous mediator\cr
    BpCm \tab binary primary predictor, continuous mediator\cr
    CpBm \tab continuous primary predictor, binary mediator\cr
    BpBm \tab binary primary predictor, binary mediator
}

The function is also generally applicable to the analogous problem of calculating sample size adequate to detect the 
effect of a primary predictor in the presence of confounding. Simply treat X2 as the primary predictor and consider X1 
the confounder.
}

\arguments{
 \item{model}{"lineari", "logisticj", "poissonk", "coxl", where i,j,k,l range from 1 to 4,5,9,9, respectively.}
 \item{opts}{A list specific to the model

   \tabular{ll}{
    b1 \tab regression coefficient for the primary predictor X1\cr
    b2 \tab regression coefficient for the mediator X2\cr
    rho \tab correlation between X1 and X2\cr
    sdx1, sdx2 \tab standard deviations (SDs) of X1 and X2\cr
    f1, f2 \tab prevalence of binary X1 and X2\cr
    sdy \tab residual SD of the outcome for the linear model\cr
    p \tab marginal prevalence of the binary outcome in the logistic model\cr
    m \tab marginal mean of the count outcome in a Poisson model\cr
    f \tab proportion of uncensored observations for the Cox model\cr
    fc \tab proportion of observations censored early\cr
    alpha \tab one-sided type-I error rate\cr
    gamma \tab type-II error rate\cr
    ns \tab number of observations to be simulated\cr
    seed \tab random number seed
    }

For linear model, the arguments are b2, rho, sdx2, sdy, alpha, and gamma. For cases CpBm and BpBm, set sdx2 = 
\eqn{\sqrt{f2(1-f2)}}{sqrt(f2*(1-f2))}. Three alternative functions are included for the linear model. These functions 
make it possible to supply other combinations of input parameters affecting mediation:

\tabular{ll}{
    b1* \tab coefficient for the primary predictor \cr
        \tab in the reduced model excluding the mediator (b1star)\cr
    b1  \tab coefficient for the primary predictor \cr
        \tab in the full model including the mediator \cr
    PTE \tab proportion of the effect of the primary predictor \cr
        \tab explained by the mediator, defined as (b1*-b1)/b1*\cr
}

These alternative functions for the linear model require specification of an extra parameter, but are provided for
convenience, along with two utility files for computing PTE and b1* from the other parameters. The required arguments
are explained in comments within the R code.

  }
 \item{alpha}{Type-I error rate, one-sided}
 \item{gamma}{Type-II error rate}
}

\details{
For linear model, a single function, linear, implements the analytic solution for all four cases, based on Hsieh et al., 
is to inflate sample size by a variance inflation factor, \eqn{1/(1-rho^2)}, where rho is the correlation of X1 and X2. 
This also turns out to be the analytic solution in cases CpCm and BpCm for the Poisson model, and underlies approximate 
solutions for the logistic and Cox models. An analytic solution is also given for cases CpBm and BpBm for the Poisson 
model. Since analytic solutions are not available for the logistic and Cox models, a simulation approach is used to 
obtain the expected information matrix instead.

For logistic model, the approximate solution due to Hsieh is implemented in the function logistic.approx, and can be 
used for all four cases. Arguments are p, b2, rho, sdx2, alpha, and gamma. For a binary mediator with prevalence f2, 
sdx2 should be reset to \eqn{\sqrt{f2(1-f2)}}{sqrt(f2*(1-f2))}. Simulating the information matrix of the logistic model 
provides somewhat more accurate sample size estimates than the Hsieh approximation. The functions for cases CpCm, BpCm, 
CpBm, and BpBm are respectively logistic.ccs, logistic.bcs, logistic.cbs, and logistic.bbs, as for the Poisson and Cox 
models. Arguments for these functions include p, b1, sdx1 or f1, b2, sdx2 or f2, rho, alpha, gamma, and ns. As in other 
functions, sdx1, sdx2, alpha, and gamma are set to the defaults listed above. These four functions call two utility 
functions, getb0 (to calculate the intercept parameter from the others) and antilogit, which are supplied.

For Poisson model, The function implementing the approximate solution based on the variance inflation factor is 
poisson.approx, and can be used for all four cases. Arguments are EY (the marginal mean of the Poisson outcome), b2, 
sdx2, rho, alpha and gamma, with sdx2, alpha and gamma set to the usual defaults; use 
sdx2=\eqn{\sqrt{f2(1-f2)}}{sqrt(f2*(1-f2))} for a binary mediator with prevalence f2 (cases CpBm and BpBm). For cases 
CpCm and BpCm (continuous mediators), the approximate formula is also the analytic solution. For these cases, we supply 
redundant functions poisson.cc and poisson.bc, with the same arguments and defaults as for poisson.approx (it's the same 
function). For the two cases with binary mediators, the functions are poisson.cb and poisson.bb. In addition to m, b2, 
f2, rho, alpha, and gamma, b1 and sdx1 or f1 must be specified. Defaults are as usual. Functions using simulation for 
the Poisson model are available: poisson.ccs, poisson.bcs, poisson.cbs, and poisson.bbs. As in the logistic case, these 
require arguments b1 and sdx1 or f1. For this case, however, the analytic functions are faster, avoid simulation error, 
and should be used. We include these functions as templates that could be adapted to other joint predictor 
distributions.

For Cox model, the function implementing the approximate solution, using the variance inflation factor and derived by 
Schmoor et al., is cox.approx, and can be used for all four cases. Arguments are b2, sdx2, rho, alpha, gamma, and f. For 
binary X2 set sdx2 = \eqn{\sqrt{f2(1-f2)}}{sqrt(f2*(1-f2))}. The approximation works very well for cases CpCm and BpCm 
(continuous mediators), but is a bit less accurate for cases CpBm and BpBm (binary mediators). We get some improvement 
for those cases using the simulation approach. This approach is implemented for all four, as functions cox.ccs, cox.bcs, 
cox.cbs, and cox.bbs. Arguments are b1, sdx1 or f1, b2, sdx2 or f2, rho, alpha, gamma, f, and ns, with defaults as 
described above. Slight variants of these functions, cox.ccs2, cox.bcs2, cox.cbs2, and cox.bbs2, make it possible to 
allow for early censoring of a fraction fc of observations; but in our experience this has virtually no effect, even 
with values of fc of 0.5. The default for fc is 0.

A summary of the argumentss is as follows, noting that additional parameter seed can be supplied for simulation-based 
method.

\tabular{lll}{
model \tab arguments \tab description\cr
\cr
linear1 \tab b2, rho, sdx2, sdy \tab linear \cr
linear2 \tab b1star, PTE, rho, sdx1, sdy \tab lineara\cr
linear3 \tab b1star, b2, PTE, sdx1, sdx2, sdy \tab linearb\cr
linear4 \tab b1star, b1, b2, sdx1, sdx2, sdy \tab linearc\cr
\cr
logistic1 \tab p, b2, rho, sdx2 \tab logistic.approx\cr
logistic2 \tab p, b1, b2, rho, sdx1, sdx2, ns \tab logistic.ccs\cr
logistic3 \tab p, b1, f1, b2, rho, sdx2, ns \tab logistic.bcs\cr
logistic4 \tab p, b1, b2, f2, rho, sdx1, ns \tab logistic.cbs\cr
logistic5 \tab p, b1, f1, b2, f2, rho, ns \tab logistic.bbs\cr
\cr
poisson1 \tab m, b2, rho, sdx2 \tab poisson.approx\cr
poisson2 \tab m, b2, rho, sdx2 \tab poisson.cc\cr
poisson3 \tab m, b2, rho, sdx2 \tab poisson.bc\cr
poisson4 \tab m, b1, b2, f2, rho, sdx1 \tab poisson.cb\cr
poisson5 \tab m, b1, f1, b2, f2, rho \tab poisson.bb\cr
poisson6 \tab m, b1, b2, rho, sdx1, sdx2, ns \tab poisson.ccs\cr
poisson7 \tab m, b1, f1, b2, rho, sdx2, ns \tab poisson.bcs\cr
poisson8 \tab m, b1, b2, f2, rho, sdx1, ns \tab poisson.cbs\cr
poisson9 \tab m, b1, f1, b2, f2, rho, ns\tab poisson.bbs\cr
\cr
cox1 \tab b2, rho, f, sdx2 \tab cox.approx\cr
cox2 \tab b1, b2, rho, f, sdx1, sdx2, ns \tab cox.ccs\cr
cox3 \tab b1, f1, b2, rho, f, sdx2, ns\tab cox.bcs\cr
cox4 \tab b1, b2, f2, rho, f, sdx1, ns\tab cox.cbs\cr
cox5 \tab b1, f1, b2, f2, rho, f, ns\tab cox.bbs\cr
cox6 \tab b1, b2, rho, f, fc, sdx1, sdx2, ns\tab cox.ccs2\cr
cox7 \tab b1, f1, b2, rho, f, fc, sdx2, ns\tab cox.bcs2\cr
cox8 \tab b1, b2, f2, rho, f, fc, sdx1, ns\tab cox.cbs2\cr
cox9 \tab b1, f1, b2, f2, rho, f, fc, ns\tab cox.bbs2\cr
}
}

\value{
A short description of model (desc, b=binary, c=continuous, s=simulation) and sample size (n). In the case of Cox model, 
number of events (d) is also indicated.
}

\section{References}{
Hsieh FY, Bloch DA, Larsen MD. A simple method of sample size calculation for linear and logistic regression. Stat 
Med 1998; 17:1623-34.

Schmoor C, Sauerbrer W, Schumacher M. Sample size considerations for the evaluation of prognostic factors in survival 
analysis. Stat Med 2000; 19:441-52.

Vittinghoff E, Sen S, McCulloch CE. Sample size calculations for evaluating mediation. Stat Med 2009; 
28:541-57.
}

\seealso{
\code{\link[gap]{ab}}
}

\source{\url{http://www.epibiostat.ucsf.edu/biostat/mediation/}}
\examples{
\dontrun{
## linear model
# CpCm
opts <- list(b2=0.5, rho=0.3, sdx2=1, sdy=1)
masize("linear1",opts)
# BpBm
opts <- list(b2=0.75, rho=0.3, f2=0.25, sdx2=sqrt(0.25*0.75), sdy=3)
masize("linear1",opts,gamma=0.1)

## logistic model
# CpBm
opts <- list(p=0.25, b2=log(0.5), rho=0.5, sdx2=0.5)
masize("logistic1",opts)
opts <- list(p=0.25, b1=log(1.5), sdx1=1, b2=log(0.5), f2=0.5, rho=0.5, ns=10000,
             seed=1234)
masize("logistic4",opts)
opts <- list(p=0.25, b1=log(1.5), sdx1=1, b2=log(0.5), f2=0.5, rho=0.5, ns=10000,
             seed=1234)
masize("logistic4",opts)
opts <- list(p=0.25, b1=log(1.5), sdx1=4.5, b2=log(0.5), f2=0.5, rho=0.5, ns=50000,
             seed=1234)
masize("logistic4",opts)

## Poisson model
# BpBm
opts <- list(m=0.5, b2=log(1.25), rho=0.3, sdx2=sqrt(0.25*0.75))
masize("poisson1",opts)
opts <- list(m=0.5, b1=log(1.4), f1=0.25, b2=log(1.25), f2=0.25, rho=0.3)
masize("poisson5",opts)
opts <- c(opts,ns=10000, seed=1234)
masize("poisson9",opts)

## Cox model
# BpBm
opts <- list(b2=log(1.5), rho=0.45, f=0.2, sdx2=sqrt(0.25*0.75))
masize("cox1",opts)
opts <- list(b1=log(2), f1=0.5, b2=log(1.5), f2=0.25, rho=0.45, f=0.2, seed=1234)
masize("cox5",c(opts, ns=10000))
masize("cox5",c(opts, ns=50000))
}
}
\keyword{misc}
