#' Sample Size Calculation for Mediation Analysis
#'
#' Computes sample size requirements for regression models in which the effect
#' of a primary predictor (\eqn{X_1}) is evaluated in the presence of a mediator
#' (\eqn{X_2}). The calculation is based on the expected information matrix of
#' the full model and a Wald test under the alternative hypothesis.
#'
#' The framework applies to linear, logistic, Poisson, and Cox regression models
#' and supports four predictor–mediator configurations:
#'
#' - CpCm: continuous X1 and continuous X2.
#' - BpCm: binary X1 and continuous X2.
#' - CpBm: continuous X1 and binary X2.
#' - BpBm: binary X1 and binary X2.
#'
#' The naming convention is:
#' - "cc" continuous-continuous.
#' - "bc" binary-continuous.
#' - "cb" continuous-binary.
#' - "bb" binary-binary.
#' - ".approx" use analytic approximations.
#' - ".s" use simulation.
#'
#' The same formulation can also be used to evaluate the effect of a predictor
#' in the presence of a confounder by swapping the roles of X1 and X2.
#'
#' Throughout, X1 denotes the primary predictor and X2 denotes the mediator.
#' Coefficients b1 and b2 correspond to X1 and X2 respectively, and f1 and f2
#' denote prevalences when variables are binary. rho denotes the correlation
#' between X1 and X2.
#'
#' The function is also applicable to the analogous problem of detecting the effect
#' of a primary predictor in the presence of a confounder by treating X2 as the
#' primary predictor and X1 as the confounder.
#'
#' @param model A character string specifying the model type:
#'  - Linear models: "linear1", "linear2", "linear3", "linear4"
#'  - Logistic models: "logistic1"–"logistic5"
#'  - Poisson models: "poisson1"–"poisson9"
#'  - Cox models: "cox1"–"cox9"
#'
#' @param opts Named list of model-specific parameters, which may include:
#' - b1 regression coefficient for X1.
#' - b1star coefficient for X1 in the reduced model.
#' - b2 regression coefficient for X2.
#' - PTE proportion of total effect explained.
#' - rho correlation between X1 and X2.
#' - sdx1, sdx2 standard deviations of X1 and X2.
#' - f1, f2 prevalences for binary X1 and X2.
#' - sdy residual SD (linear models).
#' - p marginal outcome prevalence (logistic models).
#' - m marginal outcome mean (Poisson models).
#' - f proportion of uncensored observations (Cox models).
#' - fc early censoring proportion (Cox extensions).
#' - ns simulation size for simulation-based methods.
#' - seed random seed.
#'
#' @param alpha One-sided type I error rate.
#' @param gamma Type-II error rate.
#'
#' @details
#' Let b1* (b1star) denote the coefficient of X1 in the reduced model excluding X2, and b1 the
#' coefficient in the full model. Mediation can be expressed as attenuation b1* - b1 or proportion
#' of total effect explained (PTE):
#'
#' \deqn{PTE = \frac{b_1^*-b_1}{b_1^*}}
#'
#' Testing mediation is equivalent (linear models) or approximately equivalent (GLMs and Cox models)
#' to testing the coefficient of X2 in the full model. Sample size is therefore derived from a Wald
#' test on b2 using the expected information matrix.
#'
#' For binary variables, use \eqn{\sqrt{f(1-f)}} for the corresponding standard deviation.
#'
#' Linear models parameters are based on \insertCite{hsieh98;textual}{gap} and use a variance inflation
#' factor \eqn{1/(1-\rho^2)}, where \eqn{\rho} is the correlation between X1 and X2. Logistic and Cox models
#' generally require simulation to estimate information, while Poisson models admit analytic solutions in
#' most cases.
#'
#' The quantity PTE and b1star parameter for linear model are provided as follows,
#'
#' \preformatted{
#'  getPTE <- function(b1, b2, rho, sdx1=1, sdx2=1)
#'                     b2*sdx2*rho/((b1+b2*sdx2*rho/sdx1)*sdx1)
#'  getb1star <- function(b1, b2, rho, sdx1=1, sdx2=1) b1+b2*sdx2*rho/sdx1
#'  }
#'
#' This also turns out to be the analytic solution in cases CpCm and BpCm for the Poisson model, and underlies approximate
#' solutions for the logistic and Cox models. An analytic solution is also given for cases CpBm and BpBm for the Poisson
#' model. Since analytic solutions are not available for the logistic and Cox models, a simulation approach is used to
#' obtain the expected information matrix instead.
#'
#' For logistic model, the approximate solution due to Hsieh is implemented in the function logistic.approx, and can be
#' used for all four cases. Arguments are p, b2, rho, sdx2, alpha, and gamma. Simulating the information matrix of the logistic model
#' provides somewhat more accurate sample size estimates than the Hsieh approximation. The functions for cases CpCm, BpCm,
#' CpBm, and BpBm are respectively logistic.ccs, logistic.bcs, logistic.cbs, and logistic.bbs, as for the Poisson and Cox
#' models. Arguments for these functions include p, b1, sdx1 or f1, b2, sdx2 or f2, rho, alpha, gamma, and ns. As in other
#' functions, sdx1, sdx2, alpha, and gamma are set to the defaults listed above. These four functions call the utility
#' functions and getb0 (to calculate the intercept parameter from the others), which are supplied.
#'
#' \preformatted{
#' getb0 <- function(p, X, b1, b2)
#' {
#'   stopifnot(is.numeric(p), length(p) == 1, p > 0, p < 1)
#'   eta0 <- X[,2] * b1 + X[,3] * b2
#'   f <- function(b0) mean(plogis(b0+eta0)) - p
#'   uniroot(f, interval = c(-50,50), tol = 1e-10)$root
#' }
#' }
#'
#' For Poisson model, The function implementing the approximate solution based on the variance inflation factor is
#' poisson.approx, and can be used for all four cases. Arguments are EY (the marginal mean of the Poisson outcome), b2,
#' sdx2, rho, alpha and gamma, with sdx2, alpha and gamma set to the usual defaults; For cases 
#' CpCm and BpCm (continuous mediators), the approximate formula is also the analytic solution. For these cases, we supply
#' redundant functions poisson.cc and poisson.bc, with the same arguments and defaults as for poisson.approx (it's the same
#' function). For the two cases with binary mediators, the functions are poisson.cb and poisson.bb. In addition to m, b2,
#' f2, rho, alpha, and gamma, b1 and sdx1 or f1 must be specified. Defaults are as usual. Functions using simulation for
#' the Poisson model are available: poisson.ccs, poisson.bcs, poisson.cbs, and poisson.bbs. As in the logistic case, these
#' require arguments b1 and sdx1 or f1. For this case, however, the analytic functions are faster, avoid simulation error,
#' and should be used. We include these functions as templates that could be adapted to other joint predictor distributions.
#'
#' For Cox model, the function implementing the approximate solution, using the variance inflation factor and derived by
#' \insertCite{schmoor00;textual}{gap} is cox.approx, and can be used for all four cases. Arguments are b2, sdx2, rho,
#' alpha, gamma, and f. The approximation works very well for cases CpCm and BpCm (continuous mediators), but is a bit
#' less accurate for cases CpBm and BpBm (binary mediators). #' We get some improvement for those cases using the simulation
#' approach. This approach is implemented for all four, as functions cox.ccs, cox.bcs, cox.cbs, and cox.bbs. Arguments are
#' b1, sdx1 or f1, b2, sdx2 or f2, rho, alpha, gamma, f, and ns, with defaults as described above. Slight variants of these
#' functions, cox.ccs2, cox.bcs2, cox.cbs2, and cox.bbs2, make it possible to allow for early censoring of a fraction fc
#' of observations; but in our experience this has virtually no effect, even with values of fc of 0.5. The default for fc is 0.
#'
#' A summary of the arguments is as follows, noting that additional parameter seed can be supplied for simulation-based
#' method.
#'
#' Model   | Arguments          | Description
#' --------|--------------------|--------------------------
#' Linear models: |
#' linear1 | b2, rho, sdx2, sdy  | linear
#' linear2 | b1star, PTE, rho, sdx1, sdy | lineara
#' linear3 | b1star, b2, PTE, sdx1, sdx2, sdy | linearb
#' linear4 | b1star, b1, b2, sdx1, sdx2, sdy | linearc
#' Logistic models: |
#' logistic1 | p, b2, rho, sdx2 | logistic.approx
#' logistic2 | p, b1, b2, rho, sdx1, sdx2, ns | logistic.ccs
#' logistic3 | p, b1, f1, b2, rho, sdx2, ns | logistic.bcs
#' logistic4 | p, b1, b2, f2, rho, sdx1, ns | logistic.cbs
#' logistic5 | p, b1, f1, b2, f2, rho, ns | logistic.bbs
#' Poission models: |
#' poisson1 | m, b2, rho, sdx2 | poisson.approx
#' poisson2 | m, b2, rho, sdx2 | poisson.cc
#' poisson3 | m, b2, rho, sdx2 | poisson.bc
#' poisson4 | m, b1, b2, f2, rho, sdx1 | poisson.cb
#' poisson5 | m, b1, f1, b2, f2, rho | poisson.bb
#' poisson6 | m, b1, b2, rho, sdx1, sdx2, ns | poisson.ccs
#' poisson7 | m, b1, f1, b2, rho, sdx2, ns | poisson.bcs
#' poisson8 | m, b1, b2, f2, rho, sdx1, ns | poisson.cbs
#' poisson9 | m, b1, f1, b2, f2, rho, ns | poisson.bbs
#' Cox models: |
#' cox1 | b2, rho, f, sdx2 | cox.approx
#' cox2 | b1, b2, rho, f, sdx1, sdx2, ns | cox.ccs
#' cox3 | b1, f1, b2, rho, f, sdx2, ns | cox.bcs
#' cox4 | b1, b2, f2, rho, f, sdx1, ns | cox.cbs
#' cox5 | b1, f1, b2, f2, rho, f, ns | cox.bbs
#' cox6 | b1, b2, rho, f, fc, sdx1, sdx2, ns | cox.ccs2
#' cox7 | b1, f1, b2, rho, f, fc, sdx2, ns | cox.bcs2
#' cox8 | b1, b2, f2, rho, f, fc, sdx1, ns | cox.cbs2
#' cox9 | b1, f1, b2, f2, rho, f, fc, ns | cox.bbs2
#'
#' @return A list containing:
#' - desc Description of the selected method.
#' - n Required sample size.
#' - d Required number of events (Cox models only).
#'
#' @export
#'
#' @references
#' \insertRef{hsieh98}{gap}
#'
#' \insertRef{schmoor00}{gap}
#'
#' \insertRef{vittinghoff09}{gap}
#'
#' @seealso [`ab`]
#'
#' @examples
#' \dontrun{
#' ## linear model
#' # CpCm
#' opts <- list(b2=0.5, rho=0.3, sdx2=1, sdy=1)
#' masize("linear1",opts)
#' # BpBm
#' opts <- list(b2=0.75, rho=0.3, f2=0.25, sdx2=sqrt(0.25*0.75), sdy=3)
#' masize("linear1",opts,gamma=0.1)
#'
#' ## logistic model
#' # CpBm
#' opts <- list(p=0.25, b2=log(0.5), rho=0.5, sdx2=0.5)
#' masize("logistic1",opts)
#' opts <- list(p=0.25, b1=log(1.5), sdx1=1, b2=log(0.5), f2=0.5, rho=0.5, ns=10000,
#'              seed=1234)
#' masize("logistic4",opts)
#' opts <- list(p=0.25, b1=log(1.5), sdx1=1, b2=log(0.5), f2=0.5, rho=0.5, ns=10000,
#'              seed=1234)
#' opts <- list(p=0.25, b1=log(1.5), sdx1=4.5, b2=log(0.5), f2=0.5, rho=0.5, ns=50000,
#'              seed=1234)
#' masize("logistic4",opts)
#'
#' ## Poisson model
#' # BpBm
#' opts <- list(m=0.5, b2=log(1.25), rho=0.3, sdx2=sqrt(0.25*0.75))
#' masize("poisson1",opts)
#' opts <- list(m=0.5, b1=log(1.4), f1=0.25, b2=log(1.25), f2=0.25, rho=0.3)
#' masize("poisson5",opts)
#' opts <- c(opts,ns=10000, seed=1234)
#' masize("poisson9",opts)
#'
#' ## Cox model
#' # BpBm
#' opts <- list(b2=log(1.5), rho=0.45, f=0.2, sdx2=sqrt(0.25*0.75))
#' masize("cox1",opts)
#' opts <- list(b1=log(2), f1=0.5, b2=log(1.5), f2=0.25, rho=0.45, f=0.2, seed=1234)
#' masize("cox5",c(opts, ns=10000))
#' masize("cox5",c(opts, ns=50000))
#' }
#'
#' @keywords misc

masize <- function(model,opts, alpha=0.025, gamma=0.2)
{
   for(p in c("survival")) {
      if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
         if (!requireNamespace(p, quietly = TRUE))
         warning(paste("masize needs package `", p, "' to be fully functional; please install", sep=""))
      }
   }
   linear <- paste("linear",1:4,sep="")
   logistic <- paste("logistic",1:5,sep="")
   poisson <- paste("poisson",1:9,sep="")
   cox <- paste("cox",1:9,sep="")
   model.int <- charmatch(model, c(linear,logistic,poisson,cox))
   if (is.na(model.int)) stop("Invalid model type")
   if (model.int == 0) stop("Ambiguous model type")
   if (!is.null(opts$seed)) set.seed(opts$seed)
   ## linear models
   if (model.int == 1) # (2) linear
   {
       b2 <- opts$b2       # regression coefficient for mediator in full model
       rho <- opts$rho     # correlation of primary predictor and mediator
       sdx2 <- opts$sdx2   # SD of mediator; use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
       sdy <- opts$sdy     # SD of outcome
       desc <- "linear"
       n <- (qnorm(alpha)+qnorm(gamma))^2*sdy^2/((b2*sdx2)^2*(1-rho^2))
   }
   if (model.int == 2) # (4) lineara
   {
        b1star <- opts$b1star # regression coefficient for primary predictor in reduced model
        PTE <- opts$PTE       # proportion of effect of primary predictor explained by mediator: (b1star - b1)/b1star
        rho <- opts$rho       # correlation of primary predictor and mediator
        sdx1 <- opts$sdx1     # SD of primary predictor; use sqrt(f1*(1-f1)) for binary predictor with prevalence f
        sdy <- opts$sdy       # SD of outcome
        desc <- "lineara"
        n <- (qnorm(alpha)+qnorm(gamma))^2*rho^2*sdy^2/((b1star*sdx1*PTE)^2*(1-rho^2))
   }
   if (model.int == 3) # (5) linearb
   {
        b1star <- opts$b1star # regression coefficient for primary predictor in reduced model
        b2 <- opts$b2         # regression coefficient for mediator in full model
        PTE <- opts$PTE       # proportion of effect of primary predictor explained by mediator: (b1star - b1)/b1star
        sdx1 <- opts$sdx1     # SD of primary predictor; use sqrt(f1*(1-f1)) for binary predictor with prevalence f1
        sdx2 <- opts$sdx2     # SD of mediator; use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
        sdy <- opts$sdy       # SD of outcome
        desc <- "linearb"
        n <- (qnorm(alpha)+qnorm(gamma))^2*sdy^2/((b2*sdx2)^2-(b1star*sdx1*PTE)^2)
   }
   if (model.int == 4) # (6) linearc
   {
        b1star <- opts$b1star # regression coefficient for primary predictor in reduced model
        b1 <- opts$b1         # regression coefficient for primary predictor in full model
        b2 <- opts$b2         # regression coefficient for mediator in full model
        sdx1 <- opts$sdx1     # SD of primary predictor; use sqrt(f1*(1-f1)) for binary predictor with prevalence f1
        sdx2 <- opts$sdx2     # SD of mediator; use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
        sdy <- opts$sdy       # SD of outcome
        desc <- "linearc"
        n <- (qnorm(alpha)+qnorm(gamma))^2*sdy^2/((b2*sdx2)^2-((b1star-b1)*sdx1)^2)
   }
   ## logistic model
   if (model.int == 5) # (7) logistic.approx
   {
        p <- opts$p           # marginal prevalence of outcome
        b2 <- opts$b2         # regression coefficient (log odds-ratio) for mediator
        rho <- opts$rho       # correlation of primary predictor and mediator
        sdx2 <- opts$sdx2     # SD of mediator; use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
        desc <- "logistic.approx"
        n <- (qnorm(alpha)+qnorm(gamma))^2/((b2*sdx2)^2*(1-rho^2)*p*(1-p))
   }
   if (model.int == 6) # (8) logistic.ccs
   {
        p <- opts$p           # marginal prevalence of outcome
        b1 <- opts$b1         # regression coefficient (log odds-ratio) for primary predictor
        b2 <- opts$b2         # regression coefficient (log odds-ratio) for mediator
        rho <- opts$rho       # correlation of primary predictor and mediator
        sdx1 <- opts$sdx1     # SD of primary predictor
        sdx2 <- opts$sdx2     # SD of mediator
        ns <- opts$ns         # number of observations in simulated dataset
        S <- matrix(c(sdx1^2, rho*sdx1*sdx2, rho*sdx1*sdx2, sdx2^2), 2, 2)
        X <- cbind(rep(1, ns), matrix(rnorm(2*ns), ns, 2) %*% chol(S))
        b0 <- getb0(p, X, b1, b2)
        pi <- plogis(X%*%as.matrix(c(b0, b1, b2)))
        XVX <- crossprod(X*as.vector(sqrt(pi*(1-pi))))
        desc <- "logistic.ccs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/(b2^2)
   }
   if (model.int == 7) # (8) logistic.bcs
   {
        p <- opts$p                      # marginal prevalence of outcome
        b1 <- opts$b1                    # regression coefficient (log odds-ratio) for primary predictor
        f1 <- opts$f1                    # prevalence of binary primary predictor
        b2 <- opts$b2                    # regression coefficient (log odds-ratio) for mediator
        rho <- opts$rho                  # correlation of primary predictor and mediator
        sdx2 <- opts$sdx2                # SD of mediator
        ns <- opts$ns                    # number of observations in simulated dataset
        mu0 <- -rho*sdx2*sqrt(f1/(1-f1)) # conditional mean of mediator when X1=0
        mu1 <- rho*sdx2*sqrt((1-f1)/f1)  # conditional mean of mediator when X1=1
        sdx2.1 <- sdx2*sqrt(1-rho^2)     # SD of mediator within levels of X1
        n0 <- round(ns*(1-f1))           # number of simulated observations with X1=0
        n1 <- ns-n0                      # number of simulated observations with X1=1
        X <- rbind(cbind(rep(1, n0), rep(0, n0), rnorm(n0, mean=mu0, sd=sdx2.1)),
                   cbind(rep(1, n1), rep(1, n1), rnorm(n1, mean=mu1, sd=sdx2.1)))
        b0 <- getb0(p, X, b1, b2)
        pi <- plogis(X%*%as.matrix(c(b0, b1, b2)))
        XVX <- crossprod(X*as.vector(sqrt(pi*(1-pi))))
        desc <- "logistic.bcs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/b2^2
   }
   if (model.int == 8) # (8) logistic.cbs
   {
        p <- opts$p                     # marginal prevalence of outcome
        b1 <- opts$b1                   # regression coefficient (log odds-ratio) for continuous primary predictor
        b2 <- opts$b2                   # regression coefficient (log odds-ratio) for binary mediator
        f2 <- opts$f2                   # prevalence of binary mediator
        rho <- opts$rho                 # correlation of primary predictor and mediator
        sdx1 <- opts$sdx1               # SD of continuous primary predictor
        ns <- opts$ns                   # number of observations in simulated dataset
        mu0 <- -rho*sdx1*sqrt(f2/(1-f2))
        mu1 <- rho*sdx1*sqrt((1-f2)/f2)
        sdx1.2 <- sdx1*sqrt(1-rho^2)
        n0 <- round(ns*(1-f2))
        n1 <- ns-n0
        X <- rbind(cbind(rep(1, n0), rnorm(n0, mean=mu0, sd=sdx1.2), rep(0, n0)),
                   cbind(rep(1, n1), rnorm(n1, mean=mu1, sd=sdx1.2), rep(1, n1)))
        b0 <- getb0(p, X, b1, b2)
        pi <- plogis(X%*%as.matrix(c(b0, b1, b2)))
        XVX <- crossprod(X*as.vector(sqrt(pi*(1-pi))))
        desc <- "logistic.cbs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/b2^2
   }
   if (model.int == 9) # (8) logistic.bbs
   {
        p <- opts$p                     # marginal prevalence of outcome
        b1 <- opts$b1                   # regression coefficient (log odds-ratio) for binary primary predictor
        f1 <- opts$f1                   # prevalence of primary predictor
        b2 <- opts$b2                   # regression coefficient (log odds-ratio) for binary mediator
        f2 <- opts$f2                   # prevalence of mediator
        rho <- opts$rho                 # correlation of primary predictor and mediator
        ns <- opts$ns                   # number of observations in simulated dataset
        f11 <- rho*sqrt(f1*(1-f1)*f2*(1-f2))+f1*f2 # fraction with X1=1 and X2=1
        f10 <- f1-f11                   # fraction with X1=1 and X2=0
        f01 <- f2-f11                   # fraction with X1=0 and X2=1
        n11 <- round(ns*f11)
        n10 <- round(ns*f10)
        n01 <- round(ns*f01)
        n00 <- max(1, ns-n10-n01-n11)
        X <- rbind(cbind(rep(1, n00), rep(0, n00), rep(0, n00)),
                   cbind(rep(1, n10), rep(1, n10), rep(0, n10)),
                   cbind(rep(1, n01), rep(0, n01), rep(1, n01)),
                   cbind(rep(1, n11), rep(1, n11), rep(1, n11)))
        b0 <- getb0(p, X, b1, b2)
        pi <- plogis(X%*%as.matrix(c(b0, b1, b2)))
        XVX <- crossprod(X*as.vector(sqrt(pi*(1-pi))))
        desc <- "logistic.bbs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/b2^2
   }
   ## Poisson model
   if (model.int == 10) # (9) poisson.approx
   {
        m <- opts$m                     # marginal mean of outcome
        b2 <- opts$b2                   # regression coefficient (log rate-ratio) for mediator
        rho <- opts$rho                 # correlation of primary predictor and mediator
        sdx2 <- opts$sdx2               # SD of mediator, use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
        desc <- "poisson.approx"
        n <- (qnorm(alpha)+qnorm(gamma))^2/((b2*sdx2)^2*(1-rho^2)*m)
   }
   if (model.int == 11) # (9) poisson.cc
   {
        m <- opts$m                     # marginal mean of outcome
        b2 <- opts$b2                   # regression coefficient (log rate-ratio) for mediator
        rho <- opts$rho                 # correlation of primary predictor and mediator
        sdx2 <- opts$sdx2               # SD of mediator, use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
        desc <- "poisson.cc"
        n <- (qnorm(alpha)+qnorm(gamma))^2/((b2*sdx2)^2*(1-rho^2)*m)
   }
   if (model.int == 12) # (9) poisson.bc
   {
        m <- opts$m                     # marginal mean of outcome
        b2 <- opts$b2                   # regression coefficient (log rate-ratio) for mediator
        rho <- opts$rho                 # correlation of primary predictor and mediator
        sdx2 <- opts$sdx2               # SD of mediator, use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
        desc <- "poisson.bc"
        n <- (qnorm(alpha)+qnorm(gamma))^2/((b2*sdx2)^2*(1-rho^2)*m)
   }
   if (model.int == 13) # (10) poisson.cb
   {
        m <- opts$m                      # marginal mean of outcome
        b1 <- opts$b1                    # regression coefficient (log rate-ratio) for continuous primary predictor
        b2 <- opts$b2                    # regression coefficient (log rate-ratio) for binary mediator
        f2 <- opts$f2                    # prevalence of binary mediator
        rho <- opts$rho                  # correlation of primary predictor and mediator
        sdx1 <- opts$sdx1                # SD of continuous primary predictor
        mu0 <- -rho*sdx1*sqrt(f2/(1-f2)) # mean of X1 when X2 = 0
        mu1 <- rho*sdx1*sqrt((1-f2)/f2)  # mean of X1 when X2 = 1
        sdx1.2 <- sdx1*sqrt(1-rho^2)
        B <- f2*exp(b1*mu1+b2)
        D <- (1-f2)*exp(b1*mu0)
        desc <- "poisson.cb"
        n <- (qnorm(alpha)+qnorm(gamma))^2*(((B+D)*sdx1.2)^2+B*D*(mu1-mu0)^2)/(b2^2*B*D*sdx1.2^2*m)
   }
   if (model.int == 14) # (11) poisson.bb
   {
        m <- opts$m                     # marginal mean of outcome
        b1 <- opts$b1                   # regression coefficient (log rate-ratio) for binary primary predictor
        f1 <- opts$f1                   # prevalence of primary predictor
        b2 <- opts$b2                   # regression coefficient (log rate-ratio) for binary mediator
        f2 <- opts$f2                   # prevalence of binary mediator
        rho <- opts$rho                 # correlation of primary predictor and mediator
        f11 <- f1*f2+rho*sqrt(f1*(1-f1)*f2*(1-f2))
        f10 <- f1-f11
        f01 <- f2-f11
        f00 <- 1-f01-f10-f11
        B <- f00
        C <- f10*exp(b1)
        D <- f01*exp(b2)
        E <- f11*exp(b1+b2)
        desc <- "poisson.bb"
        n <- (qnorm(alpha)+qnorm(gamma))^2*(B+C+D+E)*(B+D)*(C+E)/(b2^2*(B*C*D+B*C*E+B*D*E+C*D*E)*m)
   }
   if (model.int == 15) # (12) poisson.ccs
   {
        m <- opts$m     # marginal mean of outcome
        b1 <- opts$b1   # regression coefficient (log rate-ratio) for continuous primary predictor
        b2 <- opts$b2   # regression coefficient (log rate-ratio) for continuous mediator
        rho <- opts$rho # correlation of primary predictor and mediator
        sdx1 <- opts$sdx1 # SD of primary predictor
        sdx2 <- opts$sdx2 # SD of mediator
        ns <- opts$ns     # number of observations in simulated dataset
        S <- matrix(c(sdx1^2, rho*sdx1*sdx2, rho*sdx1*sdx2, sdx2^2), 2, 2)
        X <- matrix(rnorm(2*ns), ns, 2) %*% chol(S)
        rr <- exp(X%*%matrix(c(b1, b2)))
        XVX <- crossprod(cbind(rep(1, ns), X)*as.vector(sqrt(m*rr/mean(rr))))
        desc <- "poisson.ccs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/(b2^2)
   }
   if (model.int == 16) # (12) poisson.bcs
   {
        m <- opts$m                     # marginal mean of outcome
        b1 <- opts$b1                   # regression coefficient (log rate-ratio) for continuous primary predictor
        f1 <- opts$f1                   # prevalence of primary predictor
        b2 <- opts$b2                   # regression coefficient (log rate-ratio) for continuous mediator
        rho <- opts$rho                 # correlation of primary predictor and mediator
        sdx2 <- opts$sdx2               # SD of mediator
        ns <- opts$ns                   # number of observations in simulated dataset
        n1 <- round(ns*f1)
        n0 <- ns-n1
        mu0 <- -rho*sdx2*sqrt(f1/(1-f1)) # mean of X2 when X1 = 0
        mu1 <- rho*sdx2*sqrt((1-f1)/f1) # mean of X2 when X1 = 1
        sdx2.1 <- sdx2*sqrt(1-rho^2)
        X <- rbind(cbind(rep(0, n0), rnorm(n0, mean=mu0, sd=sdx2.1)),
                   cbind(rep(1, n1), rnorm(n1, mean=mu1, sd=sdx2.1)))
        rr <- exp(X%*%matrix(c(b1, b2)))
        XVX <- crossprod(cbind(rep(1, ns), X)*as.vector(sqrt(m*rr/mean(rr))))
        desc <- "poisson.bcs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/(b2^2)
   }
   if (model.int == 17) # (12) poisson.cbs
   {
        m <- opts$m                      # marginal mean of outcome
        b1 <- opts$b1                    # regression coefficient (log rate-ratio) for continuous primary predictor
        b2 <- opts$b2                    # regression coefficient (log rate-ratio) for continuous mediator
        f2 <- opts$f2                    # prevalence of mediator
        rho <- opts$rho                  # correlation of primary predictor and mediator
        sdx1 <- opts$sdx1                # SD of primary predictor
        ns <- opts$ns                    # number of observations in simulated dataset
        n1 <- round(ns*f2)
        n0 <- ns-n1
        mu0 <- -rho*sdx1*sqrt(f2/(1-f2)) # mean of X1 when X2 = 0
        mu1 <- rho*sdx1*sqrt((1-f2)/f2)  # mean of X1 when X2 = 1
        sdx1.2 = sdx1*sqrt(1-rho^2)
        X <-  rbind(cbind(rnorm(n0, mean=mu0, sd=sdx1.2), rep(0, n0)),
                    cbind(rnorm(n1, mean=mu1, sd=sdx1.2), rep(1, n1)))
        rr <- exp(X%*%matrix(c(b1, b2)))
        XVX <- crossprod(cbind(rep(1, ns), X)*as.vector(sqrt(m*rr/mean(rr))))
        desc <- "poisson.cbs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/(b2^2)
   }
   if (model.int == 18) # (12) poisson.bbs
   {
        m <- opts$m                      # marginal mean of outcome
        b1 <- opts$b1                    # regression coefficient (log rate-ratio) for primary predictor
        f1 <- opts$f1                    # prevalence of primary predictor
        b2 <- opts$b2                    # regression coefficient (log rate-ratio) for mediator
        f2 <- opts$f2                    # prevalence of mediator
        rho <- opts$rho                  # correlation of primary predictor and mediator
        ns <- opts$ns                    # number of observations in simulated dataset
        f11 <- rho*sqrt(f1*(1-f1)*f2*(1-f2))+f1*f2
        f10 <- f1-f11
        f01 <- f2-f11
        f00 <- 1-f01-f10-f11
        n11 <- round(f11*ns)
        n01 <- round(f01*ns)
        n10 <- round(f10*ns)
        n00 <- ns-n01-n10-n11
        X <- rbind(cbind(rep(0, n00), rep(0, n00)),
                   cbind(rep(1, n10), rep(0, n10)),
                   cbind(rep(0, n01), rep(1, n01)),
                   cbind(rep(1, n11), rep(1, n11)))
        rr <- exp(X%*%matrix(c(b1, b2)))
        XVX <- crossprod(cbind(rep(1, ns), X)*as.vector(sqrt(m*rr/mean(rr))))
        desc <- "poisson.bbs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/(b2^2)
   }
   if (model.int == 19) # (13) cox.approx
   {
        b2 <- opts$b2     # regression coefficient (log hazard ratio) for mediator
        rho <- opts$rho   # correlation of primary predictor and mediator
        f <- opts$f       # fraction of follow-up times that are uncensored
        sdx2 <- opts$sdx2 # SD of mediator; use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
        desc <- "cox.approx"
        d <- (qnorm(alpha)+qnorm(gamma))^2/((b2*sdx2)^2*(1-rho^2))
        n <- d/f
   }
   if (model.int == 20) # (14) cox.ccs
   {
        b1 <- opts$b1      # regression coefficient (log hazard ratio) for primary predictor
        b2 <- opts$b2      # regression coefficient (log hazard ratio) for mediator
        rho <- opts$rho    # correlation of primary predictor and mediator
        f <- opts$f        # fraction of follow-up times that are uncensored
        sdx1 <- opts$sdx1  # SD of primary predictor
        sdx2 <- opts$sdx2  # SD of mediator
        ns <- opts$ns      # number of observations in simulated dataset
        S <- matrix(c(sdx1^2, rho*sdx1*sdx2, rho*sdx1*sdx2, sdx2^2), 2, 2)
        X <- matrix(rnorm(2*ns), ns, 2) %*% chol(S)
        time <- rexp(ns, rate = exp(X%*%as.matrix(c(b1, b2))))
        event <- ifelse(rank(time)<=round(ns*f), 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.ccs"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int == 21) # (14) cox.bcs
   {
        b1 <- opts$b1     # regression coefficient (log hazard ratio) for primary predictor
        f1 <- opts$f1     # prevalence of primary predictor
        b2 <- opts$b2     # regression coefficient (log hazard ratio) for mediator
        rho <- opts$rho   # correlation of primary predictor and mediator
        f <- opts$f       # fraction of follow-up times that are uncensored
        sdx2 <- opts$sdx2 # SD of mediator
        ns <- opts$ns     # number of observations in simulated dataset
        mu0 <- -rho*sdx2*sqrt(f1/(1-f1))
        mu1 <- rho*sdx2*sqrt((1-f1)/f1)
        sdx2.1 <- sdx2*sqrt(1-rho^2)
        n1 <- round(f1*ns)
        n0 <- ns-n1
        X <- rbind(cbind(rep(0, n0), rnorm(n0, mean=mu0, sd=sdx2.1)),
                   cbind(rep(1, n1), rnorm(n1, mean=mu1, sd=sdx2.1)))
        time <- rexp(ns, rate = exp(X%*%as.matrix(c(b1, b2))))
        event <- ifelse(rank(time)<=round(ns*f), 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.bcs"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int == 22) # (14) cox.cbs
   {
        b1 <- opts$b1     # regression coefficient (log hazard ratio) for primary predictor
        b2 <- opts$b2     # regression coefficient (log hazard ratio) for mediator
        f2 <- opts$f2     # prevalence of mediator
        rho <- opts$rho   # correlation of primary predictor and mediator
        f <- opts$f       # fraction of follow-up times that are uncensored
        sdx1 <- opts$sdx1 # SD of primary predictor
        ns <- opts$ns     # number of observations in simulated dataset
        mu0 <- -rho*sdx1*sqrt(f2/(1-f2))
        mu1 <- rho*sdx1*sqrt((1-f2)/f2)
        sdx1.2 <- sdx1*sqrt(1-rho^2)
        n1 <- round(f2*ns)
        n0 <- ns-n1
        X <- rbind(cbind(rnorm(n0, mean=mu0, sd=sdx1.2), rep(0, n0)),
                   cbind(rnorm(n1, mean=mu1, sd=sdx1.2), rep(1, n1)))
        time <- rexp(ns, rate = exp(X%*%as.matrix(c(b1, b2))))
        event <- ifelse(rank(time)<=round(ns*f), 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.cbs"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int == 23) # (14) cox.bbs
   {
        b1 <- opts$b1   # regression coefficient (log hazard ratio) for primary predictor
        f1 <- opts$f1   # prevalence of primary predictor
        b2 <- opts$b2   # regression coefficient (log hazard ratio) for mediator
        f2 <- opts$f2   # prevalence of mediator
        rho <- opts$rho # correlation of primary predictor and mediator
        f <- opts$f     # fraction of follow-up times that are uncensored
        ns <- opts$ns   # number of observations in simulated dataset
        f11 <- rho*sqrt(f1*(1-f1)*f2*(1-f2))+f1*f2
        n11 <- round(ns*f11)
        n10 <- round(ns*(f1-f11))
        n01 <- round(ns*(f2-f11))
        n00 <- ns-n10-n01-n11
        X <- rbind(cbind(rep(0, n00), rep(0, n00)),
                   cbind(rep(1, n10), rep(0, n10)),
                   cbind(rep(0, n01), rep(1, n01)),
                   cbind(rep(1, n11), rep(1, n11)))
        time <- rexp(ns, rate = exp(X%*%as.matrix(c(b1, b2))))
        event <- ifelse(rank(time)<=round(ns*f), 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.bbs"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int == 24) # (14) cox.ccs2
   {
        b1 <- opts$b1      # regression coefficient (log hazard ratio) for primary predictor
        b2 <- opts$b2      # regression coefficient (log hazard ratio) for mediator
        rho <- opts$rho    # correlation of primary predictor and mediator
        f <- opts$f        # fraction of follow-up times that are uncensored
        fc <- opts$fc      # fraction of failure times that are censored before end of study
        sdx1 <- opts$sdx1  # SD of primary predictor
        sdx2 <- opts$sdx2  # SD of mediator
        ns <- opts$ns      # number of observations in simulated dataset
        S <- matrix(c(sdx1^2, rho*sdx1*sdx2, rho*sdx1*sdx2, sdx2^2), 2, 2)
        X <- matrix(rnorm(2*ns), ns, 2) %*% chol(S)
        re <- exp(X%*%as.matrix(c(b1, b2)))
        te <- rexp(ns, re)
        if (fc>0) {
           tc <- rexp(ns, mean(re)*fc/f)
        } else {
           tc <- rep(max(te)+1, ns)
        }
        time <- apply(cbind(te, tc), 1, min)
        event <- ifelse(rank(time)<=round(ns*(f+fc)) & te<tc, 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.ccs2"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int == 25) # (14) cox.bcs2
   {
        b1 <- opts$b1     # regression coefficient (log hazard ratio) for primary predictor
        f1 <- opts$f1     # prevalence of primary predictor
        b2 <- opts$b2     # regression coefficient (log hazard ratio) for mediator
        rho <- opts$rho   # correlation of primary predictor and mediator
        f <- opts$f       # fraction of follow-up times that are uncensored
        fc <- opts$fc     # fraction of follow-up times that are censored early
        sdx2 <- opts$sdx2 # SD of mediator
        ns <- opts$ns     # number of observations in simulated dataset
        mu0 <- -rho*sdx2*sqrt(f1/(1-f1))
        mu1 <- rho*sdx2*sqrt((1-f1)/f1)
        sdx2.1 <- sdx2*sqrt(1-rho^2)
        n1 <- round(f1*ns)
        n0 <- ns-n1
        X <- rbind(cbind(rep(0, n0), rnorm(n0, mean=mu0, sd=sdx2.1)),
                   cbind(rep(1, n1), rnorm(n1, mean=mu1, sd=sdx2.1)))
        re <- exp(X%*%as.matrix(c(b1, b2)))
        te <- rexp(ns, re)
        if (fc>0) {
           tc <- rexp(ns, mean(re)*fc/f)
        } else {
           tc <- rep(max(te)+1, ns)
        }
        time <- apply(cbind(te, tc), 1, min)
        event <- ifelse(rank(time)<=round(ns*(f+fc)) & te<tc, 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.bcs2"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int == 26) # (14) cox.cbs2
   {
        b1 <- opts$b1     # regression coefficient (log hazard ratio) for primary predictor
        b2 <- opts$b2     # regression coefficient (log hazard ratio) for mediator
        f2 <- opts$f2     # prevalence of mediator
        rho <- opts$rho   # correlation of primary predictor and mediator
        f <- opts$f       # fraction of follow-up times that are uncensored
        fc <- opts$fc     # fraction of follow-up times that are censored early
        sdx1 <- opts$sdx1 # SD of primary predictor
        ns <- opts$ns     # number of observations in simulated dataset
        mu0 <- -rho*sdx1*sqrt(f2/(1-f2))
        mu1 <- rho*sdx1*sqrt((1-f2)/f2)
        sdx1.2 <- sdx1*sqrt(1-rho^2)
        n1 <- round(f2*ns)
        n0 <- ns-n1
        X <- rbind(cbind(rnorm(n0, mean=mu0, sd=sdx1.2), rep(0, n0)),
                   cbind(rnorm(n1, mean=mu1, sd=sdx1.2), rep(1, n1)))
        re <- exp(X%*%as.matrix(c(b1, b2)))
        te <- rexp(ns, re)
        if (fc>0) {
           tc <- rexp(ns, mean(re)*fc/f)
        } else {
           tc <- rep(max(te)+1, ns)
        }
        time <- apply(cbind(te, tc), 1, min)
        event <- ifelse(rank(time)<=round(ns*(f+fc)) & te<tc, 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.cbs2"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int == 27) # (14) cox.bbs2
   {
        b1 <- opts$b1   # regression coefficient (log hazard ratio) for primary predictor
        f1 <- opts$f1   # prevalence of primary predictor
        b2 <- opts$b2   # regression coefficient (log hazard ratio) for mediator
        f2 <- opts$f2   # prevalence of mediator
        rho <- opts$rho # correlation of primary predictor and mediator
        f <- opts$f     # fraction of follow-up times that are uncensored
        fc <- opts$fc   # fraction of follow-up times that are censored early
        ns <- opts$ns   # number of observations in simulated dataset
        f11 <- rho*sqrt(f1*(1-f1)*f2*(1-f2))+f1*f2
        n11 <- round(ns*f11)
        n10 <- round(ns*(f1-f11))
        n01 <- round(ns*(f2-f11))
        n00 <- ns-n10-n01-n11
        X <- rbind(cbind(rep(0, n00), rep(0, n00)),
                   cbind(rep(1, n10), rep(0, n10)),
                   cbind(rep(0, n01), rep(1, n01)),
                   cbind(rep(1, n11), rep(1, n11)))
        re <- exp(X%*%as.matrix(c(b1, b2)))
        te <- rexp(ns, re)
        if (fc>0) {
           tc = rexp(ns, mean(re)*fc/f)
        } else {
           tc = rep(max(te)+1, ns)
        }
        time <- apply(cbind(te, tc), 1, min)
        event <- ifelse(rank(time)<=round(ns*(f+fc)) & te<tc, 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.bbs2"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int < 19) list(desc=desc, n=round(n))
   else list(desc=desc, d=round(d), n = round(n))
}

getb0 <- function(p, X, b1, b2)
{
  stopifnot(is.numeric(p), length(p) == 1, p > 0, p < 1)
  eta0 <- X[,2] * b1 + X[,3] * b2
  f <- function(b0) mean(plogis(b0+eta0)) - p
  uniroot(f, interval = c(-50,50), tol = 1e-10)$root
}
