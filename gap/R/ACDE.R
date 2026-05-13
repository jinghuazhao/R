#' Fit AE, ACE or ADE biometric mixed models to nuclear family data
#'
#' Fits classical biometric variance–decomposition models using a linear
#' mixed model formulation implemented with `nlme::lme`.
#'
#' The function estimates additive genetic (A), shared environmental (C),
#' dominance genetic (D), and unique environmental (E) variance components
#' from **nuclear family data (parents and offspring)** in long format.
#'
#' Supported family structures:
#' * Parent-child trios (one offspring)
#' * Nuclear families with **any number of siblings**
#'
#' The function automatically detects the family structure and selects the
#' correct additive-genetic parameterisation.
#'
#' Only AE, ACE and ADE models are fitted.
#'
#' @param model  Fixed-effects formula.
#' @param data   Data frame in long format (one row per individual).
#' @param type   Model type: `"AE"`, `"ACE"` or `"ADE"`.
#' @param method Estimation method: `"ML"` (default) or `"REML"`.
#'
#' @section Required columns in `data`:
#'
#' \describe{
#' \item{familyid}{Nuclear family identifier.}
#' \item{var1}{Maternal transmission coefficient.}
#' \item{var2}{Paternal transmission coefficient.}
#' \item{var3}{Offspring (Mendelian sampling) indicator.}
#' }
#'
#' These variables encode expected genetic transmission and **are not role
#' indicators**.
#'
#' Coding for a nuclear family:
#'
#' \preformatted{
#' role     var1  var2  var3
#' -------------------------
#' father     0     1     0
#' mother     1     0     0
#' child      1     1     1
#' }
#'
#' For **multiple siblings**, each offspring receives identical coding:
#'
#' \preformatted{
#' role     var1  var2  var3
#' -------------------------
#' father     0     1     0
#' mother     1     0     0
#' sib1       1     1     1
#' sib2       1     1     1
#' sib3       1     1     1
#' }
#'
#' @details
#'
#' Phenotypic variance is decomposed as
#'
#' \deqn{V_P = V_A + V_C + V_D + V_E}
#'
#' The model is fitted as a linear mixed model with family-level random
#' effects and individual residual variance.
#'
#' ## Automatic family detection
#'
#' The function detects whether families contain one or multiple offspring.
#'
#' * **Trios:** collapsed additive parameterisation.
#' * **Siblings:** transmission decomposition parameterisation.
#'
#' ## Trio parameterisation
#'
#' Additive effects represented as:
#'
#' \deqn{A = 0.5 Mother + 0.5 Father + 1 Child}
#'
#' This reproduces the expected parent–offspring covariance:
#' \deqn{Cov = 1/2 V_A}
#'
#' ## Multi-sibling parameterisation
#'
#' Additive genetic variance is decomposed into:
#'
#' * maternal transmission (\eqn{A_m})
#' * paternal transmission (\eqn{A_f})
#' * Mendelian sampling (\eqn{M_s})
#'
#' For offspring:
#' \deqn{A = A_m + A_f + M_s}
#'
#' Total additive variance:
#'
#' \deqn{V_A = 2(\sigma^2_{Am} + \sigma^2_{Af}) + \sigma^2_{Ms}}
#'
#' This produces correct covariances:
#'
#' * Parent–offspring: \eqn{1/2 V_A}
#' * Sibling–sibling:  \eqn{1/2 V_A}
#'
#' This formulation generalises to **any number of siblings**.
#'
#' ## Identifiability of C and D
#'
#' Nuclear family data cannot fully separate shared environment (C)
#' and dominance (D). ACE and ADE models should be interpreted jointly.
#'
#' @note This complements [pbsize()] and [fbsize()].
#' @return Object of class `"ACDEfit"` containing:
#' \describe{
#' \item{fit}{`nlme::lme` object}
#' \item{var}{Variance components (A,C,D,E)}
#' \item{h2}{Narrow-sense heritability}
#' \item{c2}{Shared environment (ACE only)}
#' \item{d2}{Dominance (ADE only)}
#' \item{H2}{Broad-sense heritability (ADE only)}
#' }
#'
#' @examples
#' library(nlme)
#'
#' set.seed(1)
#'
#' simulate_families <- function(n_fam = 200)
#' {
#'   VA <- 0.4; VC <- 0.2; VD <- 0.1; VE <- 0.3
#'   out <- list()
#'
#'   for(f in 1:n_fam){
#'     Cfam <- rnorm(1,0,sqrt(VC))
#'     Af <- rnorm(1,0,sqrt(VA))
#'     Am <- rnorm(1,0,sqrt(VA))
#'
#'     make_child <- function(){
#'       Mend <- rnorm(1,0,sqrt(0.5*VA))
#'       A  <- 0.5*(Af+Am)+Mend
#'       D  <- rnorm(1,0,sqrt(VD))
#'       E  <- rnorm(1,0,sqrt(VE))
#'       A + D + Cfam + E
#'     }
#'
#'     out[[f]] <- data.frame(
#'       familyid=f,
#'       role=c("father","mother","sib1","sib2","sib3"),
#'       y=c(
#'         Af + Cfam + rnorm(1,0,sqrt(VE)),
#'         Am + Cfam + rnorm(1,0,sqrt(VE)),
#'         make_child(), make_child(), make_child()
#'       )
#'     )
#'   }
#'
#'   dat <- do.call(rbind,out)
#'
#'   dat$var1 <- as.integer(dat$role=="mother")
#'   dat$var2 <- as.integer(dat$role=="father")
#'   dat$var3 <- as.integer(grepl("sib", dat$role))
#'   dat
#' }
#'
#' dat <- simulate_families()
#'
#' AE  <- ACDE(y~1, dat, "AE")
#' ACE <- ACDE(y~1, dat, "ACE")
#' ADE <- ACDE(y~1, dat, "ADE")
#'
#' anova(AE$fit, ACE$fit, ADE$fit)
#'
#' ############################################################
#' # Create rectangular variance table (important!)
#' ############################################################
#'
#' summary(AE)
#' summary(ACE)
#' summary(ADE)
#'
#' @author ChatGPT
#' @export
#'
ACDE <- function(model, data, type=c("AE","ACE","ADE"), method="ML")
{
  if(!requireNamespace("nlme", quietly=TRUE))
    stop("Package 'nlme' is required.")
  type <- match.arg(type)

  if(!all(c("familyid","var1","var2","var3") %in% names(data)))
    stop("Data must contain familyid, var1, var2, var3")

  child_count  <- tapply(data$var3, data$familyid, sum)
  has_siblings <- any(child_count > 1)

  if(!has_siblings){
    data$Acoef <- 0.5*data$var1 + 0.5*data$var2 + data$var3
    data$Fcoef <- 1
    rand <- if(type=="AE")
      list(familyid = nlme::pdDiag(~ Acoef -1))
    else
      list(familyid = nlme::pdDiag(~ Acoef + Fcoef -1))
    param <- "trio"
  } else {
    data$Amcoef <- data$var1
    data$Afcoef <- data$var2
    data$Mscoef <- data$var3
    data$Fcoef  <- 1
    rand <- if(type=="AE")
      list(familyid = nlme::pdDiag(~ Amcoef+Afcoef+Mscoef -1))
    else
      list(familyid = nlme::pdDiag(~ Amcoef+Afcoef+Mscoef+Fcoef -1))
    param <- "sibling"
  }

  fit <- nlme::lme(model, random=rand, data=data, method=method,
                   control=nlme::lmeControl(opt="optim"))

  vc <- nlme::VarCorr(fit)
  sd_vals <- as.numeric(vc[,"StdDev"])
  varE <- fit$sigma^2
  re_var <- sd_vals[1:(length(sd_vals)-1)]^2

  if(param=="trio"){
    varA <- re_var[1]
    varF <- if(type!="AE") re_var[2] else 0
  } else {
    varA <- 2*(re_var[1]+re_var[2]) + re_var[3]
    varF <- if(type!="AE") re_var[4] else 0
  }

  varC <- if(type=="ACE") varF else 0
  varD <- if(type=="ADE") varF else 0
  varP <- varA + varC + varD + varE

  h2 <- varA/varP
  c2 <- if(type=="ACE") varC/varP else NA
  d2 <- if(type=="ADE") varD/varP else NA
  H2 <- if(type=="ADE") (varA+varD)/varP else NA

  out <- list(
    call = match.call(),
    type = type,
    parameterisation = param,
    fit = fit,
    var = c(A=varA, C=varC, D=varD, E=varE),
    h2 = h2,
    c2 = c2,
    d2 = d2,
    H2 = H2
  )

  class(out) <- "ACDEfit"
  out
}

#' @export
print.ACDEfit <- function(x, digits = 3, ...)
{
  cat("\nFamily ACDE variance decomposition\n")
  cat("-----------------------------------\n")
  cat("Model:", x$type, "\n")
  cat("Parameterisation:", x$parameterisation, "\n\n")

  v <- round(x$var, digits)

  cat("Variance components\n")
  cat(" A =", v["A"], "\n")
  cat(" C =", v["C"], "\n")
  cat(" D =", v["D"], "\n")
  cat(" E =", v["E"], "\n\n")

  cat("Variance proportions\n")
  cat(" h2 =", round(x$h2, digits), "\n")

  if(!is.na(x$c2))
    cat(" c2 =", round(x$c2, digits), "\n")

  if(!is.na(x$d2))
    cat(" d2 =", round(x$d2, digits), "\n")

  if(!is.na(x$H2))
    cat(" H2 =", round(x$H2, digits), "\n")

  invisible(x)
}

#' @export
summary.ACDEfit <- function(object, ...)
{
  v <- object$var
  Vtot <- sum(v)

  out <- list(
    type = object$type,
    parameterisation = object$parameterisation,
    variance = v,
    proportions = c(
      h2 = object$h2,
      c2 = object$c2,
      d2 = object$d2,
      H2 = object$H2
    ),
    total_variance = Vtot,
    logLik = as.numeric(logLik(object$fit))
  )

  class(out) <- "summary.ACDEfit"
  out
}

#' @export
print.summary.ACDEfit <- function(x, digits = 3, ...)
{
  cat("\nSummary of", x$type, "model\n")
  cat("-----------------------------------\n")
  cat("Parameterisation:", x$parameterisation, "\n")

  cat("\nVariance components:\n")
  print(round(x$variance, digits))

  cat("\nVariance proportions:\n")
  props <- x$proportions[!is.na(x$proportions)]
  print(round(props, digits))

  cat("\nTotal phenotypic variance:", round(x$total_variance, digits), "\n")
  cat("logLik:", round(x$logLik, digits), "\n")

  invisible(x)
}
