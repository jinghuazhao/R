#' Fit AE, ACE or ADE biometric mixed models to nuclear family data
#'
#' Fits classical biometric variance-decomposition models using a
#' linear mixed model formulation implemented with `nlme::lme`.
#'
#' The function estimates additive genetic (A), shared environmental (C),
#' dominance genetic (D), and unique environmental (E) variance components
#' from nuclear family data (parents and offspring).
#'
#' The model is valid for:
#' * Parent-child trios (one offspring)
#' * Nuclear families with multiple offsprings
#'
#' Only AE, ACE and ADE models are fitted.
#' A second family-level random effect represents either shared
#' environmental (C) or dominance genetic (D) variance and,
#' in trio-only data, ACE and ADE models are statistically indistinguishable.
#'
#' @param model   Fixed-effects formula.
#' @param data    Nuclear-family data in long format (one row per individual).
#' @param type    Biometric model: "AE", "ACE" or "ADE".
#' @param method  "ML" (default) or "REML".
#'
#' Required columns in `data`:
#' \describe{
#'   \item{familyid}{Nuclear family identifier}
#'   \item{var1}{Mother indicator (1/0)}
#'   \item{var2}{Father indicator (1/0)}
#'   \item{var3}{Child indicator (1/0)}
#' }
#'
#' @return A list with:
#' \describe{
#'   \item{fit}{Fitted `nlme::lme` model}
#'   \item{var}{Variance components A, C, D, E}
#'   \item{h2}{Narrow-sense heritability (\eqn{A/V_P})}
#'   \item{c2}{Shared environmental proportion (\eqn{C/V_P})}
#'   \item{d2}{Dominance proportion (\eqn{D/V_P})}
#'   \item{H2}{Broad-sense heritability (\eqn{(A+D)/V_P}) in the ADE model;
#'             in the AE model, \code{h2} provides an upper bound on total
#'             genetic influence}
#' }
#'
#' @details
#' The function automatically detects whether the dataset contains
#' only parent–child trios or nuclear families with multiple offspring
#' and adapts the genetic parameterisation accordingly.
#'
#' \strong{Trio data (one offspring per family)}
#'
#' Additive genetic variance is represented using a single family-level
#' random effect with coefficient
#'
#' \deqn{A = 0.5\,Mother + 0.5\,Father + 1\,Child.}
#'
#' This parameterisation is sufficient when only one offspring is present.
#'
#' \strong{Sibling data (multiple offspring per family)}
#'
#' Additive genetic variance is decomposed into three independent
#' transmitted components:
#'
#' \itemize{
#'   \item transmitted maternal genes (Am)
#'   \item transmitted paternal genes (Af)
#'   \item Mendelian sampling of the offspring (Ms)
#' }
#'
#' For an offspring phenotype,
#' \deqn{A = A_m + A_f + M_s}
#'
#' and the additive genetic variance becomes
#' \deqn{Var(A) = 2(\sigma^2_{Am} + \sigma^2_{Af}) + \sigma^2_{Ms}.}
#'
#' A second family-level random effect represents either:
#' \itemize{
#'   \item shared environmental variance (C) in the ACE model, or
#'   \item dominance genetic variance (D) in the ADE model.
#' }
#'
#' With trio-only data, ACE and ADE models have identical likelihood and
#' cannot be statistically distinguished. The function nevertheless allows
#' both parameterisations to be fitted for completeness and model comparison.
#'
#' @note
#' Because shared environmental (C) and dominance genetic (D) effects are not
#' separately identifiable in nuclear family data, AE, ACE and ADE models
#' should be interpreted jointly rather than in isolation.
#'
#' The AE model typically provides an upper bound on total genetic influence.
#' The ADE model provides an estimate of broad-sense heritability when
#' dominance effects are supported by the data.
#'
#' The reported \code{H2} equals \eqn{(A + D)/V_P}. In AE models, where D is not
#' estimated, \code{h2} should be interpreted as an upper bound on total
#' genetic variance.
#'
#' This function complements \code{fbsize()} and \code{pbsize()} by enabling
#' biometric variance decomposition using linear mixed models for nuclear
#' family data.
#'
#' @references
#' \insertAllCited{}
#' \insertRef{rh08}{gap}
#'
#' @examples
#' library(gap.datasets)
#'
#' model <- bwt ~ male + first + midage + highage + birthyr
#'
#' AE  <- ACDE(model, mfblong, "AE")
#' ACE <- ACDE(model, mfblong, "ACE")
#' ADE <- ACDE(model, mfblong, "ADE")
#'
#' ACE$h2
#' tab <- anova(AE$fit, ACE$fit, ADE$fit)
#' res <- list(AE=AE, ACE=ACE, ADE=ADE, anova=tab)
#' tab <- rbind(
#'   AE  = c(res$AE$var,  h2=res$AE$h2,  c2=res$AE$c2,  d2=res$AE$d2),
#'   ACE = c(res$ACE$var, h2=res$ACE$h2, c2=res$ACE$c2, d2=res$ACE$d2),
#'   ADE = c(res$ADE$var, h2=res$ADE$h2, c2=res$ADE$c2, d2=res$ADE$d2)
#' )
#' tab <- round(tab, digits=3)
#' list(
#'   fits    = res,
#'   anova   = res$anova,
#'   summary = tab
#' )
#'
#' @author ChatGPT
#' @export
#'
ACDE <- function(model, data, type=c("AE","ACE","ADE"), method="ML")
{
  if(!requireNamespace("nlme", quietly=TRUE))
    stop("Package 'nlme' is required.")

  type <- match.arg(type)
  if(!all(c("var1","var2","var3") %in% names(data)))
    stop("Data must contain var1 (mother), var2 (father), var3 (child) indicators")

  child_count <- tapply(data$var3, data$familyid, sum)
  has_siblings <- any(child_count > 1)
  if(!has_siblings)
  {
    message("Trio data detected --> using collapsed additive model")
    data$Acoef <- 0.5*data$var1 + 0.5*data$var2 + 1*data$var3
    data$Fcoef <- 1
    if(type=="AE")
      rand <- list(familyid = nlme::pdDiag(~ Acoef - 1))
    if(type %in% c("ACE","ADE"))
      rand <- list(familyid = nlme::pdDiag(~ Acoef + Fcoef - 1))
    param <- "trio"
  } else {
    message("Sibling data detected --> using transmission decomposition")
    data$Amcoef <- data$var1
    data$Afcoef <- data$var2
    data$Mscoef <- data$var3
    data$Fcoef  <- 1
    if(type=="AE")
      rand <- list(familyid = nlme::pdDiag(~ Amcoef + Afcoef + Mscoef -1))
    if(type %in% c("ACE","ADE"))
      rand <- list(familyid = nlme::pdDiag(~ Amcoef + Afcoef + Mscoef + Fcoef -1))
    param <- "sibling"
  }
  fit <- nlme::lme(model,
                   random = rand,
                   data   = data,
                   method = method,
                   control = nlme::lmeControl(opt="optim"))
  vc <- nlme::VarCorr(fit)
  sd_vals <- as.numeric(vc[,"StdDev"])
  varE <- fit$sigma^2
  re_var <- sd_vals[1:(length(sd_vals)-1)]^2
  if(param=="trio")
  {
    varA <- re_var[1]
    varF <- if(type!="AE") re_var[2] else 0
  } else {
    varAm <- re_var[1]
    varAf <- re_var[2]
    varMs <- re_var[3]
    varF  <- if(type!="AE") re_var[4] else 0
    varA  <- 2*(varAm + varAf) + varMs
  }
  varC <- if(type=="ACE") varF else 0
  varD <- if(type=="ADE") varF else 0
  varP <- varA + varC + varD + varE
  list(
    fit = fit,
    var = c(A=varA, C=varC, D=varD, E=varE),
    h2  = varA/varP,
    c2  = varC/varP,
    d2  = varD/varP,
    H2  = (varA + varD)/varP
  )
}
