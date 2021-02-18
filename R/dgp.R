#' Create manual for data generation
#' @param p (numeric(1), function) Propensities P(trt = 1|x) (pi(x)),
#' default = 0.5 means that there is no confounding.
#' @param m (numeric(1), function) Prognostic effect (mu(x)), default 0 means that
#' there is no prognostic effect.
#' @param t (numeric(1), function) Predictive treatment effect (tau(x)),
#' default 0 means that there is no treatment effect.
#' @param sd (numeric(1)) Standard deviation of normal distribution, only
#' has an effect if model = "normal".
#' @param pi (numeric(1)) how much of predictive effect is added to prognostic effect.
#' If pi = 0, conditional mean does not depend on treatment effect.
#' @param model ("normal"|"weibull") Name of used model.
#' @return list of class gpd with entries:
#' pfct (pi(x)), mfct (mu(x)), tfct (tau(x)), sdfct (sd for normal),
#' model (model used)).
dgp <- function(p = 0.5, m = 0, t = 0, sd = 1, pi = .5,
  model = c("normal", "weibull", "binomial", "polr"), xmodel = c("normal", "unif")) {

  model <- match.arg(model)
  xmodel <- match.arg(xmodel)

  if (!is.function(pfct <- p))
    pfct <- function(x) return(rep(p, NROW(x)))

  if (!is.function(mfct <- m))
    mfct <- function(x) return(rep(m, NROW(x)))

  if (!is.function(tfct <- t))
    tfct <- function(x) return(rep(t, NROW(x)))

  if (!is.function(sdfct <- sd))
    sdfct <- function(x) return(rep(sd, NROW(x)))

  ### mean effect
  mf <- function(x) mfct(x) - pi * tfct(x)

  ret <- list(pfct = pfct, mfct = mf, tfct = tfct, sdfct = sdfct, pi = pi,
    model = model, xmodel = xmodel, ename = deparse(substitute(e)),
    mname = deparse(substitute(m)), zname = deparse(substitute(t)),
    sdname = deparse(substitute(sd)))
  class(ret) <- "dgp"
  return(ret)
}
