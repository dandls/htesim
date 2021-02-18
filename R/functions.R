#' @name treatmentfunc
#' @rdname treatmentfunc
#'
#' @title Treatment effect functions
#'
#' @details \code{tF_exp_x1_x2} corresponds to \eqn{\prod_{j = 1}^2 1 + \frac{1}{1 + \exp{(-20(x_j - \frac{1}{3}))}}}.
#' \code{tF_div_x1_x2} corresponds to \eqn{\frac{x_1 + x_2}{2}}.
#' \code{tF_log_x1_x2} to \eqn{x_1 + log(1 + \exp{x_2})}
#' and \code{tF_max_x1_x5} to \eqn{\max{x_1 + x_2 + x_3, 0} - \max{x_4 + x_5, 0}}.
#'
#' @param x (data.frame)
#'
#' @return vector of corresponding treatment effects
#' @export
tF_exp_x1_x2 <- function(x) {
  return(h_exp(x[,"X1"]) * h_exp(x[,"X2"]))
}

h_exp <- function(x) {
  return(1 + 1 / (1 + exp(-20 * (x - 1/3))))
}

#'@rdname treatmentfunc
tF_div_x1_x2 <- function(x) {
  return((x[,"X1"] + x[,"X2"]) / 2)
}

#'@rdname treatmentfunc
tF_log_x1_x2 <- function(x) {
 return(x[,"X1"] + log(1 + exp(x[,"X2"])))
}

#'@rdname treatmentfunc
tF_max_x1_x5 <- function(x){
  return(pmax(x[,"X1"] + x[,"X2"] + x[,"X3"], 0) - pmax(x[,"X4"] + x[,"X5"], 0))
}

#' @name propensityfunc
#' @rdname propensityfunc
#'
#' @title Treatment propensity functions
#'
#' @param x (data.frame)
#'
#' @return vector of corresponding propensities
#' @export
pF_x1 <- function(x) {
  return(1 / 4 * (1 + dbeta(x[,"X1"], 2, 4)))
}

#'@rdname propensityfunc
pF_x3 <- function(x) {
  return(1 / 4 * (1 + dbeta(x[,"X3"], 2, 4)))
}

#'@rdname propensityfunc
pF_x4 <- function(x) {
  return(1 / 4 * (1 + dbeta(x[,"X4"], 2, 4)))
}

#'@rdname propensityfunc
pF_sin_x3 <- function(x) {
  return(sin(2 * pi * x[,"X3"]) / 4 + .5)
}

#'@rdname propensityfunc
pF_eta_x1_x2 <- function(x, eta = 0.1) {
  return(pmax(eta, pmin(sin(pi * x[,"X1"] * x[,"X2"]), 1-eta)))
}

#'@rdname propensityfunc
pF_x2_x3 <- function(x) {
  return(1/(1 + exp(x[,"X2"] + x[,"X3"])))
}

#'@rdname propensityfunc
pF_exp_x1_x2 <- function(x) {
  return(1/(1 + exp(-x[,"X1"]) + exp(-x[,"X2"])))
}

#'@rdname treatmentfunc
tF_max_x1_x5 <- function(x){
  return(pmax(x[,"X1"] + x[,"X2"] + x[,"X3"], 0) - pmax(x[,"X4"] + x[,"X5"], 0))
}

#' @name prognosticfunc
#' @rdname prognosticfunc
#'
#' @title Prognostic effect functions
#'
#' @param x (data.frame)
#'
#' @return vector of corresponding prognostic effects
#' @export
mF_x1 <- function(x) {
  return(2 * x[,"X1"] - 1)
}

#' @rdname prognosticfunc
mF_x3 <- function(x) {
  return(2 * x[,"X3"] - 1)
}

#' @rdname prognosticfunc
mF_sin_x1_x5 <- function(x) {
  return(sin(pi * x[,"X1"] * x[,"X2"]) + 2*(x[,"X3"]-0.5)^2 + x[,"X4"] + 0.5*x[,"X5"])
}

#' @rdname prognosticfunc
mF_max_x1_x5 <- function(x) {
  return(pmax(0, x[,"X1"] + x[,"X2"], x[,"X3"]) + pmax(0, x[,"X4"] + x[,"X5"]))
}

#' @rdname prognosticfunc
mF_log_x1_x3 <- function(x) {
  return(2 * log(1 + exp(x[,"X1"] + x[,"X2"] + x[,"X3"])))
}

#' @rdname prognosticfunc
mF_max2_x1_x5 <- function(x) {
  return((pmax(x[,"X1"] + x[,"X2"] + x[,"X3"], 0) + pmax(x[,"X4"] + x[,"X5"], 0)) / 2)
}
