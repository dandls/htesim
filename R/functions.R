#' @name tF
#' @rdname tF
#'
#' @title Treatment effect functions
#'
#' @description Treatment effect functions are based on Wager and Athey (2018) and Nie and Wager (2020).
#'
#' @details
#' \code{tF_exp_x1_x2} corresponds to \deqn{\prod_{j = 1}^2 1 + \frac{1}{1 + \exp{(-20(x_j - \frac{1}{3}))}}}{\prod_{j = 1}^2 (1 +  1/(1 + exp(-20 (xj - (1/3))))).}
#' \code{tF_exp2_x1_x2} is equal to \deqn{\prod_{j = 1}^2 \frac{2}{1 + \exp{(-20(x_j - \frac{1}{2}))}}}{\prod_{j = 1}^2 (2/(1 + exp(-20 (xj - (1/2))))).}
#' \code{tF_div_x1_x2} corresponds to \deqn{\frac{x_1 + x_2}{2}}{(x1 + x2)/2,}
#' \code{tF_log_x1_x2} to \deqn{x_1 + \log{1 + \exp{x_2}}}{x1 + log(1 + x2)}
#' and \code{tF_max_x1_x5} to \deqn{\max{x_1 + x_2 + x_3, 0} - \max{x_4 + x_5, 0}}{max(x1 + x2 + x3, 0) - max(x4 + x5, 0).}
#'
#' @param x (data.frame)
#' @param h (function) Helper function.
#'
#' @return vector of corresponding treatment effects
#'
#' @references
#' Wager S, and Athey S (2018). "Estimation and Inference of Heterogeneous Treatment Effects using Random Forests". Journal of the American Statistical Association, 113(523).
#'
#' Nie X, Wager S (2020). “Quasi-Oracle Estimation of Heterogeneous Treatment Ef- fects.” Biometrika. ISSN 0006-3444. doi:10.1093/biomet/asaa076. Asaa076, https://academic.oup.com/biomet/advance-article-pdf/doi/10.1093/biomet/asaa076/33788449/asaa076.pdf, URL https://doi.org/10.1093/biomet/asaa076.
#'
#' @export
tF_exp_x1_x2 <- function(x) {
  return(h_exp(x[,"X1"]) * h_exp(x[,"X2"]))
}

h_exp <- function(x) {
  return(1 + 1 / (1 + exp(-20 * (x - 1/3))))
}

#' @rdname tF
#' @export
tF_exp2_x1_x2 <- function(x) {
  return(h2_exp(x[,"X1"]) * h2_exp(x[,"X2"]))
}

h2_exp <- function(x) {
  return(2 / (1 + exp(-20 * (x - 1/2))))
}

#' @rdname tF
#' @export
tF_div_x1_x2 <- function(x) {
  return((x[,"X1"] + x[,"X2"]) / 2)
}

#' @rdname tF
#' @export
tF_log_x1_x2 <- function(x) {
 return(x[,"X1"] + log(1 + exp(x[,"X2"])))
}

#' @rdname tF
#' @export
tF_max_x1_x5 <- function(x){
  return(pmax(x[,"X1"] + x[,"X2"] + x[,"X3"], 0) - pmax(x[,"X4"] + x[,"X5"], 0))
}

#' @name pF
#' @rdname pF
#'
#' @title Treatment propensity functions
#'
#' @description Treatment propensity functions are based on Wager and Athey (2018) and Nie and Wager (2020).
#'
#' @details
#' \code{pF_x1} corresponds to \deqn{\frac{1}{4}(1 + \beta_{2, 4}(x_1))}{1/4 (1 + \beta(x1|2, 4)).}
#' Equivalently, \code{pF_x3} corresponds to \deqn{\frac{1}{4}(1 + \beta_{2, 4}(x_3))}{1/4 (1 + \beta(x3|2, 4))}
#' and \code{pF_x4} to \deqn{\frac{1}{4}(1 + \beta_{2, 4}(x_4))}{1/4 (1 + \beta(x4|2, 4)).}
#' \code{pF_sin_x3} is equal to \deqn{\sin{2 * pi * x_3}/4 + 0.5}{sin(2 pi x3)/4 + 0.5,}
#' \code{pF_eta_x1_x2} to \deqn{\max{\eta, \min{\sin{\pi x_1 x_2}, 1 - \eta}}}{max(\eta, min(sin(\pi x1 x2), 1 - \eta)),}
#' \code{pF_x2_x3} to \deqn{1/(1 + \exp{x_2 + x_3})}{1/(1 + exp(x2 + x3))}
#' and \code{pF_exp_x1_x2} to \deqn{1/(1 + \exp{- x_1} + \exp{- x_2})}{1/(1 + exp(-x1) + exp(-x2)).}
#'
#' @param x (data.frame)
#' @param eta (numeric(1))  Trimming parameter for \code{pF_eta_x1_x2}.
#'
#' @return vector of corresponding propensities
#'
#' @references
#' Wager S, and Athey S (2018). "Estimation and Inference of Heterogeneous Treatment Effects using Random Forests". Journal of the American Statistical Association, 113(523).
#'
#' Nie X, Wager S (2020). “Quasi-Oracle Estimation of Heterogeneous Treatment Ef- fects.” Biometrika. ISSN 0006-3444. doi:10.1093/biomet/asaa076. Asaa076, https://academic.oup.com/biomet/advance-article-pdf/doi/10.1093/biomet/asaa076/33788449/asaa076.pdf, URL https://doi.org/10.1093/biomet/asaa076.
#'
#' @export
pF_x1 <- function(x) {
  return(1 / 4 * (1 + dbeta(x[,"X1"], 2, 4)))
}

#' @rdname pF
#' @export
pF_x3 <- function(x) {
  return(1 / 4 * (1 + dbeta(x[,"X3"], 2, 4)))
}

#' @rdname pF
#' @export
pF_x4 <- function(x) {
  return(1 / 4 * (1 + dbeta(x[,"X4"], 2, 4)))
}

#' @rdname pF
#' @export
pF_sin_x3 <- function(x) {
  return(sin(2 * pi * x[,"X3"]) / 4 + .5)
}

#' @rdname pF
#' @export
pF_eta_x1_x2 <- function(x, eta = 0.1) {
  return(pmax(eta, pmin(sin(pi * x[,"X1"] * x[,"X2"]), 1-eta)))
}

#' @rdname pF
#' @export
pF_x2_x3 <- function(x) {
  return(1/(1 + exp(x[,"X2"] + x[,"X3"])))
}

#' @rdname pF
#' @export
pF_exp_x1_x2 <- function(x) {
  return(1/(1 + exp(-x[,"X1"]) + exp(-x[,"X2"])))
}

#' @name mF
#' @rdname mF
#'
#' @title Prognostic effect functions
#'
#' @description Prognostic effect functions are based on Wager and Athey (2018) and Nie and Wager (2020).
#'
#' @details
#' \code{mF_x1} corresponds to \deqn{2x_1 - 1}{2 x_1 - 1.}
#' Equivalently, \code{mF_x3} corresponds to \deqn{2x_3 - 1}{2 x3 - 1.}
#' \code{mF_sin_x1_x5} is equal to \deqn{\sin{\pi x_1 x_2} + 2(x_3 - 0.5)^2 + x_4 + 0.5 x_5}{sin(\pi x1 x2) + 2 (x3 - 0.5)^2 + x4 + 0.5 x5,}
#' \code{mF_max_x1_x5} to \deqn{\max{x_1 + x_2, x_3, 0} + \max{x_4 + x_5, 0}}{max(x1 + x2 + x3, 0) + max(x4 + x5, 0),}
#' \code{mF_log_x1_x3} to \deqn{2 \log{1 + exp(x_1 + x_2 + x_3)}}{2 log(1 + x1 + x2 + x3),}
#' and \code{mF_max2_x1_x5} to \deqn{(\max{x_1 + x_2 + x_3, 0} + \max{x_4 + x_5, 0})/2}{(max(x1 + x2 + x3, 0) + max(x4 + x5, 0))/2.}
#'
#' @param x (data.frame)
#'
#' @return vector of corresponding prognostic effects
#'
#' @references
#' Wager S, and Athey S (2018). "Estimation and Inference of Heterogeneous Treatment Effects using Random Forests". Journal of the American Statistical Association, 113(523).
#'
#' Nie X, Wager S (2020). “Quasi-Oracle Estimation of Heterogeneous Treatment Ef- fects.” Biometrika. ISSN 0006-3444. doi:10.1093/biomet/asaa076. Asaa076, https://academic.oup.com/biomet/advance-article-pdf/doi/10.1093/biomet/asaa076/33788449/asaa076.pdf, URL https://doi.org/10.1093/biomet/asaa076.
#'
#' @export
mF_x1 <- function(x) {
  return(2 * x[,"X1"] - 1)
}

#' @rdname mF
#' @export
mF_x3 <- function(x) {
  return(2 * x[,"X3"] - 1)
}

#' @rdname mF
#' @export
mF_sin_x1_x5 <- function(x) {
  return(sin(pi * x[,"X1"] * x[,"X2"]) + 2*(x[,"X3"]-0.5)^2 + x[,"X4"] + 0.5*x[,"X5"])
}

#' @rdname mF
#' @export
mF_max_x1_x5 <- function(x) {
  return(pmax(0, x[,"X1"] + x[,"X2"], x[,"X3"]) + pmax(0, x[,"X4"] + x[,"X5"]))
}

#' @rdname mF
#' @export
mF_log_x1_x3 <- function(x) {
  return(2 * log(1 + exp(x[,"X1"] + x[,"X2"] + x[,"X3"])))
}

#' @rdname mF
#' @export
mF_max2_x1_x5 <- function(x) {
  return((pmax(x[,"X1"] + x[,"X2"] + x[,"X3"], 0) + pmax(x[,"X4"] + x[,"X5"], 0)) / 2)
}
