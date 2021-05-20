#' Study Setting and Results of Dandl et al. (2021)
#'
#' A dataset containing the study setting and performance results of 8 different
#' tree-based methods for estimating treatment effects based on Dandl et al. (2021).
#' Mean squared error of estimated to true treatment effect measured the performance.
#' Study settings were based on Wager and Athey (2018) and Nie and Wager (2020).
#'
#' @format A data frame with 9600 rows and 23 variables:
#' \describe{
#'   \item{setup}{Identifier for data generating process defined by p, m, t, sd, pi, model, xmodel.}
#'   \item{repl}{Replication index (1 to 100).}
#'   \item{p, m, t}{Propensity score, prognostic effect and treatment effect
#'   function character names. Input for \code{dgp}.}
#'   \item{sd, pi, model, xmodel}{Standard deviation of
#'   normal distributed outcome, extend of treatment effect added to prognostic
#'   effect and names of used model for outcome and covariates.
#'   Input for \code{dgp}.}
#'   \item{nsim, dim, nsimtest, seed}{Number of observations, number of
#'   covariates, number of observations for test set and seed. Input for
#'   \code{simulate}.}
#'   \item{m4y, m4yhonest, hybrid, hybridhonest, equalized, equalizedhonest,
#'   cf, cfhonest, m4ycf, m4ycfhonest}{Mean squared error of different
#'   methods for estimating treatment effects. See Dandl et al. (2021) for details.}
#' }
#' @references
#' Dandl S, Torsten H, Heidi S (2021). Divide or Unite? A Comparison of Two
#' Strategies to Random Forest-type Heterogeneous Treatment Effect Estimation.
#'
#' Wager S, and Athey S (2018). "Estimation and Inference of Heterogeneous Treatment Effects using Random Forests". Journal of the American Statistical Association, 113(523).
#'
#' Nie X, Wager S (2020). “Quasi-Oracle Estimation of Heterogeneous Treatment Effects.” Biometrika. ISSN 0006-3444. doi:10.1093/biomet/asaa076. Asaa076, https://academic.oup.com/biomet/advance-article-pdf/doi/10.1093/biomet/asaa076/33788449/asaa076.pdf, URL https://doi.org/10.1093/biomet/asaa076.
#'
"dandl"
