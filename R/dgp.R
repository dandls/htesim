#' Create manual for data generation
#' @param p (numeric(1), character(1), function) Propensities P(trt = 1|x) (pi(x)),
#' default = 0.5 means that there is no confounding.
#' @param m (numeric(1), character (1), function) Prognostic effect (mu(x)), default 0 means that
#' there is no prognostic effect.
#' @param t (numeric(1), character(1), function) Predictive treatment effect (tau(x)),
#' default 0 means that there is no treatment effect.
#' @param sd (numeric(1), character(1), function) Standard deviation of normal distribution, only
#' has an effect if model = "normal".
#' @param pi (numeric(1)) how much of predictive effect is added to prognostic effect.
#' The default is pi = 0 which means that the conditional mean does not depend on treatment effect.
#' @param model ("normal"|"weibull"|"binomial"|"polr") Name of used model to simulate outcome y.
#' @param xmodel ("normal") Name of used model to simulate covariates x.
#' @return list of class gpd with entries:
#' pfct (pi(x)), mfct (mu(x)), tfct (tau(x)), sdfct (sd for normal model),
#' model (model used)).
#' @seealso \code{\link{pF}}, \code{\link{mF}} and \code{\link{tF}}
#' for predefined functions for \code{p}, \code{m} and \code{t}.
#' @references
#' Wager S, and Athey S (2018). "Estimation and Inference of Heterogeneous Treatment Effects using Random Forests". Journal of the American Statistical Association, 113(523).
#'
#' Nie X, Wager S (2020). “Quasi-Oracle Estimation of Heterogeneous Treatment Ef- fects.” Biometrika. ISSN 0006-3444. doi:10.1093/biomet/asaa076. Asaa076, https://academic.oup.com/biomet/advance-article-pdf/doi/10.1093/biomet/asaa076/33788449/asaa076.pdf, URL https://doi.org/10.1093/biomet/asaa076.
#'
#' @examples
#' # Wager and Athey (2018) - first experiment
#' dgp1 <- dgp(p = pF_x1, m = mF_x1, t = 0, model = "normal", xmodel = "unif")
#'
#' # Wager and Athey (2018) - second experiment
#' dgp2 <- dgp(p = 0.5, m = 0, t = tF_exp_x1_x2, model = "normal", xmodel = "unif")
#'
#' # Wager and Athey (2018) - third experiment
#' dgp3 <- dgp(p = 0.5, m = 0, t = tF_exp2_x1_x2, model = "normal", xmodel = "unif")
#'
#' # Nie and Wager (2020) - Setup A
#' dgpA <- dgp(p = pF_eta_x1_x2, m = mF_sin_x1_x5, t = tF_div_x1_x2, model = "normal", xmodel = "unif")
#'
#' # Nie and Wager (2020) - Setup B
#' dgpB <- dgp(p = 0.5, m = mF_max_x1_x5, t = tF_log_x1_x2, model = "normal", xmodel = "normal")
#'
#' # Nie and Wager (2020) - Setup C
#' dgpC <- dgp(p = pF_x2_x3, m = mF_log_x1_x3, t = 1, model = "normal", xmodel = "normal")
#'
#' # Nie and Wager (2020) - Setup D
#' dgpD <- dgp(p = pF_exp_x1_x2, m = mF_max2_x1_x5, t = tF_max_x1_x5, model = "normal", xmodel = "normal")
#'
#' @export
dgp <- function(p = 0.5, m = 0, t = 0, sd = 1, pi = 0, model = c("normal", "weibull", "binomial", "polr"), xmodel = c("normal", "unif")) {

  cl <- match.call()

  # sanity checks
  assert_true(is.function(p) || is.numeric(p) || is.character(p))
  assert_true(is.function(m) || is.numeric(m) || is.character(m))
  assert_true(is.function(t) || is.numeric(t) || is.character(t))
  assert_number(sd, lower = 0)
  assert_number(pi)

  model <- tryCatch({match.arg(model)},
    error = function(e) {
      stop("Assertion on 'model' failed: Must be element of set {'normal', 'weibull', 'binomial', 'polr'}")
    })
  xmodel <- tryCatch({match.arg(xmodel)},
    error = function(e) {
      stop("Assertion on 'xmodel' failed: Must be element of set {'normal', 'unif'}")
    })

  pfct <- sanitize_fct(p, cl$p)
  mfct <- sanitize_fct(m, cl$m)
  tfct <- sanitize_fct(t, cl$t)
  sdfct <- sanitize_fct(sd, cl$sd)

  ### mean effect
  mf <- function(x) mfct(x) - pi * tfct(x)

  ret <- list(pfct = pfct, mfct = mf, tfct = tfct, sdfct = sdfct, pi = pi,
    model = model, xmodel = xmodel, pname = deparse(substitute(p)),
    mname = deparse(substitute(m)), tname = deparse(substitute(t)),
    sdname = deparse(substitute(sd)))
  class(ret) <- "dgp"
  return(ret)
}

# Helper function to properly handle propensity/prognostic/treatment effects as
# numerics, characters or functions
sanitize_fct <- function(fct, call) {

  # proper handling of characters
  # stop if function not found
  if (is.character(fct)) {
    tryCatch({fct <- eval(parse(text = fct))},
      error = function(e) {
        stop(paste("function ", call, " not found", sep = "'"))
      })
  }

  # proper handling of numerics
  if (is.numeric(fct)) {
    fct <- return(function(x) return(rep(fct, NROW(x))))
  }

  return(fct)
}


#' @title Simulate from a given dgp object
#' @description Simulate data from a given \code{dgp} object.
#' @param object (list/gpd) Manual for data generation.
#' @param nsim (numeric(1)) Number of observations (n), default 1.
#' @param seed (numeric(1)) Seed for data generation, default NULL.
#' @param dim (numeric(1)) Dimension, number of predictors (p), default 4.
#' @param nsimtest (numeric(1)) Number of observations (n) for test dataset.
#' Default NULL means that no test dataset is generated.
#' @seealso \code{\link{dgp}}
#' @return data.frame of class simdpg with columns: x, y and trt and attributes
#'  * \code{truth}: the  \code{object}
#'  * \code{testxdf}: the test data with \code{nsimtest} rows and \code{dim} columns
#'  * \code{runseed}: an additional seed (e.g., to run a method)
#' @examples
#'
#' # Wager and Athey (2018) - first experiment
#' dgp1 <- dgp(p = pF_x1, m = mF_x1, t = 0, model = "normal", xmodel = "unif")
#' sim1 <- simulate(dgp1, nsim = 500L, d = 2L, nsimtest = 1000L) # in paper d = {2, 5, 10, 15, 20, 30}
#' head(sim1)
#'
#' # Nie and Wager (2020) - Setup A
#' dgpA <- dgp(p = pF_eta_x1_x2, m = mF_sin_x1_x5, t = tF_div_x1_x2,
#'  model = "normal", xmodel = "unif")
#' simA <- simulate(dgpA, nsim = 500L, d = 6L, nsimtest = 500L)
#' head(simA)
#'
#' # Get test data
#' testdf <- attr(simA, "testxdf")
#' head(testdf)
#'
#' @export
simulate.dgp <- function(object, nsim = 1, dim = 4, nsimtest = NULL, seed = NULL) {

  ###  input checks
  assertIntegerish(nsim, lower = 1, len = 1, any.missing = FALSE)
  assertIntegerish(dim, lower = 1, len = 1, any.missing = FALSE)
  assertIntegerish(nsimtest, lower = 1, len = 1, any.missing = FALSE, null.ok = TRUE)
  assertNumber(seed, null.ok = TRUE)

  ### mlt & tram packages required for Weibull models
  if (object$model == "weibull") {
    if (!requireNamespace("mlt", quietly = TRUE)) {
      stop("Package \"mlt\" needed for this function to work. Please install it.",
        call. = FALSE)
    }
    if (!requireNamespace("tram", quietly = TRUE)) {
      stop("Package \"tram\" needed for this function to work. Please install it.",
        call. = FALSE)
    }

    if (!requireNamespace("survival", quietly = TRUE)) {
      stop("Package \"survival\" needed for this function to work. Please install it.",
        call. = FALSE)
    }
  }

  ### set and re-set seed
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  ### prognostic and predictive variables
  ### generate dedicated test set here
  if (object$xmodel == "normal") {
    x <- matrix(rnorm(nsim * dim), nrow = nsim, ncol = dim)
    if (!is.null(nsimtest)) {
      testx <- matrix(rnorm(nsimtest * dim), nrow = nsimtest, ncol = dim)
    }
  } else if (object$xmodel == "unif") {
    x <- matrix(runif(nsim * dim), nrow = nsim, ncol = dim)
    if (!is.null(nsimtest)) {
      testx <- matrix(runif(nsimtest * dim), nrow = nsimtest, ncol = dim)
    }
  }
  colnames(x) <- paste("X", 1:ncol(x), sep = "")

  tX <- tryCatch({object$tfct(x)},
    error = function(e) {
      stop("treatment effects function operates on variables out of bounds, increase dim")
    })
  pX <- tryCatch({object$pfct(x)},
    error = function(e) {
      stop("propensity score function operates on variables out of bounds, increase dim")
    })

  mX <- tryCatch({object$mfct(x)},
    error = function(e) {
      stop("prognostic function operates on variables out of bounds, increase dim")
    })
  sd <- tryCatch({object$sdfct(x)},
    error = function(e) {
      stop("standard deviation function operates on variables out of bounds, increase dim")
    })
  model <- object$model

  ### sample treatments Wi|Xi ~Bernoulli (e(Xi))
  trt <- rbinom(nsim, size = 1, prob = pX)

  ### effect function
  ### NOTE: this is mfct(x) - .5 * tfct(x) + trt * tfct(x) for pi = .5
  ### and thus the mean always depends on BOTH mX and tX
  efct <- mX + trt * tX

  ### same from models
  y <- switch(model,
    "normal"  = rnorm(nsim, mean = efct, sd = sd),
    ### U = 1 - exp(-exp(2 * log(time) - (mX + (trt - 0.5) * tX)))
    ### NOTE: Survreg has negative = TRUE
    "weibull" = exp((log(-log(1 - runif(nsim))) + efct) / 2),
    ### NOTE: Polr has negative = TRUE
    "polr" = cut(qlogis(runif(nsim)) + efct,
      breaks = c(-Inf, qlogis(1:3/4), Inf), ordered = TRUE),
    "binomial" = factor(rbinom(nsim, size = 1, prob = plogis(efct)),
      labels = 0:1, levels = 0:1)
  )

  if (model == "weibull") {
    ### sample censoring times from the _marginal_ distribution
    ### => Prob(cens < y) ~ .5
    cens <- simulate(mlt::as.mlt(tram::Coxph(y ~ 1, data = data.frame(y = y),
      log_first = TRUE)))
    cens <- trtf:::.R2vec(cens)
    y <- survival::Surv(ifelse(y < cens, y, cens), ifelse(y < cens, 1, 0))
  }
  df <- data.frame(x, y = y, trt = factor(trt))
  attributes(df)$truth <- object
  ### test data set (predictor variables only)
  if (!is.null(nsimtest)) {
    colnames(testx) <- paste("X", 1:ncol(x), sep = "")
    testxdf <- as.data.frame(testx)
    testxdf$trt <- factor(c(0, 1))[2]
    attributes(df)$testxdf <- testxdf
  } else {
    attributes(df)$testxdf <- NULL
  }
  ### seed for model fitting
  attributes(df)$runseed <- round(runif(1) * 100000)
  class(df) <- c("simdgp", class(df))
  return(df)
}


#' @nane predict
#' @rdname predict
#' @title Predict ground truth effects for new data
#' @param object (dgp/simdgp) Manual for data generation (\code{dgp}) or generated dataset (\code{simdgp}).
#' @param newdata (data.frame) New data to predict on.
#' @return data.frame with true values of
#' pfct (pi(x)), mfct (mu(x)), tfct (tau(x)) and sdfct (sd for normal model).
#' @examples
#' # Wager and Athey (2018) - first experiment
#' dgp1 <- dgp(p = pF_x1, m = mF_x1, t = 0, model = "normal", xmodel = "unif")
#' sim1 <- simulate(dgp1, nsim = 500L, d = 2, nsimtest = 1000L)
#' newdata <- attr(sim1, "testxdf")
#' pred1 <- predict(dgp1, newdata)
#' head(pred1)
#' pred2 <- predict(sim1, newdata)
#' all.equal(pred1, pred2)
#' @export
predict.dgp <- function(object, newdata) {
  assert_data_frame(newdata)
  atr <- object[sapply(object, is.function)]
  ret <- sapply(atr, function(f) f(newdata))
  return(data.frame(ret))
}

#' @rdname predict
#' @export
predict.simdgp <- function(object, newdata) {
  objectdgp <- attributes(object)$truth
  predict.dgp(objectdgp, newdata)
}

