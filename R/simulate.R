#' Simulate from a given dgp object
#' @param object (list/gpd) Manual for data generation.
#' @param nsim (numeric(1)) Number of observations (n), default 1.
#' @param seed (numeric(1)) Seed for data generation, default NULL.
#' @param dim (numeric(1)) Dimension, number of predictors (p), default 4.
#' @param nsimtest (numeric(1)) Number of observations (n) for test dataset, default 1000.
#' @return data.frame of class simdpg with columns: x, y and trt/w.
#' @export
simulate.dgp <- function(object, nsim = 1, seed = NULL, dim = 4,
  nsimtest = 1000) {

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
    testx <- matrix(rnorm(nsimtest * dim), nrow = nsimtest, ncol = dim)
  } else if (object$xmodel == "unif") {
    x <- matrix(runif(nsim * dim), nrow = nsim, ncol = dim)
    testx <- matrix(runif(nsimtest * dim), nrow = nsimtest, ncol = dim)
  }
  colnames(x) <- colnames(testx) <- paste("X", 1:ncol(x), sep = "")
  testxdf <- as.data.frame(testx)
  testxdf$trt <- factor(c(0, 1))[2]

  tX <- object$tfct(x)
  pX <- object$pfct(x)
  mX <- object$mfct(x)
  sd <- object$sdfct(x)
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
    y <- Surv(ifelse(y < cens, y, cens), ifelse(y < cens, 1, 0))
  }
  df <- data.frame(x, y = y, trt = factor(trt))
  attributes(df)$truth <- object
  ### test data set (predictor variables only)
  attributes(df)$testxdf <- testxdf
  ### seed for model fitting
  attributes(df)$runseed <- round(runif(1) * 100000)
  class(df) <- c("simdgp", class(df))
  return(df)
}
