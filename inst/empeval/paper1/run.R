### number of trees in forest (grf default)
NumTrees <- 500
min_size_group <- 7L
min_node_size <- min_size_group*2L
### for splitting, no early stopping, large trees
ctrl <- ctree_control(testtype = "Univ", minsplit = 2,
  minbucket = 2, splittry = 20, minprob = 0, mincriterion = 0,
  lookahead = TRUE, saveinfo = FALSE)
ctrl$converged = function(mod, data, subset) {
  # The following should be added to avoid nodes with only treated or control samples
  if (length(unique(data$trt[subset])) == 1) return(FALSE)
  ### W - W.hat is not a factor
  if (!is.factor(data$trt)) {
    return(all(table(data$trt[subset] > 0) >= min_size_group))
  } else {
    ### at least min_size_group obs in both treatment arms
    return(all(table(data$trt[subset]) >= min_size_group))
  }
}

prt <- list(replace = FALSE, fraction = .5)
prt_honest <- list(replace = FALSE, fraction = c(0.25, 0.25))

run <- function(
  ### data (with benefits, such as the ground truth)
  d,
  ### fit propensities only when e != .5
  propensities = !identical(unique(predict(d, newdata = d)[, "pfct"]),
    .5),
  marginal_mean = FALSE,
  prognostic_effect = TRUE,
  ### use grf::causal_forest
  causal_forest = FALSE,
  ### use Weibull models for Weibull DGP by default
  Cox = FALSE,
  ### calculate honest trees
  honesty = FALSE,
  ### see progress bar
  TRACE = TRUE,
  ### return object not MSE
  object = FALSE,
  ### stabilize splits in grf::causal_forests
  stabilize.splits = TRUE,
  ### return ATE (trt effect of base model)
  ATE = FALSE,
  ...
) {

  seed <- attributes(d)$runseed

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

  # type <- match.arg(type)
  if (causal_forest) stopifnot(attributes(d)$truth$mod %in% c("normal", "weibull")) # SD

  ### standard choice in grf
  mtry <- min(floor(sqrt(ncol(d) - 2) + 20), ncol(d) - 2)

  if (!causal_forest & honesty) {
    perturb <- prt_honest
  } else {
    perturb <- prt
  }

  ### update splittry
  ctrl$splittry <- ncol(d[, grep("^X", colnames(d))])

  ### initialize
  offset <- NULL
  y.hat <- NULL
  fixed <- NULL
  nm <- nmt <- "trt1"

  ### replace treatment indicator by in randomized trial
  W.hat <- if (identical(unique(predict(d, newdata = d)[, "pfct"]), .5)) .5 else NULL
  trt <- d$trt

  if ((causal_forest && !prognostic_effect) | propensities | marginal_mean) {
    ### call causal_forest to compute W.hat and Y.hat
    ### make sure all forests center in the same way
      cf <- causal_forest(X = as.matrix(d[, grep("^X", colnames(d))]),
        Y = d$y, W = (0:1)[d$trt], W.hat = W.hat, # W.hat = W.hat if necessary!
        stabilize.splits = stabilize.splits,
        min.node.size = min_size_group, sample.fraction = prt$fraction,
        mtry = mtry, ci.group.size = 1,
        num.trees = NumTrees, honesty = honesty, num.threads = 32L)
      myW.hat <- cf$W.hat
      myY.hat <- cf$Y.hat
  } else if (causal_forest && prognostic_effect) {
    cf <- multi_arm_causal_forest(X = as.matrix(d[, grep("^X", colnames(d))]),
      Y = d$y, W = d$trt, W.hat = W.hat, split.on.intercept = TRUE,
      stabilize.splits = stabilize.splits,
      min.node.size = min_size_group, sample.fraction = prt$fraction,
      mtry = mtry,
      num.trees = NumTrees, honesty = honesty, num.threads = 32L)
  }

  if (propensities & !causal_forest) {
    ### use W.hat of causal forest to center trt
    W.hat <- myW.hat
    d$trt <- (0:1)[d$trt] - W.hat
    ### this is the name of the parameter we care for
    nm <- nmt <- "trt"
  }

  ### estimate marginal mean E(Y|X = x)
  if (marginal_mean & !causal_forest) {
    ### use Y.hat of causal forest to center y
    Y.hat <- myY.hat
    d$y <- d$y - Y.hat
  }

  ### ground truth
  testxdf <- attributes(d)$testxdf
  tau <- predict(d, newdata = testxdf)

  ### Setup up models according to the underlying model
  ### normally distributed response
    if (causal_forest) {
      ret <- predict(cf,
        newdat = testxdf[grep("^X", colnames(testxdf))])$predictions
      return(mean((tau[, "tfct"] - ret)^2))
    } else {
      ### set-up linear model
      if (!prognostic_effect) {
        if (!propensities) {
          d$trt <- (0:1)[d$trt]
        }
        m <- lm(y ~ trt, data = d)
        class(m) <- c("nomu", class(m))
      } else {
        m <- lm(y ~ trt, data = d)
      }
    }


  ### fit forest and partition wrt to BOTH intercept and treatment effect
    rf <- model4you::pmforest(m, data = d, ntree = NumTrees, perturb = perturb,
      mtry = mtry, control = ctrl, trace = TRACE)

  if (object) return(rf)

  ### estimate model coefficients on test set

  cf <- model4you::pmodel(rf, newdata = testxdf)
  mod <- attributes(d)$truth$mod
  if (prognostic_effect) {
    ret <- cf[, nmt]
  } else {
    ret <- c(cf)
  }
  mse <- mean((tau[, "tfct"] - ret)^2)
  return(mse)
}

coef.nomu <- function(object, ...) {
  class(object) <- class(object)[-1L]
  coef(object)["trt"]
}

estfun.nomu <- function(object, ...) {
  class(object) <- class(object)[-1L]
  ef <- tryCatch({estfun(object)[,"trt", drop = FALSE]},
    error = function(e) {
      return(rep(0, length(object$model$y)))
    })
  return(ef)
}

update.nomu <- function(object, ...) {
  class(object) <- class(object)[-1L]
  ret <- update(object, ...)
  class(ret) <- c("nomu", class(ret))
  ret
}


# Specify functions for batchtools

fun.cf <- function(instance, ...) {
  run(instance, causal_forest = TRUE,
    prognostic_effect = FALSE, honesty = FALSE)
}

fun.cfhonest <- function(instance, ...) {
  run(instance, causal_forest = TRUE,
    prognostic_effect = FALSE, honesty = TRUE)
}

fun.cfmob <- function(instance, ...) {
  run(instance, causal_forest = TRUE, propensities = FALSE,
    prognostic_effect = TRUE, honesty = FALSE)
}

fun.cfmobhonest <- function(instance, ...) {
  run(instance, causal_forest = TRUE, propensities = FALSE,
    prognostic_effect = TRUE, honesty = TRUE)
}

fun.mob <- function(instance, ...) {
  mob <- run(instance, propensities = FALSE, causal_forest = FALSE,
    honesty = FALSE, ...)
  return(mob)
}

fun.mobhonest <- function(instance, ...) {
  mob <- run(instance, propensities = FALSE,
    causal_forest = FALSE, honesty = TRUE, ...)
  return(mob)
}

fun.hybrid <- function(instance, ...) {
  mob <- run(instance, propensities = TRUE, causal_forest = FALSE,
    honesty = FALSE, ...)
  return(mob)
}

fun.hybridhonest <- function(instance, ...) {
  mob <- run(instance,  propensities = TRUE, causal_forest = FALSE,
      honesty = TRUE, ...)
  return(mob)
}

fun.equalized <- function(instance, ...) {
  mob <- run(instance, propensities = TRUE, marginal_mean = TRUE,
    causal_forest = FALSE, honesty = FALSE, ...)
  return(mob)
}


fun.equalizedhonest <- function(instance, ...) {
  mob <- run(instance, propensities = TRUE, marginal_mean = TRUE,
    causal_forest = FALSE, honesty = TRUE, ...)
  return(mob)
}

fun.mobcf <- function(instance, ...) {
    mob <- run(instance, propensities = TRUE, marginal_mean = TRUE,
      prognostic_effect = FALSE, causal_forest = FALSE, honesty = FALSE,
     ...)
    return(mob)
  }


fun.mobcfhonest <- function(instance, ...) {
  mob <- run(instance, propensities = TRUE, marginal_mean = TRUE,
    prognostic_effect = FALSE, causal_forest = FALSE, honesty = TRUE,
   ...)
  return(mob)
}

# Average Treatment Effect Methods: get coef from base model
# original treatment indicator
fun.bm <- function(instance, ...) {
 run(instance, propensities = FALSE, causal_forest = FALSE,
    honesty = FALSE, ATE = TRUE, ...)
}

# orthogonalized treatment indicator
fun.bmhybrid <-  function(instance, ...) {
  run(instance, propensities = TRUE, causal_forest = FALSE,
  honesty = FALSE, ATE = TRUE, ...)
}
