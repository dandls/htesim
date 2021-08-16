### number of trees in forest (grf default)
NumTrees <- 500
min_size_group <- 7L
min_node_size <- min_size_group*2L
### for splitting, no early stopping, large trees
ctrl <- ctree_control(testtype = "Univ", minsplit = 2,
  minbucket = min_node_size,
  mincriterion = 0, saveinfo = FALSE)
ctrl$converged <- function(mod, data, subset) {
  ### W - W.hat is not a factor
  if (!is.factor(data$trt)) {
    return(all(table(data$trt[subset] > 0) > min_size_group))
  } else {
    ### at least min_size_group obs in both treatment arms
    return(all(table(data$trt[subset]) > min_size_group))
  }
}
## try to update model parameters for nodes with more than
## 20 observations. Reuse parameters in smaller nodes
min_update <- 20

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
  ### return tau as restricted mean survival time (only for survival data)
  RMST = FALSE,
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

  ### make a copy of the data but remove treatment indicator
  # tmp <- d
  # tmp$trt <- NULL
  offset <- NULL
  y.hat <- NULL
  fixed <- NULL
  nm <- nmt <- "trt1"

  ### replace treatment indicator by treatment indicator - propensities
  W.hat <- if (!propensities && causal_forest) .5 else NULL
  trt <- d$trt

  if ((causal_forest && !prognostic_effect) | propensities | marginal_mean) {
    ### call causal_forest to compute W.hat and Y.hat
    ### make sure all forests center in the same way
    if (attributes(d)$truth$mod == "normal") {
      cf <- causal_forest(X = as.matrix(d[, grep("^X", colnames(d))]),
        Y = d$y, W = (0:1)[d$trt], W.hat = W.hat, # W.hat = W.hat if necessary!
        stabilize.splits = stabilize.splits,
        min.node.size = min_size_group, sample.fraction = prt$fraction,
        mtry = mtry, # ci.group.size = 1,
        num.trees = NumTrees, honesty = honesty)
      myW.hat <- cf$W.hat
      myY.hat <- cf$Y.hat
    } else {
      forest.W <- grf::regression_forest(X = as.matrix(d[, grep("^X", colnames(d))]),
        Y = (0:1)[d$trt], num.trees = max(50, NumTrees/4),
        mtry = mtry, honesty = TRUE, min.node.size = 5L,
        sample.fraction = prt$fraction,
        ci.group.size = 1)
      myW.hat <- predict(forest.W)$predictions
    }
  } else if (causal_forest && prognostic_effect) {
    cf <- multi_arm_causal_forest(X = as.matrix(d[, grep("^X", colnames(d))]),
      Y = d$y, W = d$trt, W.hat = W.hat, split.on.intercept = TRUE,
      stabilize.splits = stabilize.splits,
      min.node.size = min_size_group, sample.fraction = prt$fraction,
      mtry = mtry,
      num.trees = NumTrees, honesty = honesty)
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
  if (attributes(d)$truth$mod == "normal") {
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
  } else if (attributes(d)$truth$mod == "weibull") {
    if (causal_forest) { ##SD
      survobj <- as.data.frame(as.matrix(d$y))
      cf <- causal_survival_forest(X = as.matrix(d[, grep("^X", colnames(d))]),
        Y = survobj$time, W = (0:1)[trt], D = survobj$status,
        min.node.size = min_node_size, sample.fraction = perturb$fraction,
        mtry = mtry, # ci.group.size = 1,
        num.trees = NumTrees, honesty = honesty)
      ret <- predict(cf,
        newdat = testxdf[grep("^X", colnames(testxdf))])$predictions
      ### tau is restricted mean survival time difference in this case
      return(mean((tau[, "rmst"] - ret)^2))
    }
    else {
      if (Cox) {
        ### set-up Cox model
        m <- coxph(y ~ trt, data = d)
      } else {
        ### set-up Weibull model
        m <- as.mlt(Survreg(y ~ trt, data = d, offset = offset))
      }
    }
  }
  else if (attributes(d)$truth$mod == "polr") {
    ### set-up ordinal regression model
    m <- polr(y ~ trt, data = d)
  } else {
    m <- glm(y ~ trt, data = d, family = binomial)
  }
  
  if (ATE) {
    if(attributes(d)$truth$mod == "weibull" & !Cox) {
      ret <- m$coef[nmt]
    } else {
      ret <- m$coefficients[nmt]
    } 
    if (Cox) {
      ret <- -ret
    }
    mse <- mean((tau[, "tfct"] - ret)^2)
    return(mse)
  }
  
  ### fit forest and partition wrt to BOTH intercept and treatment effect
  if (attributes(d)$truth$mod == "weibull" & !Cox) {
    rf <- trtf::traforest(m, formula = y | trt ~ ., data = d, ntree = NumTrees,
      perturb = perturb, mtry = mtry, control = ctrl, parm = nm, min_update = min_update,
      mltargs = list(fixed = fixed, offset = offset),
      trace = TRACE)
  } else {
    rf <- model4you::pmforest(m, data = d, ntree = NumTrees, perturb = perturb,
      mtry = mtry, control = ctrl, trace = TRACE)
  }

  rmstdiff <- NA

  if (object) return(rf)

  ### estimate model coefficients on test set
  if (attributes(d)$truth$mod == "weibull" & !Cox) {
    cf <- predict(rf, newdata = testxdf, type = "coef")
    cf <- do.call("rbind", cf)

    ### extract treatment effect (on log-hazard ratio scale)
    ret <- cf[, nmt]
    mse <- mean((tau[, "tfct"] - ret)^2)

    if (RMST) {
      time <- max(d$y[d$y[,2] > 0, 1])
      ### compute restricted mean survival time
      mX <- cf[, "(Intercept)"]
      sc <- cf[, "log(y)"]
      tX <- cf[, nmt]
      tm <- matrix(seq(from = 0, to = time, length.out = 100),
        nrow = length(mX), ncol = 100, byrow = TRUE)
      rmstdiff <- rowSums(exp(-exp(sc * tm - (mX + tX)))) * (time / ncol(tm)) -
        rowSums(exp(-exp(sc * tm - mX))) * (time / ncol(tm))
      rmstmse <- (mean((tau[, "rmst"] - rmstdiff)^2))
      attr(mse, "rmst") <- rmstmse
    }
    return(mse)
  } else {
    cf <- model4you::pmodel(rf, newdata = testxdf)
  }
  mod <- attributes(d)$truth$mod
  if (mod == "weibull" & Cox) {
    ret <- -c(cf)
  } else if (mod == "polr" | (mod == "normal" & !prognostic_effect)) {
    ret <- c(cf)
  } else if (mod %in% c("weibull", "binomial", "normal")) {
    ret <- cf[, nmt]
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
  estfun(object)[,"trt", drop = FALSE]
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

fun.mob <- function(instance, RMST = FALSE, ...) {
  mob <- run(instance, propensities = FALSE, causal_forest = FALSE,
             honesty = FALSE, RMST = RMST, ...)
  return(mob)
}

fun.mobhonest <- function(instance, RMST = FALSE, ...) {
  mob <- run(instance, propensities = FALSE,
     causal_forest = FALSE, honesty = TRUE, RMST = RMST, ...)
  return(mob)
}

fun.hybrid <- function(instance, RMST = FALSE, ...) {
  mob <- run(instance, propensities = TRUE, causal_forest = FALSE,
      honesty = FALSE, RMST = RMST, ...)
  return(mob)
}

fun.hybridhonest <- function(instance, RMST = FALSE, ...) {
  mob <- run(instance,  propensities = TRUE, causal_forest = FALSE,
      honesty = TRUE, RMST = RMST, ...)
  return(mob)
}

fun.equalized <- function(instance, RMST = FALSE, ...) {
  mob <- run(instance, propensities = TRUE, marginal_mean = TRUE,
    causal_forest = FALSE, honesty = FALSE, RMST = RMST, ...)
  return(mob)
}


fun.equalizedhonest <- function(instance, RMST = FALSE, ...) {
  mob <- run(instance, propensities = TRUE, marginal_mean = TRUE,
    causal_forest = FALSE, honesty = TRUE, RMST = RMST, ...)
  return(mob)
}

fun.mobcf <- function(instance, RMST = FALSE, ...) {
    mob <- run(instance, propensities = TRUE, marginal_mean = TRUE,
      prognostic_effect = FALSE, causal_forest = FALSE, honesty = FALSE,
      RMST = RMST, ...)
    return(mob)
  }


fun.mobcfhonest <- function(instance, RMST = FALSE, ...) {
  mob <- run(instance, propensities = TRUE, marginal_mean = TRUE,
    prognostic_effect = FALSE, causal_forest = FALSE, honesty = TRUE,
    RMST = RMST, ...)
  return(mob)
}

# Special treatments for cox model!
fun.mobcox <- function(instance, ...) {
  run(instance, propensities = FALSE, causal_forest = FALSE, Cox = TRUE, ...)
}

fun.mobcoxhonest <- function(instance, ...) {
  run(instance, propensities = FALSE, causal_forest = FALSE, Cox = TRUE,
     honesty = TRUE, ...)
}

fun.hybridcox <- function(instance, ...) {
  run(instance, propensities = TRUE, causal_forest = FALSE, Cox = TRUE,
     honesty = FALSE, ...)
}

fun.hybridcoxhonest <- function(instance, ...) {
  run(instance, propensities = TRUE, causal_forest = FALSE, Cox = TRUE,
     honesty = TRUE, ...)
}


# Average Treatment Effect Methods: get coef from base model
# original treatment indicator
fun.bm <- function(instance, RMST = FALSE, ...) {
 run(instance, propensities = FALSE, causal_forest = FALSE,
    honesty = FALSE, RMST = RMST, ATE = TRUE, ...)
}

fun.bmcox <- function(instance, ...) {
  run(instance, propensities = FALSE, causal_forest = FALSE, Cox = TRUE,
    ATE = TRUE, ...)
}

# orthogonalized treatment indicator
fun.bmhybrid <-  function(instance, RMST = FALSE, ...) {
  run(instance, propensities = TRUE, causal_forest = FALSE,
  honesty = FALSE, RMST = RMST, ATE = TRUE, ...)
}

fun.bmhybridcox <- function(instance, ...) {
  run(instance, propensities = TRUE, causal_forest = FALSE, Cox = TRUE,
    honesty = FALSE, ATE = TRUE, ...)
}
