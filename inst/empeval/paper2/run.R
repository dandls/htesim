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
min_update <- 20L

prt <- list(replace = FALSE, fraction = .5)
prt_honest <- list(replace = FALSE, fraction = c(0.25, 0.25))

run <- function(
  ### data (with benefits, such as the ground truth)
  d,
  ### fit propensities only when e != .5
  ### use What to center W 
  propensities = !identical(unique(predict(d, newdata = d)[, "pfct"]),
    .5),
  ### use Yhat to center Y
  marginal_mean = FALSE,
  ### use Weibull models for Weibull DGP by default
  Cox = FALSE,
  ### calculate honest trees
  honesty = FALSE,
  ### see progress bar
  TRACE = TRUE,
  ### return object not MSE
  object = FALSE,
  ### compute What/pihat (if prognostic = TRUE) and Yhat (if marginal_mean = TRUE) based on paper of Gao and Hastie (2021)
  GAO = TRUE,
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
  
  ### standard choice in grf
  mtry <- min(floor(sqrt(ncol(d) - 2) + 20), ncol(d) - 2)
  
  ### set up splitting proportions
  if (honesty) {
    perturb <- prt_honest
  } else {
    perturb <- prt
  }
  
  ### initialize
  offset <- NULL
  d$offset <- 0
  y.hat <- NULL
  nm <- nmt <- "trt1"
  
  ### ground truth
  testxdf <- attributes(d)$testxdf
  tau <- predict(d, newdata = testxdf)

  ### replace treatment indicator by ground truth in randomized trial
  W.hat <- if (identical(unique(predict(d, newdata = testxdf)[, "pfct"]), .5)) .5 else NULL
  trt <- d$trt
  
  ## Compute myWhat/a
  if (propensities) {
    ### call causal_forest to compute W.hat and Y.hat
    # needs matrix input
    X <- as.matrix(d[, grep("^X", colnames(d))])
    W <- (0:1)[d$trt]
    if (attributes(d)$truth$mod != "normal") {
      Y <- rep(0, nrow(X))
    } else {
      Y <- as.numeric(d$y)
    }
    cf <- causal_forest(X = X,
      Y = Y, W = W, W.hat = W.hat,
      stabilize.splits = TRUE,
      min.node.size = min_size_group, sample.fraction = prt$fraction,
      mtry = mtry, ci.group.size = 1,
      num.trees = NumTrees, honesty = honesty)
    myW.hat <- cf$W.hat
    a <- myW.hat
    # update a if strategy by Gao and Hastie
    # even if strategy by Gao and Hastie is not applied,
    # similar calculations might be necessary for centering Y (marginal_mean = TRUE)
    if (GAO || marginal_mean) {
      d0 <- d[d$trt == 0,]
      d0$trt <- NULL
      d1 <- d[d$trt == 1,]
      d1$trt <- NULL
      if (attributes(d)$truth$mod == "binomial") {
        d0$y <- c(FALSE, TRUE)[d0$y]
        d1$y <- c(FALSE, TRUE)[d1$y]
        boost.Y0 <- gbm(y ~ ., data = d0, distribution = "bernoulli", interaction.depth = 2L)
        boost.Y1 <- gbm(y ~ ., data = d1, distribution = "bernoulli", interaction.depth = 2L) 
        if (GAO) {
          p0 <- predict(boost.Y0, newdata = d, type = "response")
          p1 <- predict(boost.Y1, newdata = d, type = "response")
          p01 <- p0 * (1 - p0) / (p1 * (1 - p1))
          a <- myW.hat / (myW.hat + (1 - myW.hat) * p01)
        }
      } else if (attributes(d)$truth$mod == "weibull") {
        if (Cox) {
          if (GAO) {
            # get censoring dataset
            dcens0 <- d0
            dcens0$cens <- as.logical(dcens0$y[,2])
            dcens0$y <- NULL
            dcens1 <- d1
            dcens1$cens <- as.logical(dcens1$y[,2])
            dcens1$y <- NULL
            # boosting for censoring probability
            boost.cens0 <- gbm(cens ~ ., data = dcens0, distribution = "bernoulli", interaction.depth = 2L)
            boost.cens1 <- gbm(cens ~ ., data = dcens1, distribution = "bernoulli", interaction.depth = 2L)
            prob.cens0 <- predict(boost.cens0, newdata = d, type = "response")
            prob.cens1 <- predict(boost.cens1, newdata = d, type = "response")
            a <- (myW.hat*prob.cens1)/(myW.hat*prob.cens1+(1-myW.hat)*prob.cens0)
          }
        } 
      } 
    } 
    # center treatment indicator
    d$trt <- (0:1)[d$trt] - a
    ### this is the name of the parameter we care for
    nm <- nmt <- "trt"
  }
  
  ## Compute Yhat
  if (marginal_mean) {
    if (attributes(d)$truth$mod == "normal") {
      ### use Y.hat of causal forest to center y
      Y.hat <- cf$Y.hat
      d$y <- d$y - Y.hat
    } else {
      if (attributes(d)$truth$mod == "binomial") {
        lp0 <- as.numeric(predict(boost.Y0, newdata = d, type = "link"))
        lp1 <- as.numeric(predict(boost.Y1, newdata = d, type = "link"))
        if (GAO) {
          # use a for weighting two linear predictors
          nu <- a * lp1 + (1-a) * lp0
        } else {
          # use propensity score for weighting two linear predictors =
          # "naive" Robinson
          nu <- myW.hat * lp1 + (1 - myW.hat) * lp0
        }
      } else if (attributes(d)$truth$mod == "weibull") {
          boost.Y0 <- gbm(y ~., data = d0, distribution = "coxph", interaction.depth = 2L)
          boost.Y1 <- gbm(y ~., data = d1, distribution = "coxph", interaction.depth = 2L)
          lp0 <- predict(boost.Y0, newdata = d, type = "link")
          lp1 <- predict(boost.Y1, newdata = d, type = "link")
          if (Cox & GAO) {
            # use a for weighting two linear predictors
            a0 <- ((1-myW.hat)*prob.cens0)/(myW.hat*prob.cens1+(1-myW.hat)*prob.cens0)
            nu <- a * lp1 + a0 * lp0
          } else {
            # use propensity score for weighting two linear predictors
            # "naive" Robinson
            nu <- myW.hat * lp1 + (1 - myW.hat) * lp0
          }
      } else if (attributes(d)$truth$mod == "polr") {
        boost.Y0 <- blackboost(y ~ ., data = d0, family = PropOdds())
        boost.Y1 <- blackboost(y ~ ., data = d1, family = PropOdds())
        lp0 <- as.numeric(predict(boost.Y0, newdata = d, type = "link"))
        lp1 <- as.numeric(predict(boost.Y1, newdata = d, type = "link"))
        # use propensity score for weighting two linear predictors
        # "naive" Robinson
        nu <- myW.hat * lp1 + (1 - myW.hat) * lp0
      }
      # define nu as offset
      d$offset <- offset <- nu
    }
  }
  

  ### Setup up models according to the underlying model
  ### normally distributed response
  if (attributes(d)$truth$mod == "normal") {
    d$offset <- NULL
      ### set-up linear model
        m <- lm(y ~ trt, data = d)
      } else if (attributes(d)$truth$mod == "weibull") {
      if (Cox) {
        ### set-up Cox model
          m <- coxph(y ~ trt + offset(offset), data = d)
      } else {
        ### set-up Weibull model
        m <- as.mlt(Survreg(y ~ trt, data = d, offset = offset))
      }
    }
  else if (attributes(d)$truth$mod == "polr") {
    ### set-up ordinal regression model
      m <- polr(y ~ trt + offset(offset), data = d)
  } else {
        m <- glm(y ~ trt + offset(offset), data = d, family = binomial)
  }
  
  ### fit forest and partition wrt to BOTH intercept and treatment effect
  if (attributes(d)$truth$mod == "weibull" & !Cox) {
    nm <- c("(Intercept)", nm)
    # remove offset from dataframe to not use it as a splitting variable
    # for equalized == TRUE, the offset is still available in vector offset
    d$offset <- NULL
    rf <- trtf::traforest(m, formula = y | trt ~ ., data = d, ntree = NumTrees,
      perturb = perturb, mtry = mtry, control = ctrl, parm = nm, min_update = min_update,
      mltargs = list(offset = offset), trace = TRACE)
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
    return(mse)
  } else {
    cf <- model4you::pmodel(rf, newdata = testxdf)
  }
  mod <- attributes(d)$truth$mod
  if (mod == "weibull" & Cox) {
    ret <- -c(cf)
  } else if (mod == "polr") {
    ret <- c(cf)
  } else if (mod %in% c("weibull", "binomial", "normal")) {
    ret <- cf[, nmt]
  }
  mse <- mean((tau[, "tfct"] - ret)^2)
  return(mse)
}

# Specify functions for batchtools
# No centering, original mob
fun.mob <- function(instance, ...) {
  mob <- run(instance, propensities = FALSE, 
    honesty = FALSE, ...)
  return(mob)
}

fun.mobhonest <- function(instance, ...) {
  mob <- run(instance, propensities = FALSE,
     honesty = TRUE,  ...)
  return(mob)
}

fun.hybridgao <- function(instance,  ...) {
  mob <- run(instance, propensities = TRUE, 
    honesty = FALSE, GAO = TRUE, ...)
  return(mob)
}

# centering W based on Gao and Hastie (2021)
fun.hybridgaohonest <- function(instance, ...) {
  mob <- run(instance,  propensities = TRUE, 
    honesty = TRUE, GAO = TRUE, ...)
  return(mob)
}

# centering W based on Robinson
fun.hybrid <- function(instance, ...) {
  mob <- run(instance, propensities = TRUE, 
    honesty = FALSE, GAO = FALSE, ...)
  return(mob)
}

fun.hybridhonest <- function(instance, ...) {
  mob <- run(instance,  propensities = TRUE, 
    honesty = TRUE, GAO = FALSE, ...)
  return(mob)
}

# Centering Y and W based on Gao and Hastie (2021)
fun.equalizedgao <- function(instance, ...) {
  mob <- run(instance, propensities = TRUE, marginal_mean = TRUE,
     honesty = FALSE, GAO = TRUE, ...)
  return(mob)
}

fun.equalizedgaohonest <- function(instance, ...) {
  mob <- run(instance, propensities = TRUE, marginal_mean = TRUE,
     honesty = TRUE, GAO = TRUE,  ...)
  return(mob)
}

# Centering Y and W based on Robinson
fun.equalized <- function(instance, ...) {
  mob <- run(instance, propensities = TRUE, marginal_mean = TRUE,
     honesty = FALSE, GAO = FALSE, ...)
  return(mob)
}

fun.equalizedhonest <- function(instance, ...) {
  mob <- run(instance, propensities = TRUE, marginal_mean = TRUE,
     honesty = TRUE, GAO = FALSE, ...)
  return(mob)
}

# Special treatments for cox model!
fun.mobcox <- function(instance, ...) {
  run(instance, propensities = FALSE,  Cox = TRUE, 
    GAO = FALSE, ...)
}

fun.mobcoxhonest <- function(instance, ...) {
  run(instance, propensities = FALSE,  Cox = TRUE,
    honesty = TRUE, GAO = FALSE, ...)
}

fun.hybridcox <- function(instance, ...) {
  run(instance, propensities = TRUE,  Cox = TRUE,
    honesty = FALSE, GAO = FALSE, ...)
}

fun.hybridcoxhonest <- function(instance, ...) {
  run(instance, propensities = TRUE,  Cox = TRUE,
    honesty = TRUE, GAO = FALSE, ...)
}

fun.equalizedcox <- function(instance, ...) {
  run(instance, propensities = TRUE, marginal_mean = TRUE,  Cox = TRUE,
    honesty = FALSE, GAO = FALSE, ...)
}

fun.equalizedcoxhonest <- function(instance, ...) {
  run(instance, propensities = TRUE, marginal_mean = TRUE,  Cox = TRUE,
    honesty = TRUE, GAO = FALSE, ...)
}

# Centering W based on Gao and Hastie (2021)
fun.hybridgaocox <- function(instance, ...) {
  run(instance, propensities = TRUE,  Cox = TRUE,
    honesty = FALSE, GAO = TRUE, ...)
}

# Centering Y and W based on Gao and Hastie (2021)
fun.equalizedgaocox <- function(instance,  ...) {
  mob <- run(instance, propensities = TRUE, marginal_mean = TRUE,  
    Cox = TRUE, honesty = FALSE, GAO = TRUE, ...)
  return(mob)
}