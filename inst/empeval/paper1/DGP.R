#--- Define data generating process ----

set.seed(290875)

create_setups <- function(model) {

  ### set-up dgp objects. Note that tau(x) always depends on X1 and X2
  setups <- vector(mode = "list", length = 0)

  setups <- c(setups,
    ### Setup Nie (2020) - conditional mean depends on treatment effect t
    ### Setup A
    list(htesim::dgp(p = pF_eta_x1_x2, m = mF_sin_x1_x5, t = tF_div_x1_x2, ol = .5, model = model, xmodel = "unif")),
    ### Setup B
    list(htesim::dgp(p = 0.5, m = mF_max_x1_x5, t = tF_log_x1_x2, ol = .5, model = model, xmodel = "normal")),
    ### Setup C
    list(htesim::dgp(p = pF_x2_x3, m = mF_log_x1_x3, t = 1, ol = .5, model = model, xmodel = "normal")),
    ### Setup D
    list(htesim::dgp(p = pF_exp_x1_x2, m = mF_max2_x1_x5, t = tF_max_x1_x5, ol = .5, model = model, xmodel = "normal"))
  )

  return(setups)
}

if (is.null(REPL)) {
  REPL <- 100L
}
DIM <- c(10, 20)
NSIM <- c(800, 1600)

create_args <- function(setups) {
  args <- expand.grid(setup = 1:length(setups),
    ### sample sizes
    nsim = NSIM,
    ### number of covariates
    dim = DIM,
    repl = 1:REPL)
  ### we store the seed and can thus reproduce each dataset whenever needed
  args$seed <- round(runif(nrow(args)) * 1e8)
  args$model <- sapply(setups, function(x) x$model)[args$setup]
  args$p <- sapply(setups, function(x) x$pname)[args$setup]
  args$m <- sapply(setups, function(x) x$mname)[args$setup]
  args$t <- sapply(setups, function(x) x$tname)[args$setup]
  args$sd <- sapply(setups, function(x) x$sdname)[args$setup]
  args$ol <- sapply(setups, function(x) x$ol)[args$setup]
  return(args)
}
