rm(list = ls())
#----
# 0) Load helper functions & libraries
#----

# Setup
source("def.R")
paste(NumTrees, CORES,  REPL, sep = ", ")

# Load packages and helper functions
source("../libs.R")
source("DGP.R")
source("run.R")

#-------
# 1) Setup DGP
# Result: lookup table resdf
#-------
setups <- create_setups(model = "normal")

## Code to use older setups (seeds)
if (!is.null(oldresdf)) {
  resdf <- readRDS(oldresdf)
  onemethod <- unique(resdf$algorithm)[1]
  resdf <- resdf[resdf$algorithm == onemethod,
    c("setup", "nsim", "dim", "repl", "seed", "model", "p", "m", "t", "sd", "ol")]
  if (REPL < max(resdf$repl)) {
    resdf <- resdf[resdf$repl %in% 1:REPL,]
  }
} else {
  resdf <- create_args(setups)
}
if (!is.null(StudyIDs)) {
  resdf <- resdf[resdf$setup %in% StudyIDs, ]
}
resdf$rows <- 1:nrow(resdf)
message("nrows of resdf is: ", nrow(resdf))


#-----
# 2) Create study environment (TEST/NO TEST)
# Result: experimental registry
#-----
OVERWRITE <- TRUE

# Create registry
if (file.exists(registry_name)) {
  if (OVERWRITE) {
    unlink(registry_name, recursive = TRUE)
    reg = makeExperimentRegistry(file.dir = registry_name, packages = packages,
      seed = 123L)
  } else {
    reg = loadRegistry(registry_name, writeable = TRUE)
  }
} else {
  reg = makeExperimentRegistry(file.dir = registry_name, packages = packages,
    seed = 123L)
}

# Specifications
reg$default.resources = list(
  ntasks = 1L,
  ncpus = 1L,
  nodes = 1L,
  clusters = "serial")
reg$cluster.functions = makeClusterFunctionsMulticore(CORES)

#----
# 3) Add problem = training dataset based on lookup table resdf
#----
fun = function(job, data, setup, nsim, dim, seed, ...) {
  simulate(setups[[setup]], nsim = nsim,
    dim = dim, nsimtest = 1000L, seed = seed)
}

addProblem("normal", fun = fun, reg = reg)
prob.designs <- list()
prob.designs$normal <- resdf

#-----
# 4) Add algorithms
# Functions are defined in run.R
#-----
algo.designs <- list()
for (method in methods) {
  addAlgorithm(method, fun = eval(parse(text = paste("fun", method, sep = "."))), reg = reg)
  algo.designs[[method]] <- data.frame()
}

#----
# 5) add experiments
#----
addExperiments(prob.designs, algo.designs, reg = reg)

#----
# 6) check setup
#----
# check what has been created
summarizeExperiments(reg = reg)
unwrap(getJobPars(reg = reg))
# testJob(3L)

#-----
# 7) submit jobs
#----
submitJobs()
waitForJobs()

#----
# 8) Save results as rds
#----
res <- ijoin(
  getJobPars(),
  reduceResultsDataTable(fun = function(x) list(res = x))
)
res <- unwrap(res, sep = ".")
names(res) <- stringr::str_remove_all(names(res), "prob.pars.")

if (!is.null(oldresdf)) {
  resold <- readRDS(oldresdf)
  if (names(resold) == names(res)) res <- rbind(res, resold)
}

saveRDS(res, file = resname)





