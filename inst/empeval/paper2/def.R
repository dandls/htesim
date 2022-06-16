## directory to save results
dir.create("results", showWarnings = FALSE)
resname <- "results/results_study.rds"

## parameters for batchtools
CORES <- 30L
registry_name <- "reg_study"

## forest settings
# number of trees
NUMTREES <- 500L

## customization of study
# number of repetitions of each study setting x N x P
REPL <- 100L
# only use subset of problems/DGP? --> remove/add model names from/to following vector
dgps <- c("normal", "binomial", "polr", "weibull")
# only use subset of methods? --> remove/add method names from/to following vector
methods <- c("mob", "hybrid", "equalized", "hybridgao", "equalizedgao", "mobcox", 
  "hybridcox", "equalizedcox", "hybridgaocox", "equalizedgaocox") 
# specific binomial methods
binomialmethods <- methods[grep("gao$", methods)]
# specify cox methods 
coxmethods <- methods[grep("cox$", methods)]
# use older study setting? --> provide path to study results of previous study
oldresdf <- NULL # e.g. oldresdf <- "results/results_study.rds"
# only run subset of studies? --> provide study IDs (based on DGP.R)
StudyIDs <- NULL  # e.g. c(1, 2)

