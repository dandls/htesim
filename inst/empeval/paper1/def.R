## directory to save results
dir.create("results", showWarnings = FALSE)
resname <- "results/results_study.rds"

## parameters for batchtools
CORES <- 30L
registry_name <- "reg_normal"

## forest settings
# number of trees
NumTrees <- 500L

## customization of study
# number of repetitions of each study setting x N x P
REPL <- 100L
# only use subset of methods? --> remove/add method names from/to following vector
# e.g. "cf","cfhonest", "mob","mobhonest", "hybrid", "hybridhonest", "equalized", "equalizedhonest", "mobcf", "mobcfhonest"
methods <- c("cf","cfhonest", "mob","mobhonest", "hybrid", "hybridhonest", "equalized", "equalizedhonest", "mobcf", "mobcfhonest")
# use older study setting? --> provide path to study results of previous study
oldresdf <- NULL # e.g. oldresdf <- "results/results_study.rds"
# only run subset of studies? --> provide study IDs (based on DGP.R)
StudyIDs <- NULL  # e.g. c(1, 2)

