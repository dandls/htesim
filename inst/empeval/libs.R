### Requirements

#---- packages

packages = c("survival", "tram", "MASS", "trtf", "partykit", "sandwich",
  "ggplot2", "batchtools", "parallelMap", "stringr", "checkmate")
new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(packages, require, character.only = TRUE)

#--- RForge repos
if ("model4you" %in% installed.packages() & packageVersion("model4you") >= "0.9.5") {
  library("model4you")
} else {
  install.packages("model4you", repos = "http://R-Forge.R-project.org")
  library("model4you")
}

if ("partykit" %in% installed.packages() & packageVersion("partykit") >= "1.2.14") {
  library("partykit")
} else {
  install.packages("partykit", repos = "http://R-Forge.R-project.org")
  library("partykit")
}

if ("grf" %in% installed.packages() & packageVersion("grf") >= "1.2.0") {
  library("grf")
} else {
  devtools::install_github("grf-labs/grf", subdir = "r-package/grf")
  library("grf")
}

devtools::install_github("susanne-207/htesim", auth_token = "56699e45d01135ca797b7e14209eb23d4eeee4b8")
library("htesim")

packages <- c(packages, "htesim", "model4you", "partykit", "grf")

