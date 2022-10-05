#### Plot results #####

source("setup.R")

#---- Load data ----
res.folder = "results"
resall <- readRDS(file.path(res.folder, "results_study.rds"))

#---- Setup by Nie and Wager (2020) ----
res <- resall[resall$setup %in% 1:4,]
HONESTY <- FALSE # set to TRUE to show also the results for honest forests
source("helpers.R")

#---- Plots ----
plot_results(normalB, scB, ylim = c(-.1, 1.61), cexstrip = 1)

#---- Results table ----
lev <- c("pF_eta_x1_x2.mF_sin_x1_x5.tF_div_x1_x2" = "Setup A",
  "0.5.mF_max_x1_x5.tF_log_x1_x2" = "Setup B",
  "pF_x2_x3.mF_log_x1_x3.1" = "Setup C",
  "pF_exp_x1_x2.mF_max2_x1_x5.tF_max_x1_x5" = "Setup D")

create_table(dataA = normalB, dataB = NULL, lev, lev)
