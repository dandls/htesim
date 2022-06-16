#### Plot results #####

source("setup.R")

#---- Load data ----
res.folder = "results"
resall <- readRDS(file.path(res.folder, "results_study.rds"))

#---- Setup by Athey et al. (2019) ----
res <- resall[resall$setup %in% c(1:16),]
HONESTY <- FALSE # set to TRUE to show also the results for honest forests
source("helpers.R")

#---- Plots ----
plot_results(normalA, scA, ylim = c(-.1, 0.66))
plot_results(normalB, scB, ylim = c(-.1, 0.66))

#---- Results table ----
levA <- c("pF_x3.mF_x3.0" = "$\\mu(x_3) + 0 \\cdot W(x_3)$",
  "0.5.0.tF_exp_x1_x2" = "$\\tau(x_1, x_2) W$",
  "pF_x1.mF_x1.tF_exp_x1_x2" = "$\\mu(x_1) + \\tau(x_1, x_2) W(x_1)$",
  "0.5.mF_x1.tF_exp_x1_x2" = "$\\mu(x_1) + \\tau(x_1, x_2) W$",
  "0.5.mF_x3.tF_exp_x1_x2" = "$\\mu(x_3) + \\tau(x_1, x_2) W$",
  "pF_x3.mF_x3.tF_exp_x1_x2" = "$\\mu(x_3) + \\tau(x_1, x_2) W(x_3)$",
  "pF_x3.0.tF_exp_x1_x2" = "$\\tau(x_1, x_2) W(x_3)$",
  "pF_x4.mF_x3.tF_exp_x1_x2" = "$\\mu(x_3) + \\tau(x_1, x_2) W(x_4)$")


levB <- c("pF_x3.mF_x3.0" = "$\\mu(x_3) + 0 \\cdot (W(x_3) - 0.5)$ ",
  "0.5.0.tF_exp_x1_x2" = "$\\tau(x_1, x_2) (W - 0.5)$",
  "pF_x1.mF_x1.tF_exp_x1_x2" = "$\\mu(x_1) + \\tau(x_1, x_2) (W(x_1) - 0.5)$",
  "0.5.mF_x1.tF_exp_x1_x2" = "$\\mu(x_1) + \\tau(x_1, x_2) (W - 0.5)$",
  "0.5.mF_x3.tF_exp_x1_x2" = "$\\mu(x_3) + \\tau(x_1, x_2) (W - 0.5)$",
  "pF_x3.mF_x3.tF_exp_x1_x2" = "$\\mu(x_3) + \\tau(x_1, x_2) (W(x_3) - 0.5)$",
  "pF_x3.0.tF_exp_x1_x2" = "$\\tau(x_1, x_2) (W(x_3) - 0.5)$",
  "pF_x4.mF_x3.tF_exp_x1_x2" = "$\\mu(x_3) + \\tau(x_1, x_2)(W(x_4) - 0.5)$")

lr <- create_table(normalA, normalB, levA, levB)
lr

#---- Setup by Nie and Wager (2020) ----
res <- resall[resall$setup %in% c(17:24),]
HONESTY <- FALSE # set to TRUE to show also the results for honest forests
source("helpers.R")

#---- Plots ----
plot_results(normalA, scA, ylim = c(-.1, 1.61), cexstrip = 1)
plot_results(normalB, scB, ylim = c(-.1, 1.61), cexstrip = 1)

#---- Results table ----
lev <- c("pF_eta_x1_x2.mF_sin_x1_x5.tF_div_x1_x2" = "Setup A",
  "0.5.mF_max_x1_x5.tF_log_x1_x2" = "Setup B",
  "pF_x2_x3.mF_log_x1_x3.1" = "Setup C",
  "pF_exp_x1_x2.mF_max2_x1_x5.tF_max_x1_x5" = "Setup D")

lr <- create_table(normalA, normalB, lev, lev)
lr
