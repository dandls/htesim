#### Results #####

source("setup.R")

#---- Load data ----
res.folder = "results"
resall <- readRDS(file.path(res.folder, "results_study.rds"))

#---- Section 4: Study based on Nie and Wager (2020) ----
res <- resall[resall$setup %in% c(22, 24:26),]
res <- res[- grep("gao", res$algorithm), ]
HONESTY <- FALSE
source("helpers.R")

### Plots
plot_results_all(partB, cex = 1)

### Tables
lev <- c("pF_eta_x1_x2.mF_sin_x1_x5.tF_div_x1_x2" = "Setup A",
  "0.5.mF_max_x1_x5.tF_log_x1_x2" = "Setup B",
  "pF_x2_x3.mF_log_x1_x3.1" = "Setup C",
  "pF_exp_x1_x2.mF_max2_x1_x5.tF_max_x1_x5" = "Setup D")

iv_normal <- get_intervals(dataA = NULL, partB, model = "Normal", gao = FALSE, lev, lev)
iv_binomial <- get_intervals(dataA = NULL, partB, model = "Binomial", gao = FALSE, lev, lev)
iv_polr <- get_intervals(dataA = NULL, partB, model = "Multinomial", gao = FALSE, lev, lev)
iv_weibull <- get_intervals(dataA = NULL, partB, model = "Weibull", gao = FALSE, lev, lev)
iv_cox <- get_intervals(dataA = NULL, partB, model = "Cox", gao = FALSE, lev, lev)

# RQ 1
colnam <- "equalized:mob"
descnam <- names(iv_normal)[1:3]
x <- cbind(iv_normal[, c(descnam, colnam)], iv_binomial[, colnam], iv_polr[, colnam],
  iv_weibull[, colnam], iv_cox[, colnam])
create_table(x, tabcolnams = c("Normal", "Binomial", "Multinomial", "Weibull", "Cox"),
  levA = NULL, levB = lev, hypothesisnam = "\\\\textbf{RQ 1}: Robinson vs. Naive")

# RQ 2
colnam <- "hybrid:mob"
x <- cbind(iv_normal[, c(descnam, colnam)], iv_binomial[, colnam], iv_polr[, colnam],
  iv_weibull[, colnam], iv_cox[, colnam])
create_table(x, tabcolnams = c("Normal", "Binomial", "Multinomial", "Weibull", "Cox"),
  levA = NULL, levB = lev, hypothesisnam = "\\\\textbf{RQ 2}: Robinson$\\_{\\\\hat{W}}$ vs. Naive")

# RQ 3
colnam <- "equalized:hybrid"
x <- cbind(iv_normal[, c(descnam, colnam)], iv_binomial[, colnam], iv_polr[, colnam],
  iv_weibull[, colnam], iv_cox[, colnam])
create_table(x, tabcolnams = c("Normal", "Binomial", "Multinomial", "Weibull", "Cox"),
  levA = NULL, lev, hypothesisnam = "\\\\textbf{RQ 3}: Robinson vs. Robinson$\\_{\\\\hat{W}}$")

# ---- Appendix A: Noncollapsibility ----
res <- resall[resall$setup %in% c(22:26),]
res <- res[res$problem == "binomial" | (res$problem == "weibull" & grepl("cox", res$algorithm)),]
HONESTY <- FALSE
source("helpers.R")

### Plots
plot_results_all(partB, cex = 1)

### Tables
lev <- c("pF_eta_x1_x2.mF_sin_x1_x5.tF_div_x1_x2" = "Setup A",
  "pF_eta_x1_x2.mF_sin_x1_x5.tF_div_x1_x2.X3" = "Setup A'",
  "0.5.mF_max_x1_x5.tF_log_x1_x2" = "Setup B",
  "pF_x2_x3.mF_log_x1_x3.1" = "Setup C",
  "pF_exp_x1_x2.mF_max2_x1_x5.tF_max_x1_x5" = "Setup D")

iv_binomial <- get_intervals(dataA = NULL, partB, model = "Binomial", gao = TRUE, lev, lev)
iv_cox <- get_intervals(dataA = NULL, partB, model = "Cox", gao = TRUE, lev, lev)
descnam <- names(iv_binomial)[1:3]

# RQ 4:
colnam <- "equalizedgao:equalized"
x <- cbind(iv_binomial[, c(descnam, colnam)], iv_cox[, colnam])
create_table(x, tabcolnams = c("Binomial", "Cox"), levA = NULL, levB = lev,
  hypothesisnam = "RQ 4: Gao vs. Robinson")

# RQ 5:
colnam <- "hybridgao:hybrid"
x <- cbind(iv_binomial[, c(descnam, colnam)], iv_cox[, colnam])
create_table(x, tabcolnams = c("Binomial", "Cox"), levA = NULL, levB = lev,
  hypothesisnam = "RQ 5: Gao$\\_\\{\\\\hat{W}}$ vs. Robinson$\\_{\\\\hat{W}}$")

# ---- Appendix B: Study based on Wager and Athey (2018) ----
res <- resall[resall$setup %in% c(1:16),]
res <- res[- grep("gao", res$algorithm), ]
HONESTY <- FALSE
source("helpers.R")

### Plots
plot_results_all(partA)
plot_results_all(partB, sc = scB)

### Tables
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

iv_normal <- get_intervals(partA, partB, model = "Normal", gao = FALSE, levA, levB)
iv_binomial <- get_intervals(partA, partB, model = "Binomial", gao = FALSE, levA, levB)
iv_polr <- get_intervals(partA, partB, model = "Multinomial", gao = FALSE, levA, levB)
iv_weibull <- get_intervals(partA, partB, model = "Weibull", gao = FALSE, levA, levB)
iv_cox <- get_intervals(partA, partB, model = "Cox", gao = FALSE, levA, levB)

### RQ 1
colnam <- "equalized:mob"
descnam <- names(iv_normal)[1:3]
x <- cbind(iv_normal[, c(descnam, colnam)], iv_binomial[, colnam], iv_polr[, colnam], iv_weibull[, colnam], iv_cox[, colnam])
create_table(x, tabcolnams = c("Normal", "Binomial", "Multinomial", "Weibull", "Cox"), levA, levB,
  hypothesisnam = "RQ 1: Robinson vs. Naive")

### RQ 2
colnam <- "hybrid:mob"
x <- cbind(iv_normal[, c(descnam, colnam)], iv_binomial[, colnam], iv_polr[, colnam],
  iv_weibull[, colnam], iv_cox[, colnam])
create_table(x, tabcolnams = c("Normal", "Binomial", "Multinomial", "Weibull", "Cox"),
  levA = NULL, levB = lev, hypothesisnam = "RQ 2: Robinson$\\_{\\\\hat{W}}$ vs. Naive")


### RQ 3
colnam <- "equalized:hybrid"
x <- cbind(iv_normal[, c(descnam, colnam)], iv_binomial[, colnam], iv_polr[, colnam],
  iv_weibull[, colnam], iv_cox[, colnam])
create_table(x, tabcolnams = c("Normal", "Binomial", "Multinomial", "Weibull", "Cox"),
  levA, levB, hypothesisnam = "RQ 3: Robinson vs. Robinson$\\_{\\\\hat{W}}$")

