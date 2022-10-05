#--- Libraries ----
library("model4you") # to fit forests
library("ggplot2") # for plotting
library("grf") # to get propensity scores
library("fastDummies") # create dummy encoded features for grf
library("gbm") # to estimate offset for coxph
library("MASS") # for polr model
library("survival") # for coxph model
library("mboost") # to estimate offset for polr
library("colorspace") # colors for plotting
library("stringr") # for label generation on plot
library("ggpubr") # arrange ggplots in columns
library("gridExtra")
library("survminer")

set.seed(2007)
dir.create("results", showWarnings = FALSE)

#--- Setups for trees/forests ----
NumTrees <- 500L
min_size_group <- 7L
min_node_size <- min_size_group*2L
ctrl <- ctree_control(testtype = "Univ", minsplit = 2,
  minbucket = min_node_size,
  mincriterion = 0, saveinfo = FALSE)
# stopping criterion: min number of observations for one treatment group
converged_crit <- function(data, weights, control) {
  function(subset, weights)  {
    trt <- data$data[, data$variables$x]
    if (length(unique(trt)) == 2 | is.factor(trt)) {
      return(all(table(trt[subset]) > min_node_size))
    } else {
      # hybrid: trt not 0/1 but between -1 and 1
      return(all(table(trt[subset] > 0) > min_size_group))
    }
  }
}

min_update <- 20L
prt <- list(replace = FALSE, fraction = .5)

fit_forest <- function(bm, data, zformula) {
  frst <- pmforest(bm, data = data, zformula = zformula, ntree = NumTrees, perturb = prt, mtry = mtry,
    control = ctrl, trace = TRUE)
  pm <- pmodel(frst)
  return(list(frst = frst, pm = pm))
}


#---- Setups for plots ----
mytheme <- theme_bw()

cols <- c("#DBB267", "#57967A")

dp_plot <- function(var, data, ymin, ymax, y0) {
  p <- ggplot(data, aes_string(y = "tau", x = var))

  if(is.factor(data[, var])) {
    p <- p + geom_boxplot(varwidth = TRUE, fill = "lightgrey") +
      stat_summary(fun = mean, aes(colour = "mean"), geom = "point",
        shape = 18, size = 3) +
      scale_color_manual(values = c("mean" = "#619CFF")) +
      theme(legend.position = "none") +
      scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "),
        width = 20))
  } else {
    p <- p + geom_point(alpha = 0.2) +
      geom_smooth(se = FALSE, aes(color = "smooth curve")) +
      scale_color_manual(values = c("smooth curve" = "#619CFF")) +
      theme(legend.position = "none")
  }
    p <- p + ylim(c(ymin, ymax)) +
      mytheme +
      theme(legend.position = "none")
    if (y0) {
      p <- p +  ylab(expression(hat(tau)(x, y[0])))
    } else {
      p <- p +  ylab(expression(hat(tau)(x)))
    }
   return(p)
}

plot_side <- function(xvar, dp1, dp2, dp3 = NULL, ymin, ymax, y0 = FALSE) {
  p1 <- dp_plot(xvar, dp1, ymin = ymin, ymax = ymax, y0 = y0) +
    xlab("Naive")
  p2 <- dp_plot(xvar, dp2, ymin = ymin, ymax = ymax, y0 = y0) +
    xlab("Robinson")
  if (!is.null(dp3)) {
    p3 <- dp_plot(xvar, dp3, ymin = ymin, ymax = ymax, y0 = y0) +
      xlab("Gao")
    return(grid.arrange(p1, p2, p3, ncol = 3))
  }
  return(grid.arrange(p1, p2, ncol = 2))
}
