### Postpartum blood loss data example

#---- Setup ----
library("model4you")
library("survival")
library("ggplot2")
library("trtf")
library("tram")
library("multcomp")
library("lattice")
library("tidyr")
library("stringr")
library("grf")
library("fastDummies")
library("colorspace")
library("htesim")
library("ggpubr")

# Load data
data(blood)

# Helper functions
# Draw dependency plot
dp_plot <- function(var, data, i, xlim = NULL, beta, ytxt = "odds ratio") {
  p <- ggplot(data, aes_string(y = beta, x = var)) +
    ylab(ytxt) +
    theme
  #
  if(is.factor(data[, var])) {
    p <- p + geom_boxplot(varwidth = TRUE)
  } else {
    p <- p + geom_point(alpha = 0.2) + geom_smooth(se = FALSE)
  }
  if (!is.null(xlim)) {
    p <- p + xlim(xlim)
  }
  return(p)
}

#--- Data Processing ----

# set NA in Dauer to 0
blood$DAUER.ap[is.na(blood$DAUER.ap)] <- 0

# Define binary treatment variable VCmode (vaginal delivery vs. cesarean section)
blood$mode <- with(blood, SECTIO.prim == "yes" | SECTIO.sek == "yes" |
    SECTIO.not == "yes") + 1
blood$mode[blood$SECTIO.sek == "yes" | blood$SECTIO.not == "yes"] <- 3
blood$mode <- factor(blood$mode, levels = 1:3,
  labels = c("Vaginal delivery", "Planned Cesarean", "Unplanned Cesarean"))
blood$VCmode <- blood$mode
levels(blood$VCmode) <- c("vaginal", "cesarean", "cesarean")

# VCmode as dummy variable (necessary for causal forests)
blood$VCmodedummy <- c(0, 1)[blood$VCmode]

# create intervals of blood loss
### interval length: 50 for MBL < 1000; 100 for MBL > 1000
off <- 25
tm1 <- with(blood, ifelse(MBL < 1000, MBL - off, MBL - 2 * off)) # SD: lower interval limit
tm2 <- with(blood, ifelse(MBL >= 1000, MBL + 2 * off, MBL + off)) # SD: upper interval  limit
blood$MBLsurv <- Surv(time = tm1, time2 = tm2, type = "interval2")

# define which patient characteristics used as splitting variables
x <- c("GA", "AGE", "MULTIPAR", "BMI", "TWIN", "NW", "IOL", "AIS")
xfm <- paste(x, collapse = "+") # prepartum variables

# define formula
xfm_MBL <- as.formula(paste("MBLsurv ~ VCmodecenter |", xfm))

# Remove cases with missing values in considered variables
xMBL_cc <- complete.cases(blood[, c(x, "MBLsurv", "VCmode")])
blood <- blood[xMBL_cc,]

# Remove case with extreme MBL value
mid <- which.max(blood$MBL)
blood <- blood[-mid, ]


#--- Setups for plots ---
# limits of MBL
qy  <- 0:max(blood$MBL)
MBLlim <- c(0, 2700)

# ggplot theme
theme <- theme_classic()

# colors
tcols <- diverge_hcl(50, h = c(246, 40), c = 96, l = c(65, 90), alpha = .5)
cols <- qualitative_hcl(3, palette = "Harmonic")

#--- Setups for trees/forests ----
mtry <- min(floor(sqrt(length(x)) + 20), length(x))
NumTrees <- 500
min_size_group <- 7L
min_node_size <- min_size_group*2L

# Additional setups for model-based forest to be equal to causal forest
min_update <- min_node_size
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
prt <- list(replace = FALSE, fraction = .5)

#--- Causal forest ----
# Dummy encode splitting variables
x_dummy <- dummy_cols(blood[, x], remove_first_dummy = TRUE, remove_selected_columns = TRUE)

#--- Personalized Model ---
set.seed(1234L)
cf <- causal_forest(X = as.matrix(x_dummy),
  Y = blood$MBL, W = blood$VCmodedummy,
  min.node.size = min_size_group,
  mtry = mtry,
  num.trees = NumTrees, honesty = FALSE)
W.hat <- cf$W.hat

# Center treatment indicator
blood$VCmodecenter <- blood$VCmodedummy - W.hat

#--- Base Model ---
lin_MBL <- lm(MBL ~ VCmodecenter, data = blood)
summary(lin_MBL)
tau_lin <- coef(lin_MBL)["VCmodecenter"]
alpha_lin <- coef(lin_MBL)["(Intercept)"]

# save(W.hat, Y.hat, file = "What_Yhat.rda")
pred <- predict(cf, estimate.variance = TRUE)  # OUT-OF-BAG
pred_df <- data.frame(pred)

p_cf <- ggplot(pred_df, aes(predictions)) +
  geom_density() +
  theme +
  xlab(expression(hat(tau)(x[i]))) +
  geom_vline(aes(xintercept=tau_lin, linetype = "base model"), colour = "red") +
  scale_linetype_manual(name = "", values = 2, guide = guide_legend(override.aes = list(color = c("red")))) +
  theme(legend.position = c(0.25, 0.8), legend.title = element_blank())
p_cf

#--- Dependence plots ----
# Replace . by _ in splitting variable names
var_shown <- x
propnams <- str_replace(var_shown, pattern = "\\.", replacement = "_")
nvar <- length(var_shown)

dpcf <- cbind(pred, blood)

ps <- lapply(seq_len(nvar), function(i) dp_plot(var_shown[i], dpcf, i = "", beta = "predictions", ytxt = expression(hat(tau))))
propnams <- str_replace(var_shown, pattern = "\\.", replacement = "_")
ps

#---- Plot estimated propensity scores and estimated centered treatments-----
#
# Estimate centered treatments
update.What <- causal_forest(X = as.matrix(x_dummy),
  Y = blood$MBL, W = blood$VCmodecenter,
  min.node.size = min_size_group,
  mtry = mtry,
  num.trees = NumTrees, honesty = FALSE)$W.hat

# plot density of propensity scores
pwhat <- ggplot(blood, aes(x = W.hat, color = VCmode)) +
  stat_density(geom = "line", position = "identity") +
  theme +
  theme(legend.position = "none", legend.title = element_blank()) +
  ylim(c(0, 9)) + xlab(expression(hat(pi)(x[i]))) +
  scale_color_manual(values = cols)

# plot density of estimated centered treatments
pwhatup <- ggplot(blood, aes(x = update.What, color = VCmode)) +
  # geom_density(show_guide = FALSE) +
  stat_density(geom = "line", position = "identity") +
  theme + theme(legend.title = element_blank()) +
  ylim(c(0, 9)) + xlab(expression(paste("estimated ", w[i]-hat(pi)(x[i])))) +
  scale_color_manual(values = cols)

ppi <- ggarrange(pwhat, pwhatup, common.legend = TRUE)
ppi

#--- Model-based forests ----

### Fit base model & analyse treatment effect estimate
m_MBL <- BoxCox(MBLsurv ~ VCmodecenter, data = blood,
  bounds = c(0, Inf), support = c(250, 2000))
logLik(m_MBL)
coef(m_MBL)
confint(m_MBL)

# Distribution of measured blood loss (mL) stratified by
# mode of delivery. Rugs indicate measured blood loss observations,
# stratified by mode of delivery.
nd <- data.frame(VCmodecenter = unique(blood$VCmodedummy))
par(las = 1)
plot(as.mlt(m_MBL), newdata = nd,
  q = qy, type = "distribution", col = cols, lwd = 3, xlim = MBLlim,
  xlab = "MBL", ylab = "Probability", ylim = c(-.05, 1.05))
rug(blood$MBL[blood$VCmodedummy == 0], lwd = 2, col = cols[1])
rug(blood$MBL[blood$VCmodedummy == 1], side = 3, lwd = 2, col = cols[2])
legend("bottomright", lwd = 2, col = cols, legend = levels(blood$VCmode), bty = "n")


# Fit personalized model
set.seed(290875)
rf <- traforest(as.mlt(m_MBL), formula = xfm_MBL, data = blood, ntree = NumTrees,
  perturb = prt, trace = TRUE, converged = converged_crit,
  mtry = mtry, control = ctrl, min_update = min_update
)

# Predict on training data
pml <- predict(rf, type = "coef", OOB = TRUE)
pm <- do.call("rbind", pml)
weights <- predict(rf, type = "weights", OOB = TRUE)

### Distribution of coefficients:
pmd <- as.data.frame(pm)
names(pmd)[length(names(pmd))] <- "tau"

ptau <- ggplot(pmd, aes(tau)) +
  geom_density() +
  theme +
  xlab(expression(hat(tau)(x[i]))) +
  geom_vline(aes(xintercept=coef(m_MBL), linetype = "base model"), colour = "red") +
  scale_linetype_manual(name = "", values = 2, guide = guide_legend(override.aes = list(color = c("red")))) +
  theme(legend.position = c(0.2, 0.8), legend.title = element_blank())
ptau

# Get likelihood
lk <- logLik(rf, OOB = FALSE) # -3539.598

#---- Dependence Plots of treatment effect tau----
dp <- cbind(pmd, blood)

print(var_shown)
ps <- lapply(seq_len(nvar), function(i) dp_plot(var_shown[i], dp, i = "", beta = "tau",
  ytxt = expression(hat(tau))))
ps


#--- Dependency plots of prognostic effect/intercept alpha ----
mnd <- data.frame(trt = 0)
names(mnd) <- "VCmodecenter"

# get median MBL for each observation
med <- predict(rf, newdata = blood, OOB = TRUE,
  mnewdata = mnd, type = "quantile", prob = 0.5)
dp$median <- c(do.call(rbind, med))

ps_alpha <- lapply(seq_len(nvar), function(i) dp_plot(var_shown[i], data = dp, i = "",
  beta = "median", ytxt = "Median(MBL|w = vaginal)"))
propnams <- str_replace(var_shown, pattern = "\\.", replacement = "_")
ps_alpha