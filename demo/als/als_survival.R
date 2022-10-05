####### ALS data example - survival time ######

source("setup.R")

#---- Data preparation ----
# Data must be downloaded from https://nctu.partners.org/ProACT
# Please run the corresponding demo from the R package TH.data (`demo("PROACT", package = "TH.data")`).
# Store the data in a subfolder called /data.
load("data/RALSfinal.rda")

del <- grepl(".halfYearAfter", names(data))
ALSsurvdata <- data[ , -which(del)]

### delete all rows where there is no survival info
delete <- sapply(ALSsurvdata$survival.time, is.na)
table(delete)
ALSsurvdata <- ALSsurvdata[ -which(delete) , ]


### delete columns with more than 50% NAs
### except scores, t.onsettrt, Riluzole
keepvarnames <- c("survival.time", "cens", "t.onsettrt", "Riluzole")
keepvars <- grepl(paste0(keepvarnames, collapse = "|"), names(ALSsurvdata))


pNA <- 0.5 * nrow(ALSsurvdata)
ALSsurvdata <- ALSsurvdata[ , keepvars | (colSums(is.na(ALSsurvdata)) < pNA)]
names(ALSsurvdata)

keepvars1 <- grepl(paste0(keepvarnames, collapse = "|"), names(ALSsurvdata))
dim(na.omit(ALSsurvdata[ , !keepvars1]))

### delete all of non-useful variables
delete <- grepl("Delta|delta|SubjectID|Unit|Onset|treatment.group", names(ALSsurvdata))
ALSsurvdata <- ALSsurvdata[ , -which(delete)]

ALSsurvdata <- ALSsurvdata[complete.cases(ALSsurvdata[, c("survival.time", "cens", "Riluzole")]),]
ALSsurvdata$survival.time[ALSsurvdata$survival.time == 0] <- 0.1

### dummy encode Riluzole treatment variable
ALSsurvdata$Riluzole <- c(0, 1)[as.numeric(ALSsurvdata$Riluzole)]

# ### Rename variables for plotting
names(ALSsurvdata)[names(ALSsurvdata) == "t.onsettrt"] <- "time_onset_treatment"
Z <- names(ALSsurvdata)[!(names(ALSsurvdata) %in% c("cens", "survival.time", "Riluzole"))]
names(ALSsurvdata)[names(ALSsurvdata) %in% Z] <- tolower(names(ALSsurvdata)[names(ALSsurvdata) %in% Z])
names(ALSsurvdata) <- gsub("fam.hist.", "family_history_", names(ALSsurvdata))

### Censor observations at t = 750
id.750 <- ALSsurvdata$survival.time >= 750
sum(id.750) # 5 obs
ALSsurvdata[id.750, "cens"] <- 0
ALSsurvdata[id.750, "survival.time"] <- 750

rm(data)
dim(ALSsurvdata)
save(ALSsurvdata, file = "data/ALSsurvdata.rda")

#---- Fit Naive and Robinson forest ----

### Fit Naive forest
# Get base model
bm_alssurv <- coxph(Surv(survival.time, cens) ~ Riluzole, data = ALSsurvdata)
summary(bm_alssurv)
logLik(bm_alssurv) # nocens: -6834.855 #cens: -6830.705

# Forest setup
xnams <- names(ALSsurvdata)
xnams <- xnams[!xnams %in% all.vars(bm_alssurv$terms)]
zsurv <- as.formula(paste("~", paste(xnams, collapse =  " + ")))
# set mtry based on number of splitting variables
mtry <- min(ceiling(sqrt(length(xnams)) + 20), length(xnams))

# Build forest
set.seed(2007)
frst_alssurv <- fit_forest(bm_alssurv, data = ALSsurvdata, zformula = zsurv)

### Fit Robinson forest
# Get propensity scores
X <- dummy_cols(ALSsurvdata[, xnams], remove_first_dummy = TRUE, remove_selected_columns = TRUE)

set.seed(2007)
cf <- causal_forest(X = X, Y = rep(0, nrow(ALSsurvdata)), W = ALSsurvdata$Riluzole,
  min.node.size = min_size_group, sample.fraction = prt$fraction,
  mtry = mtry, # ci.group.size = 1,
  num.trees = NumTrees, honesty = FALSE)
What <- cf$W.hat

# Get natural parameter estimates
ALSsurvdata1 <- ALSsurvdata[ALSsurvdata$Riluzole == 1,]
ALSsurvdata1$Riluzole <- NULL
ALSsurvdata0 <- ALSsurvdata[ALSsurvdata$Riluzole == 0,]
ALSsurvdata0$Riluzole <- NULL

set.seed(2007)
boost.Y0 <- gbm(Surv(survival.time, cens) ~ ., data = ALSsurvdata0, distribution = "coxph")
boost.Y1 <- gbm(Surv(survival.time, cens) ~ ., data = ALSsurvdata1, distribution = "coxph")
lp0 <- predict(boost.Y0, newdata = ALSsurvdata, type = "link")
lp1 <- predict(boost.Y1, newdata = ALSsurvdata, type = "link")
nu <- What *lp1 + (1-What) * lp0

# Update base model
origRiluzole <- ALSsurvdata$Riluzole
ALSsurvdata$trt <- origRiluzole - What
ALSsurvdata$Riluzole <- NULL
ALSsurvdata$nu <- nu
bm_alssurv_robinson <- coxph(Surv(survival.time, cens) ~ trt + offset(nu), data = ALSsurvdata)
summary(bm_alssurv_robinson)
logLik(bm_alssurv_robinson) # without: -6632.395 # with: -6624.68

# Build forest
set.seed(2007)
frst_alssurv_robinson <- fit_forest(bm_alssurv_robinson, data = ALSsurvdata, zformula = zsurv)

#--- Visualizations ----

### Propensity scores
ALSsurvdata$Riluzole <- origRiluzole

pwhat <- ggplot(ALSsurvdata, aes(x = What, fill = as.factor(Riluzole),
  color = as.factor(Riluzole))) +
  stat_density(alpha = 0.2, position = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(c(0, 9)) +
  xlim(c(-0.001, 1.001)) +
  xlab(expression(hat(pi)(x))) +
  scale_color_manual(values = cols, name = "Riluzole", labels = c("No", "Yes")) +
  scale_fill_manual(values = cols, name = "Riluzole", labels = c("No", "Yes")) +
  labs(fill = "Riluzole")
pwhat

trtrobinson <- origRiluzole - What
set.seed(2007)
update.What <- causal_forest(X = as.matrix(X),
  Y = rep(0, nrow(X)), W = trtrobinson,
  min.node.size = min_size_group,
  mtry = mtry,
  num.trees = NumTrees, honesty = FALSE)$W.hat

pwhatup <- ggplot(ALSsurvdata, aes(x = update.What, fill = as.factor(Riluzole),
  color = as.factor(Riluzole))) +
  stat_density(alpha = 0.2, position = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(c(0, 9)) +
  xlab(expression(paste("estimated ", w-hat(pi)(x)))) +
  scale_color_manual(values = cols, name = "Riluzole", labels = c("No", "Yes")) +
  scale_fill_manual(values = cols, name = "Riluzole", labels = c("No", "Yes")) +
  labs(fill = "Riluzole")
pwhatup

ppi <- ggarrange(pwhat, pwhatup, common.legend = TRUE, ncol = 2)
ppi


### Kaplan-Meier curves
m <- survfit(Surv(survival.time, cens) ~ Riluzole, data = ALSsurvdata)
gsur <- ggsurvplot(m, data = ALSsurvdata, size = 1,                 # change line size
  palette =  cols,# custom color palettes
  # conf.int = TRUE,          # Add confidence interval
  ggtheme = mytheme,      # Change ggplot2 theme
  legend.title = "Riluzole",
  legend.labs = c("No", "Yes"),
  ylab = "survival probability",
  xlim = c(0, 750)
)
gsur

### Distribution/density plot of treatment effects
dp_alssurv <- cbind(data.frame(tau = c(frst_alssurv$pm)), ALSsurvdata)
dp_alssurv_robinson <- cbind(data.frame(tau = c(frst_alssurv_robinson$pm)), ALSsurvdata)
taus_surv <- data.frame(Naive = dp_alssurv$tau, Robinson = dp_alssurv_robinson$tau)
naive_mean <- mean(taus_surv$Naive)
robinson_mean <- mean(taus_surv$Robinson)
p_alssurv <- ggplot(taus_surv) +
  geom_density(aes(x = Robinson, fill = "Robinson") , alpha = 0.4) +
  geom_density(aes(x = Naive, fill = "Naive"), alpha = 0.4) +
  mytheme +
  ylab("density") +
  xlim(-0.8, 0.8) +
  xlab(expression(hat(tau)(x))) +
  scale_fill_manual(values = c(Naive = "darkorange", Robinson = "darkblue")) +
  theme(legend.title = element_blank()) +
  geom_vline(mapping = aes(xintercept = naive_mean), linetype = "dashed", col = "darkorange") +
  geom_vline(mapping = aes(xintercept = robinson_mean), linetype = "dashed", col = "darkblue") +
  geom_text(aes(x = naive_mean, label = "\nNaive", y = 0.5), colour = "black", angle = 90) +
  geom_text(aes(x = robinson_mean, label = "\nRobinson", y = 0.9), colour = "black", angle = 90)

p_alssurv

### Dependence plots
xnams <- names(ALSsurvdata)[c(2:19)]
all_dp <- c(dp_alssurv_robinson$tau, dp_alssurv$tau)
taumax <- max(all_dp)
taumin <- min(all_dp)
for (nam in xnams) {
  plot_side(nam, dp1 = dp_alssurv, dp2 = dp_alssurv_robinson, ymin = taumin, ymax = taumax)
}

