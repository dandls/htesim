####### ALS data example - survival time ######

source("setup.R")

#---- Data preparation ----
# Data must be downloaded from https://nctu.partners.org/ProACT
# Please run the corresponding demo from the R package TH.data (`demo("PROACT", package = "TH.data")`).
# Store the data in a subfolder called /data.
load("data/RALSfinal.rda")
ALSFRSdata <- data

# extract target = handwriting ALSFR score = V_1217
hya <- grep(".halfYearAfter", names(ALSFRSdata), value = TRUE)
del <- hya[-grep("1217", hya)]
st <- grep(".Start", names(ALSFRSdata), value = TRUE)
del <- c(del, st[-grep("1217", hya)])
del <- c(del,"survival.time", "cens")
ALSFRSdata <- ALSFRSdata[ , !(names(ALSFRSdata) %in% del)]
ALSFRSdata <- ALSFRSdata[complete.cases(ALSFRSdata[, c("V_1217.halfYearAfter",
  "V_1217.Start", "Riluzole")]),]


### delete all of non-useful variables
delete <- grepl("Delta|delta|SubjectID|Unit|Onset|treatment.group", names(ALSFRSdata))
ALSFRSdata <- ALSFRSdata[ , -which(delete)]

### delete Basophil count and Hematocrit (non-explicable values)
ALSFRSdata$Value_Absolute_Basophil_Count <- NULL
sm1 <- ALSFRSdata$Value_Hematocrit < 1 & !is.na(ALSFRSdata$Value_Hematocrit)
ALSFRSdata$Value_Hematocrit[sm1] <- ALSFRSdata$Value_Hematocrit[sm1] * 100

### delete extremely high Phosphorus value
ALSFRSdata$Value_Phosphorus[!is.na(ALSFRSdata$Value_Phosphorus) &
    ALSFRSdata$Value_Phosphorus > 5] <- NA

### delete columns with more than 50% NAs
### except scores, t.onsettrt, Riluzole
keepvarnames <- c("V_1217.halfYearAfter", "V_1217.Start", "t.onsettrt", "Riluzole")
keepvars <- grepl(paste0(keepvarnames, collapse = "|"), names(ALSFRSdata))

pNA <- 0.5 * nrow(ALSFRSdata)
ALSFRSdata <- ALSFRSdata[ , keepvars | (colSums(is.na(ALSFRSdata)) < pNA)]

### dummy encode Riluzole treatment variable
ALSFRSdata$Riluzole <- c(0, 1)[as.numeric(ALSFRSdata$Riluzole)]

### Rename variables for plotting
names(ALSFRSdata)[names(ALSFRSdata) == "t.onsettrt"] <- "time_onset_treatment"

### Rename target variable to be more descriptive
names(ALSFRSdata)[names(ALSFRSdata) == "V_1217.Start"] <- "writing.Start"
names(ALSFRSdata)[names(ALSFRSdata) == "V_1217.halfYearAfter"] <- "writing.halfYearAfter"
names(ALSFRSdata)[names(ALSFRSdata) == "SubjectLiters_fvc"] <- "FVC"

Z <- names(ALSFRSdata)[!(names(ALSFRSdata) %in% c("writing.halfYearAfter",
  "Riluzole",
  "writing.Start"))]
names(ALSFRSdata)[names(ALSFRSdata) %in% Z] <-
  tolower(names(ALSFRSdata)[names(ALSFRSdata) %in% Z])
names(ALSFRSdata) <- gsub("fam.hist.", "family_history_", names(ALSFRSdata))

dim(ALSFRSdata)
save(ALSFRSdata, file = "data/ALSFRSdata.rda")


#---- Fit Naive and Robinson forest ----
### Fit Naive forest
# Get base model
bm_alsfrs <- polr(writing.halfYearAfter ~ Riluzole,
  data = ALSFRSdata)
summary(bm_alsfrs)
logLik(bm_alsfrs)

# Forest setup
xnams <- names(ALSFRSdata)
xnams <- xnams[!xnams %in% c(all.vars(bm_alsfrs$terms))]
zalsfrs <- as.formula(paste("~", paste(xnams, collapse =  " + ")))
# set mtry based on number of splitting variables
mtry <- min(floor(sqrt(length(xnams)) + 20), length(xnams))

# Build forest
set.seed(2007)
frst_alsfrs <- fit_forest(bm_alsfrs, data = ALSFRSdata, zformula = zalsfrs)

### Fit Robinson forest
# Get propensity scores
X <- ALSFRSdata[, xnams]
X$writing.Start <- factor(X$writing.Start, ordered = FALSE)
X <- dummy_cols(X, remove_first_dummy = TRUE, remove_selected_columns = TRUE)

set.seed(2007)
cf <- causal_forest(X = X, Y = rep(0, nrow(ALSFRSdata)), W = ALSFRSdata$Riluzole,
  min.node.size = min_size_group, sample.fraction = prt$fraction,
  mtry = mtry, # ci.group.size = 1,
  num.trees = NumTrees, honesty = FALSE)
What <- cf$W.hat

# Get natural parameter estimates
ALSFRSdata1 <- ALSFRSdata[ALSFRSdata$Riluzole == 1,]
ALSFRSdata1$Riluzole <- NULL
ALSFRSdata0 <- ALSFRSdata[ALSFRSdata$Riluzole == 0,]
ALSFRSdata0$Riluzole <- NULL

set.seed(2007)
boost.Y0 <- blackboost(writing.halfYearAfter ~ ., data = ALSFRSdata0, family = PropOdds())
boost.Y1 <- blackboost(writing.halfYearAfter ~ ., data = ALSFRSdata1, family = PropOdds())
lp0 <- as.numeric(predict(boost.Y0, newdata = ALSFRSdata, type = "link"))
lp1 <- as.numeric(predict(boost.Y1, newdata = ALSFRSdata, type = "link"))
nu <- What * lp1 + (1-What) * lp0

# Update base model
origRiluzole <- ALSFRSdata$Riluzole
ALSFRSdata$trt <- ALSFRSdata$Riluzole - What
ALSFRSdata$Riluzole <- NULL
ALSFRSdata$nu <- nu
bm_alsfrs_robinson <- polr(writing.halfYearAfter ~ trt + offset(nu),
  data = ALSFRSdata)
summary(bm_alsfrs_robinson)
logLik(bm_alsfrs_robinson) # -2858.673 (df = 5)

# Build forest
set.seed(2007)
frst_alsfrs_robinson <- fit_forest(bm_alsfrs_robinson, data = ALSFRSdata, zformula = zalsfrs)

#--- Visualizations ----
library("reshape")
library("plyr")
library("ggalluvial")

### Propensity scores
ylim <- c(0, 12.01)
ALSFRSdata$Riluzole = origRiluzole

pwhat <- ggplot(ALSFRSdata, aes(x = What, fill = as.factor(Riluzole), color = as.factor(Riluzole))) +
  stat_density(alpha = 0.2, position = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(ylim) +
  xlim(c(-0.001, 1.001)) +
  xlab(expression(hat(pi)(x, y[0]))) +
  scale_color_manual(values = cols, name = "Riluzole", labels = c("No", "Yes")) +
  scale_fill_manual(values = cols, name = "Riluzole", labels = c("No", "Yes")) +
  labs(fill = "Riluzole")
pwhat

set.seed(2007)
update.What <- causal_forest(X = as.matrix(X),
  Y = rep(0, nrow(X)), W = ALSFRSdata$trt,
  min.node.size = min_size_group,
  mtry = mtry,
  num.trees = NumTrees, honesty = )$W.hat

pwhatup <- ggplot(ALSFRSdata, aes(x = update.What, fill = as.factor(Riluzole),
  color = as.factor(Riluzole))) +
  stat_density(alpha = 0.2, position = "identity") +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(ylim) +
  xlab(expression(paste("estimated ", w-hat(pi)(x, y[0])))) +
  scale_color_manual(values = c("0" = cols[1], "1" = cols[2]), name = "Riluzole", labels = c("No", "Yes")) +
  scale_fill_manual(values = c("0" = cols[1], "1" = cols[2]), name = "Riluzole", labels = c("No", "Yes")) +
  labs(fill = "Riluzole")
pwhatup

ppi <- ggarrange(pwhat, pwhatup, common.legend = TRUE)
ppi

### Distribution of handwriting ability scores
df.new <- ddply(ALSFRSdata, .(Riluzole), summarise,
  prop = prop.table(table(writing.halfYearAfter)),
  writing.score = names(table(writing.halfYearAfter)))
df.new$Riluzole <- as.factor(df.new$Riluzole)

galsfrs1 <- ggplot(df.new,aes(writing.score, prop, fill = Riluzole))+
  geom_bar(stat="identity",position='dodge') +
  scale_fill_manual(values = cols, name = "Riluzole", labels = c("No", "Yes")) +
  ylab("relative frequency") +
  xlab(expression(paste("Handwriting ability score ",(Y[6])))) +
  ylim(c(0, 0.6)) +
  mytheme
galsfrs1

# Alluvial plots
avdat <- as.data.frame(ftable(xtabs(~ Riluzole + writing.Start + writing.halfYearAfter,
  data=ALSFRSdata)))
ggplot(avdat,
  aes(y = Freq, axis1 = writing.halfYearAfter, axis2 = writing.Start)) +
  geom_alluvium(aes(fill = Riluzole), width = 1/12, alpha = 0.75) +
  geom_stratum(width = 1/12) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum), label.size = 0.1)) +
  scale_x_discrete(limits = c("6 months", "Start"), expand = c(.05, .05)) +
  scale_fill_manual(values = cols, name = "Riluzole", labels = c("No", "Yes")) +
  coord_flip() +
  mytheme


### Distribution/density plot of treatment effects
dp_alsfrs <- cbind(data.frame(tau = frst_alsfrs$pm), ALSFRSdata)
dp_alsfrs_robinson <- cbind(data.frame(tau = frst_alsfrs_robinson$pm), ALSFRSdata)

taus <- data.frame(Naive = dp_alsfrs$tau, Robinson = dp_alsfrs_robinson$tau)
mean_alsfrs <- mean(taus$Naive)
mean_alsfrs_robinson <- mean(taus$Robinson)
p_alsfrs <- ggplot(taus) +
  geom_density(aes(x = Robinson, fill = "Robinson"), alpha = 0.4) +
  geom_density(aes(x = Naive, fill = "Naive"), alpha = 0.4) +
  mytheme +
  ylab("density") +
  xlab(expression(hat(tau)(x, y[0]))) +
  scale_fill_manual(values = c("darkorange", "darkblue")) +
  theme(legend.title = element_blank()) +
  geom_vline(aes(xintercept = mean_alsfrs), linetype = "dashed", col = "darkorange") +
  geom_vline(aes(xintercept = mean_alsfrs_robinson), linetype = "dashed", col = "darkblue") +
  geom_text(aes(x = mean_alsfrs, label = "\nNaive", y = 1), colour = "black", angle = 90) +
  geom_text(aes(x = mean_alsfrs_robinson - 0.06, label = "\nRobinson", y = 1), colour = "black", angle = 90)

p_alsfrs

### Dependence plot
taumax <- max(c(dp_alsfrs_robinson$tau, dp_alsfrs$tau))
taumin <- min(c(dp_alsfrs_robinson$tau, dp_alsfrs$tau))

for (nam in xnams) {
  plot_side(nam, dp_alsfrs, dp_alsfrs_robinson, ymin = taumin, ymax = taumax, y0 = TRUE)
}

