# Translate method names to formula (expression) if necessary
lookup <- c(
  "mob" = "mob",
  "hybrid" = expression("mob"(widehat(W))),
  "equalized" = expression("mob"(widehat(W), widehat(Y))),
  "mobcf" = "mobcf",
  "cf" = "cf"
  )

if (HONESTY) {
  lookup2 <- c(
    "mob-honest" = "h-mob",
    "hybrid-honest" = expression("h-mob"(widehat(W))),
    "equalized-honest" = expression("h-mob"(widehat(W), widehat(Y))),
    "mobcf-honest" = "h-mobcf",
    "cf-honest" = "h-cf"
  )
  lookup <- c(lookup, lookup2)
}

# Color scheme
c(
  "mob" = "snow", "mob-honest" = "#888888", # white
  "hybrid" = "#E69F00",  "hybrid-honest" = "#D55E00", #red
  "equalized" = "#CC79A7", "equalized-honest" = "#882255", #rosa
  "mobcf" = "#999933", "mobcf-honest" = "#117733", # green
  "cf" = "#56B4E9", "cf-honest" = "#0072B2" # blue
)

# Naming convention
exA <- c("pF_x3.mF_x3.0" = expression(mu(x[3]) + 0 %.% W(x[3])),
  "0.5.0.tF_exp_x1_x2" = expression(tau(x[1], x[2]) * W),
  "pF_x1.mF_x1.tF_exp_x1_x2" = expression(mu(x[1]) + tau(x[1], x[2]) * W(x[1])),
  "0.5.mF_x1.tF_exp_x1_x2" = expression(mu(x[1]) + tau(x[1], x[2]) * W),
  "0.5.mF_x3.tF_exp_x1_x2" = expression(mu(x[3]) + tau(x[1], x[2]) * W),
  "pF_x3.mF_x3.tF_exp_x1_x2" = expression(mu(x[3]) + tau(x[1], x[2]) * W(x[3])),
  "pF_x3.0.tF_exp_x1_x2" = expression(tau(x[1], x[2])* W(x[3])),
  "pF_x4.mF_x3.tF_exp_x1_x2" = expression(mu(x[3]) + tau(x[1], x[2]) * W(x[4])),
  "pF_eta_x1_x2.mF_sin_x1_x5.tF_div_x1_x2" = "Setup A",
  "0.5.mF_max_x1_x5.tF_log_x1_x2" = "Setup B",
  "pF_x2_x3.mF_log_x1_x3.1" = "Setup C",
  "pF_exp_x1_x2.mF_max2_x1_x5.tF_max_x1_x5" = "Setup D")

exB <- c("pF_x3.mF_x3.0" = expression(mu(x[3]) + 0 %.% (W(x[3]) - .5)),
  "0.5.0.tF_exp_x1_x2" = expression(tau(x[1], x[2]) * (W - .5)),
  "pF_x1.mF_x1.tF_exp_x1_x2" = expression(mu(x[1]) + tau(x[1], x[2]) * (W(x[1]) - .5)),
  "0.5.mF_x1.tF_exp_x1_x2" = expression(mu(x[1]) + tau(x[1], x[2]) * (W - .5)),
  "0.5.mF_x3.tF_exp_x1_x2" = expression(mu(x[3]) + tau(x[1], x[2]) * (W - .5)),
  "pF_x3.mF_x3.tF_exp_x1_x2" = expression(mu(x[3]) + tau(x[1], x[2]) * (W(x[3]) - .5)),
  "pF_x3.0.tF_exp_x1_x2" = expression(tau(x[1], x[2])* (W(x[3]) - .5)),
  "pF_x4.mF_x3.tF_exp_x1_x2" = expression(mu(x[3]) + tau(x[1], x[2]) * (W(x[4]) - .5)),
  "pF_eta_x1_x2.mF_sin_x1_x5.tF_div_x1_x2" = "Setup A",
  "0.5.mF_max_x1_x5.tF_log_x1_x2" = "Setup B",
  "pF_x2_x3.mF_log_x1_x3.1" = "Setup C",
  "pF_exp_x1_x2.mF_max2_x1_x5.tF_max_x1_x5" = "Setup D")

# Elements of plot
scA <- function(which.given, ..., factor.levels) {
  if (which.given == 2) {
    strip.default(which.given = which.given, ..., factor.levels)
  } else {
    strip.default(which.given = 1, ...,
      factor.levels = exA[factor.levels])
  }
}

scB <- function(which.given, ..., factor.levels) {
  if (which.given == 2) {
    strip.default(which.given = which.given, ..., factor.levels)
  } else {
    strip.default(which.given = 1, ...,
      factor.levels = exB[factor.levels])
  }
}

mypanel <- function(x, y, groups, subscripts, ...) {
  fill <- cols[intersect(levels(x), unique(x))]
  panel.bwplot(x = x, y = y, fill = fill, ...)
  if (exists("use_diff") && use_diff) {
    panel.abline(h=0, lty = "dotted")
  }
  tapply(1:length(y), groups[subscripts], function(i) {
    xi <- 1:nlevels(x)
    yi <- y[i][match(xi, unclass(x[i]))]
    llines(x = xi, y = yi,
      col = rgb(0.1, 0.1, 0.1, 0.03))
  })
}

# Preprocess data
refactor <- function(df, method.nams) {
  elev <- c("pF_x1", "pF_x3", "pF_x4", "pF_eta_x1_x2", "0.5", "pF_x2_x3", "pF_exp_x1_x2")
  mlev <- c("0", "mF_x1", "mF_x3", "mF_sin_x1_x5", "mF_max_x1_x5", "mF_log_x1_x3", "mF_max2_x1_x5")
  tlev <- c("0", "tF_exp_x1_x2", "tF_div_x1_x2", "tF_log_x1_x2", "1", "tF_max_x1_x5")
  nlev <- c(800, 1600)
  dlev <- c(10, 20)

  df$p <- factor(df$p, labels = elev, levels = elev)
  df$m <- factor(df$m, labels = mlev, levels = mlev)
  df$t <- factor(df$t, labels = tlev, levels = tlev)
  df$nsim <- factor(df$nsim, labels = nlev, levels = nlev)
  levels(df$nsim) <- paste("N = ", levels(df$nsim))
  df$dim <- factor(df$dim, labels = dlev, levels = dlev)
  levels(df$dim) <- paste("P = ", levels(df$dim))

  df$i <- with(df, interaction(p, m, t))[, drop = TRUE]
  df$i <- factor(df$i, levels = names(exA))
  df$i <- droplevels(df$i)
  df$nd <- with(df, interaction(dim, nsim, sep = ", "))
  lev <- levels(df$nd)
  levels(df$nd) <- sapply(strsplit(lev, ", "), function(x) paste(x[2], ", ", x[1]))

  df <- df[, c("problem", "algorithm", "repl", "i", "nd", "ol", "seed", "result.res")]

  df$algorithm <- factor(as.character(df$algorithm), levels = method.nams,
    labels = gsub("honest", "-honest", method.nams))
  df <- df[!is.na(df$algorithm),]
  names(df)[names(df) == "result.res"] <- "value"
  return(df)
}

# Names of methods to show
methodnams <- c(
  "mob", "hybrid",
  "equalized", "mobcf",
  "cf"
  )
colornams <- lookup[methodnams]

# also show honest versions?
if (HONESTY) {
  methodnams <- c(
    "mob",  "hybrid", "equalized", "mobcf",
    "cf",
    "mobhonest", "hybridhonest", "equalizedhonest", "mobcfhonest",
    "cfhonest"
    )

  colornams <- lookup[names(cols)]
}

Nna <- sum(is.na(res$result.res))
res <- refactor(res, method.nams =  methodnams)

ylab <- expression(paste(frac(1, 1000), "",
  sum((tau(x[i]) - hat(tau)(x[i]))^2, i == 1, 1000)))

# If use_diff = TRUE, display MSE difference to cf: MSE(...) - MSE(cf)
if (exists("use_diff") && use_diff) {
  if (HONESTY) {
    res$honesty <- grepl("honest", res$algorithm)
    res$algo <- gsub("-honest", "", res$algorithm)
    res <- res %>%
      group_by(problem, repl, i, nd, ol, seed, honesty) %>%
      arrange(algo) %>%
      mutate(value = `value` - `value`[algo == "cf"]) %>%
      ungroup()
    res$algo <- NULL
    res$honesty <- NULL
    ylab <- "MSE(...) - MSE((h-)cf)"
  } else {
    res <- res %>%
      group_by(problem, repl, i, nd, ol, seed) %>%
      arrange(algorithm) %>%
      mutate(value = `value` - `value`[algorithm == "cf"])
    ylab <- "MSE(...) - MSE(cf)"
  }
  stopifnot(all(res$value[res$algorithm == "cf"] == 0))
  stopifnot(all(res$value[res$algorithm == "cf-honest"] == 0))
	res <- res[!res$algorithm %in% c("cf", "cf-honest"), ]
	res$algorithm <- droplevels(res$algorithm)
}

### w/o overlap: W
normalA <- subset(res, ol == 0)
### w overlap: W - .5
normalB <- subset(res, ol == 0.5)

## Omit color coding in figures for methods not shown
colsub <- cols[as.character(levels(res$algorithm))]
if (HONESTY) colsub <- cols[c("mob",  "mob-honest", "hybrid",  "hybrid-honest",
  "equalized", "equalized-honest", "mobcf", "mobcf-honest", "cf",  "cf-honest")]

nrcols <- 5L

if (exists("use_diff") && use_diff) {
  if (!HONESTY) {
    colsub <- colsub[-5]
    colsub <- colsub[-5]
    colornams <- colornams[-5]
    nrcols <- 4L
  } else {
    colsub <- colsub[c(-9, -10)]
    colornams <- colornams[c(-9, -10)]
    lookup <- lookup[c(-5, -10)]
    nrcols <- 4L
  }
}

mykey <- list(space="top", rectangles=list(col=colsub),
  text=list(colornams), columns = nrcols)

### plots
plot_results <- function(pltdata, sc, ylim, cexstrip = 0.5, rot = 60) {
  plt <- bwplot(value ~ algorithm | i + nd, data = pltdata,
    ylab = list(ylab), ylim = ylim,
    groups = repl, panel = mypanel, layout = c(4, 4),
    as.table = TRUE, strip = sc,
    key = mykey,
    scales = list(relation = "same", x = list(rot = rot,
      labels = lookup)),
    par.strip.text = list(cex = cexstrip),
    par.settings = list(layout.heights = list(strip = 1.5)))
  useOuterStrips(plt, strip = sc)
}

### tables
frmt2 <- function(x, math = FALSE) {
  if (!is.numeric(x)) return(x)
  ret <- formatC(round(x, 3), digits = 3, format = "f")
  if (math) ret <- paste("$", ret, "$")
  if (is.matrix(x)) {
    ret <- matrix(ret, nrow = nrow(x))
    dimnames(ret) <- dimnames(x)
  }
  ret
}

f <- function(x) {

  prm <- do.call("rbind", strsplit(rownames(x), ":"))
  prm <- cbind(prm[,1], do.call("rbind", strsplit(prm[,2], ", ")))

  txt <- paste(frmt2(x[, "Estimate"]),
    " (", frmt2(x[, "lwr"]), ", ", frmt2(x[, "upr"]), ")", sep = "")
  sup <- x[, "lwr"] > 1
  inf <- x[, "upr"] < 1
  txt[sup] <- paste("\\textbf{", txt[sup], "}", sep = "")
  txt[inf] <- paste("\\textit{", txt[inf], "}", sep = "")

  prm <- cbind(prm, txt)
  prm
}

calculate_ci <- function(data, lev) {

  data$ID <- with(data, interaction(i, nd))
  tmp <- glm(value ~ 0 + i : nd : algorithm,
             data = data, family = gaussian(link = "log"))
  mA <- glmmTMB(value ~ 0 + i : nd : algorithm + (1 | ID:repl),
                data = data, family = gaussian(link = "log"), start = list(beta = coef(tmp)))
  ### for multcomp::glht
  class(mA) <- c("my", class(mA))
  coef.my <- function(object) fixef(object)$cond
  vcov.my <- function(object) {
      class(object) <- class(object)[-1L]
      vcov(object)$cond
  }

  cf <- coef(mA)
  KA <- diag(length(cf))
  rownames(KA) <- names(cf)
  # Research Question 1:
  Kcf_mob <- KA[grep("algorithmcf$", names(cf)),] - KA[grep("algorithmmob$", names(cf)),]
  rownames(Kcf_mob) <- gsub(":algorithmcf", "", rownames(Kcf_mob))

  # Research Question 2:
  Kmc_eq <- KA[grep("algorithmmobcf$", names(cf)),] - KA[grep("algorithmequalized$", names(cf)),]

  # Research Question 3:
  Kcf_mc <- KA[grep("algorithmcf$", names(cf)),] - KA[grep("algorithmmobcf$", names(cf)),]

  # Research Question 4:
  Khb_mob <- KA[grep("algorithmhybrid$", names(cf)),] - KA[grep("algorithmmob$", names(cf)),]

  # Research Question 5:
  Khb_mc <- KA[grep("algorithmhybrid$", names(cf)),] - KA[grep("algorithmmobcf$", names(cf)),]

  # Research Question 5:
  Khb_eq <-  KA[grep("algorithmhybrid$", names(cf)),] - KA[grep("algorithmequalized$", names(cf)),]

  ci <- confint(glht(mA, linfct = rbind(Kcf_mob, Kmc_eq, Kcf_mc, Khb_mob, Khb_mc, Khb_eq)))
  ci1 <- exp(ci$confint[1:nrow(Kcf_mob),])
  ci2 <- exp(ci$confint[1:nrow(Kmc_eq) + nrow(Kcf_mob),])
  ci3 <- exp(ci$confint[1:nrow(Kcf_mc) + nrow(Kcf_mob) + nrow(Kmc_eq),])
  ci4 <- exp(ci$confint[1:nrow(Khb_mob) + nrow(Kcf_mob) + nrow(Kmc_eq) + nrow(Kcf_mc),])
  ci5 <- exp(ci$confint[1:nrow(Khb_mc) + nrow(Kcf_mob) + nrow(Kmc_eq) + nrow(Kcf_mc) + nrow(Khb_mob),])
  ci6 <- exp(ci$confint[1:nrow(Khb_eq) + nrow(Khb_mc) + nrow(Kcf_mob) + nrow(Kmc_eq) + nrow(Kcf_mc) + nrow(Khb_mob),])

  colnams <- c("cf:mob", "mobcf:equalized", "cf:mobcf", "hybrid:mob", "hybrid:mobcf", "hybrid:equalized")

  tabA <- cbind(f(ci1), f(ci2)[,4], f(ci3)[,4], f(ci4)[,4], f(ci5)[,4], f(ci6)[,4])

  colnames(tabA) <- c("LP", "N", "P", colnams)
  tabA[, "N"] <- gsub("ndN =  | ", "", tabA[, "N"])
  tabA[, "P"] <- gsub(" P =  ", "", tabA[, "P"])
  tabA <- as.data.frame(tabA)
  tabA$LP <- gsub("^i", "", tabA$LP)
  tabA$LP <- factor(tabA$LP, lev = names(lev),
    labels = lev)
  tabA$N <- factor(tabA$N, levels = c(800, 1600), labels = c(800, 1600))
  tabA$P <- factor(tabA$P, levels = c(10, 20), labels = c(10, 20))

  ndA <- expand.grid(P = unique(tabA$P), N = unique(tabA$N), LP = unique(tabA$LP))[,3:1]

  ndA <- join(ndA, tabA, by = c("LP", "N", "P"))

  return(ndA)
}

create_table <- function(dataA, dataB = NULL, levA, levB, honesty = FALSE) {
  if (honesty) {
    dataA <- dataA[stringr::str_detect(dataA$algorithm, "-honest"),]
    dataA$algorithm <- stringr::str_remove_all(dataA$algorithm, "-honest")
    if (!is.null(dataB)) {
      dataB <- dataB[stringr::str_detect(dataB$algorithm, "-honest"),]
      dataB$algorithm <- stringr::str_remove_all(dataB$algorithm, "-honest")
    }
  }

  # Get results latex table
  ndA <- calculate_ci(dataA, levA)
  if (!is.null(dataB)) {
    ndB <- calculate_ci(dataB, levB)
    x <- rbind(ndA, ndB)
  } else {
    x <- ndA
  }

  # Combine both
  all(unique(x$LP) == levels(x$LP))

  # Create results table
  hypcolnams <- c("(RQ 1)", "(RQ 2)", "(RQ 3)", "(RQ 4)", "(RQ 5)", "(RQ 5)")
  tabcolnams <- c("\\shortstack{(RQ 1) \\\\ cf vs. mob \\textcolor{white}{$\\hat{0}$}}",
    "\\shortstack{(RQ 2) \\\\ mobcf vs. mob($\\hat{W}, \\hat{Y}$)}",
    "\\shortstack{(RQ 3) \\\\ cf vs. mobcf\\textcolor{white}{$\\hat{0}$}}",
    "\\shortstack{(RQ 4) \\\\ mob($\\hat{W}$) vs. mob}",
    "\\shortstack{(RQ 5) \\\\ mob($\\hat{W}$) vs. mobcf}",
    "\\shortstack{(RQ 5) \\\\ mob($\\hat{W}$) vs. mob($\\hat{W}, \\hat{Y}$)}")
  nrowtab <- length(levA) - 1

  ftab <- as.matrix(x[, 4:9])
  class(ftab) <- "ftable"
  if (!is.null(dataB)) {
    attributes(ftab)$row.vars <- list("\\rule{0pt}{6ex} Part \\hspace{.4cm} DGP" =
        paste(c("A \\hspace{.8cm}",
          replicate(nrowtab, "\\hspace{1.1cm}"),
          "\\hline B \\hspace{.8cm}",
          replicate(nrowtab, "\\hspace{1.1cm}")),
          levels(x$LP)),
      "N" = levels(x$N),
      "P" = levels(x$P))
  } else {
    attributes(ftab)$row.vars <- list("DGP" = levels(x$LP),
      "N" = levels(x$N),
      "P" = levels(x$P))
  }
  attributes(ftab)$col.vars <- list("Mean squared error ratio" = tabcolnams)
  rownames(ftab) <- NULL

  lt <- toLatex(ftab, useDcolum = FALSE, useBook = FALSE, toprule = "\\hline", bottomrule = "\\hline")
  class(lt) <- "Latex"
  return(lt)
}
