#' Blood Loss Dataset
#'
#' The blood loss dataset is based on a study by the University Hospital Zurich from October 2015 to
#' November 2016.
#' The outcome is defined as the measured postpartum blood loss `MBL` in mL.
#' Details are given in Haslinger et al. (2020).
#'
#' @name blood
#'
#' @format A data frame with 1309 and 35 variables:
#' \describe{
#'   \item{Hb.prae}{prepartum hemoglobin}
#'   \item{MBL}{measured blood loss}
#'   \item{EC}{erythrocyte concentrate}
#'   \item{HOSP.7d}{length of hospitalization after delivery}
#'   \item{COLLOIDS}{number of colloids}
#'   \item{SPONTANEOUS}{spontaneous delivery}
#'   \item{VE}{vacuum delivery}
#'   \item{SECTIO.prim}{elective cesarean delivery}
#'   \item{SECTIO.sek}{unplanned cesarean delivery}
#'   \item{SECTIO.not}{emergency cesarean delivery}
#'   \item{GA}{gestational age}
#'   \item{AGE}{maternal age}
#'   \item{MULTIPAR}{multiparity}
#'   \item{BMI}{body mass index}
#'   \item{DAUER.ap}{duration of second stage labor}
#'   \item{TWIN}{multiple fetus pregnancy}
#'   \item{IOL}{induction of labor}
#'   \item{IOL.48h}{induction of labor $>$ 48 hours}
#'   \item{AIS}{chorioamnionitis}
#'   \item{NW}{neonatal weight}
#'   \item{RUPTUR}{uterine rupture}
#'   \item{ATONY}{uterine atony}
#'   \item{PLAC.RET}{retained placenta}
#'   \item{PLAC.REST}{retained placental material}
#'   \item{PLAC.MAP}{morbidly adherent placenta}
#'   \item{PLAC.PRAE}{placenta praevia}
#'   \item{RISS}{bleeding from laceration}
#'   \item{PLAC.LSG}{placental abruption}
#'   \item{F2.prae}{prepartum F. II}
#'   \item{F1.prae}{prepartum F. I}
#'   \item{F13.Akt.prae}{prepartum F. XIII}
#'   \item{SUBST.F13}{substitution F. XIII}
#'   \item{SUBST.F1}{substitution F. I}
#'   \item{F13.Ag.prae}{prepartum F. XIII (Ag)}
#'   \item{DD.prae}{prepartum D-dimer}
#'   }
#' @references
#' Dandl S, Torsten H, Heidi S et al. (2021). Divide or Unite? Theoretical and Empirical Insights into
#' Two Strategies to Random Forest-type Heterogeneous Treatment Effect Estimation.
#'
#' Haslinger C, Korte W, Hothorn T, Brun R, Greenberg C, Zimmermann R (2020). “The impact
#' of prepartum factor XIII activity on postpartum blood loss.” Journal of Thrombosis and
#' Haemostasis, 18. doi:10.1111/jth.14795.
#'
#'
NULL




#' Study Setting and Results of Dandl et al. (2021)
#'
#' A dataset containing the study setting and performance results of 8 different
#' tree-based methods for estimating treatment effects based on Dandl et al. (2021).
#' Mean squared error of estimated to true treatment effect measured the performance.
#' Study settings were based on Wager and Athey (2018) and Nie and Wager (2020).
#'
#' @name dandl
#'
#' @format A data frame with 9600 rows and 23 variables:
#' \describe{
#'   \item{setup}{Identifier for data generating process defined by p, m, t, sd, pi, model, xmodel.}
#'   \item{repl}{Replication index (1 to 100).}
#'   \item{p, m, t}{Propensity score, prognostic effect and treatment effect
#'   function character names. Input for \code{dgp}.}
#'   \item{sd, ol, model, xmodel}{Standard deviation of
#'   normal distributed outcome, fraction of treatment effect added to prognostic
#'   effect and names of used model for outcome and covariates.
#'   Input for \code{dgp}.}
#'   \item{nsim, dim, nsimtest, seed}{Number of observations, number of
#'   covariates, number of observations for test set and seed. Input for
#'   \code{simulate}.}
#'   \item{mob, mobhonest, hybrid, hybridhonest, equalized, equalizedhonest,
#'   cf, cfhonest, mobcf, mobcfhonest}{Mean squared error of different
#'   methods for estimating treatment effects. See Dandl et al. (2021) for details.}
#' }
#' @references
#' Dandl S, Torsten H, Heidi S et al. (2021). Divide or Unite? A Comparison of Two
#' Strategies to Random Forest-type Heterogeneous Treatment Effect Estimation.
#'
#' Wager S, and Athey S (2018). "Estimation and Inference of Heterogeneous Treatment Effects using Random Forests". Journal of the American Statistical Association, 113(523).
#'
#' Nie X, Wager S (2020). “Quasi-Oracle Estimation of Heterogeneous Treatment Effects.” Biometrika. ISSN 0006-3444. doi:10.1093/biomet/asaa076. Asaa076, https://academic.oup.com/biomet/advance-article-pdf/doi/10.1093/biomet/asaa076/33788449/asaa076.pdf, URL https://doi.org/10.1093/biomet/asaa076.
#'
NULL

