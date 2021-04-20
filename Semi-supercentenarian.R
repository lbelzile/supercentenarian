#Set directory when sourcing file
wd_source <- try(setwd(utils::getSrcDirectory()[1]))
#Set directory from Rstudio to current document
wd_rstudio <- try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
if(!isTRUE(is.character(wd_source)) && !isTRUE(is.character(wd_rstudio))){
  warning("Could not set working directory")
}

source("00-setup.R")
# Functions for fitting the various model
source("Semi-supercentenarian_fn.R")
###################################################
## Analysis of the Italian semi-supercentenarian ##
###################################################
# This database must be purchased from Istat
load("italcent.rda")
# 105 years threshold (criterion for inclusion in dataset)
u <- 38351L
# Calendar date at which individual reaches 105 years
xcal <- italcent$birth + u
# Calendar time for sampling frame
c1 <- lubridate::dmy("01-01-2009")
c2 <- lubridate::dmy("01-01-2016")

# Lower truncation level, zero if individual reached 105 between c1 and c2
slow <- as.numeric(pmax(0, c1 - xcal))
# Figure 1: Lexis diagram
source("01-Lexis_diagram_italcent.R")
# MLE of Generalized Pareto (GP) and exponential for Istat
# Table 1 (top panel)
source("02-fit_italcent_gpd_exp.R")
# Figure 6: Exponential QQ-plot for Istat
source("03-exponential_QQ_plot.R")
# Additional analyses of SM, section E
# Power for detecting difference in survival time per gender
source("04-power_gender_italcent.R")
# Figure 2 (top panel): parameter stability plots
source("05-param_stability_italcent.R")
# Figure 5: parameter stability plots, by cohort
source("06-param_stability_cohort_italcent.R")
# Figure 4: local hazard, by cohort
source("07-local_hazard_cohort_italcent.R")
# LRT for difference per gender for Istat
# These are not quoted, but are mentioned in the Discussion
source("08-lrt_gender_gp_exp_italcent.R")
# Figure 7: local hazard of GP, fitted using splines
source("09-local_gp_spline_italcent.R")
# Table 3 (top panel) and 4: comparisons with the Gompertz model
source("10-gompertz_extendedgp_italcent.R")

##################################################
## Analysis of the French semi-supercentenarian ##
##################################################
# Data extracted October 2019
load("francent.rda")
u <- 38350L
xcal <- francent$birth + u
# Calendar time for sampling frame
c1 <- lubridate::dmy("01-01-1987")
c2 <- lubridate::dmy("31-12-2017")

# Lower truncation level, zero if individual reached 105 between c1 and c2
francent$slow <- as.numeric(pmax(0, c1 - xcal))
francent$supp <- as.numeric(c2 - xcal)


# MLE of Generalized Pareto (GP) and exponential for French
# Table 1 (bottom panel) + LRT test for gender (Appendix E)
source("11-fit_francent_gpd_exp.R")
# Parameter stability plots for France
# Figure 2 (bottom panel)
source("12-param_stability_francent.R")
# Gompertz fit to French data
# Table 3 (bottom panel)
source("13-gompertz_extendedgp_francent.R")

##################################################

# Figure 3: Power analysis for test of finite endpoint 
# and power study for test of exponential distribution 
source("14-power_shape.R")
source("15-power_endpoint.R")
source("16-power_plots.R")
# Table 2: exponential distribution
source("17-exponential.R")