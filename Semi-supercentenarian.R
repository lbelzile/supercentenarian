rm(list=ls())
# Two files are used, italcent.rda and francent.rda
# The first can be bought for a small fee from the
# National Institute of Statistics by registering at the Contact Center
# (https://contact.istat.it) and mentioning the Semisupercentenarian
# Survey and Marco Marsili as contact person.
# 
# The second data can be downloaded by registering on 
# http://www.supercentenarians.org/
# Load libraries and packages used in the analyses
library(boot)  
library(evd)
library(mev)
library(xts)
library(lubridate)
library(survival)
library(xtable)
library(tikzDevice)
# The following package is not available on the CRAN
# devtools::install_github("OpenIntroStat/openintro-r-package", subdir = "OIsurv")
library(OIsurv)
source(Supercentenarian_fn.R)
dwidth <- 4
dheight <- 2.5

#Set directory when sourcing file
# setwd(utils::getSrcDirectory()[1])
#Set directory from Rstudio to current document
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
fig_dir <- paste0(substr(getwd(),start = 1, stop = nchar(getwd())-4), "figure")
code_dir <- getwd()
table_dir <- paste0(substr(getwd(),start = 1, stop = nchar(getwd())-4), "tables")
# To save Figures and Tables (.tex format), change the following to TRUE
figures <- tables <- FALSE
# This database must be purchased from Istat
load("italcent.rda")

########################################################################
########################################################################
# Figure 1: Lexis diagram

# Split sample by sex and right-censoring status
menunc <- which(as.logical((as.numeric(italcent$gender)-1) * !italcent$rightcens))
mencens <- which(as.logical((as.numeric(italcent$gender)-1) * italcent$rightcens))
womunc <- which(as.logical(((as.numeric(italcent$gender)-1)+1) %%2 * !italcent$rightcens))
womcens <- which(as.logical(((as.numeric(italcent$gender)-1)+1) %%2 * italcent$rightcens))

if(figures){
  fig <- "Fig1.tex"
  setwd(fig_dir)
  tikz(fig, width = dwidth, height = dheight, standAlone = TRUE)
}
par(mar = c(4,4,1,3))
plot(x = (italcent$birth + italcent$numdays)[womunc], 
     y = italcent$numdays[womunc]/365.25,
     col = scales::alpha("black", alpha = 0.5),
     bty = "l", pch = 20, cex = 0.8, ylab = "Age at death (in years)", xaxt = 'n',
     xlab = "year",yaxs = "i", ylim = c(105, 117), xlim = c(14244, 16825), xaxs = "i")
axis(side=1,at=seq(14244, 16825, by = 365.25), labels = 2009:2016, tick = TRUE)
points(x = (italcent$birth + italcent$numdays)[menunc], 
       y = italcent$numdays[menunc]/365.25,  pch = 4, cex = 0.8,
       col = scales::alpha("red", alpha = 0.5),lwd = 2)
rug(italcent$numdays[mencens]/365.25, side = 4, line =  1, lwd = 0.04, col = scales::alpha(2, 0.5), ticksize = 0.05)
rug(italcent$numdays[womcens]/365.25, side = 4, line =  2, lwd = 0.04, col = scales::alpha(1, 0.5), ticksize = 0.05)
# Number of death per year
deathcohorts_wom <- table(substr((italcent$birth + italcent$numdays)[womunc], 1, 4))
deathcohorts_men <- table(substr((italcent$birth + italcent$numdays)[menunc], 1, 4))
# Add counts to top of figure
text(labels = paste0(deathcohorts_wom, "/\\textcolor{red}{", deathcohorts_men,"}"), 
     x =as.Date(paste0(2009:2015, "-06-15")), y = 116, cex = 0.8)
if(figures){
  dev.off()
  system(command = paste0("pdflatex ", getwd(), "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}
########################################################################
########################################################################
# Fit generalized Pareto distribution to the data for multiple threshold
########################################################################
########################################################################
# Threshold for 105 years (criterion for inclusion in dataset
u <- 38351L
# Calendar date at which individual reaches 105 years
xcal <- italcent$birth + u
# Calendar time for sampling frame
c1 <- lubridate::dmy("01-01-2009")
c2 <- lubridate::dmy("01-01-2016")

# Lower truncation level, zero if individual reached 105 between c1 and c2
slow <- as.numeric(pmax(0, c1 - xcal))
# Exceedances
dat <- italcent$numdays - u
rightcens <- italcent$rightcens
expo <- FALSE

# Sequence of thresholds, from 105 years (in days), increments of half years
thresh <- seq(from = 0, by = 365.25/2, length = 14L)

# Fit a generalized Pareto model
param_gpd <- matrix(0, nrow = length(thresh), ncol = 9L)
colnames(param_gpd) <- c("loc","scale","shape","deviance", "scale.stderror", "shape.stderror","endpoint", "endpoint.stderror", "nu")
for(i in 1:length(thresh)){
  #For each threshold, compute new threshold exceedances
  datu <- dat - thresh[i]
  # Keep only exceedances
  ind <- which(datu > 0)
  # Compute maximum likelihood estimates numerically
  vals <- optim(par = c(5, 0.01), 
                fn = gpd_cens, 
                method = "N", 
                dat = (datu[ind])/365.25,
                rightcens = rightcens[ind], 
                slow = (pmax(0, slow[ind]-thresh[i]))/365.25, 
                expo = FALSE, 
                control = list(fnscale = c(1,0.1), reltol = 1e-12, maxit= 1e5),
                hessian = TRUE)
  # Compute standard errors
  endpt_stderror <- sqrt(solve(obs.infomat(c(- vals$par[1]/vals$par[2], vals$par[2]), 
                                           dat = (datu[ind])/365.25,
                                           rightcens = rightcens[ind], 
                                           slow = (pmax(0, slow[ind]-thresh[i]))/365.25))[1,1])
  # Bundle estimates in a vector
  param_gpd[i,] <- c(loc = (u + thresh[i])/365.25, vals$par, -2*vals$value, sqrt(diag(solve(vals$hessian))),
                     ifelse(vals$par[2] >0, Inf, (u + thresh[i])/365.25 - vals$par[1]/vals$par[2]), endpt_stderror,  length(ind))
}


# Exponential model
param_exp <- matrix(0, nrow = length(thresh), ncol = 4L)
colnames(param_exp) <- c("loc","scale","deviance", "scale.stderror")
for(i in 1:length(thresh)){
  datu <- dat - thresh[i]
  ind <- which(datu > 0)
  vals <- optim(par = c(500), fn = gpd_cens, method = "Brent", lower = 0, upper = 800, dat = (datu[ind])/365.25,
                rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25, expo = TRUE, hessian = TRUE)
  param_exp[i,] <- c(u + thresh[i], vals$par, -2*vals$value, sqrt(solve(vals$hessian)))
}
# difference in deviance between nested models
param_gpd <- cbind(param_gpd, 
                   pchisq(param_gpd[,"deviance"] - param_exp[,"deviance"], 1, lower.tail = FALSE),
                   pnorm(sign(param_gpd[,"shape"])*sqrt(param_gpd[,'deviance'] - param_exp[,'deviance']))
)

if(tables){
  table.res <- matrix(nrow = 6, ncol = 7)
  ind <- seq(1L, 13L, by = 2L)
  table.res[1,] <- paste0("$", round(param_gpd[ind,"nu"],0),"$")
  table.res[2,] <- paste0("$", round(param_gpd[ind, "scale"],2), "\\; (",round(param_gpd[ind, "scale.stderror"],2),")$")
  table.res[3,] <- paste0("$", round(param_gpd[ind, "shape"],2), "\\; (",round(param_gpd[ind, "shape.stderror"],2),")$")
  table.res[4,] <- paste0("$", round(param_exp[ind, "scale"],2), "\\; (",round(param_exp[ind, "scale.stderror"],2),")$")
  table.res[5,] <- paste0("$", round(param_gpd[ind,10],2),"$")
  table.res[6,] <- paste0("$", round(param_gpd[ind,11],2),"$")
  rownames(table.res) <- c("$n_u$", "$\\sigma$", "$\\gamma$", "$\\sigma_e$", "$p$-value","$p_{\\infty}$")
  xtab1 <- xtable::xtable(table.res, 
                          caption = "Estimates (standard errors) of scale  and shape parameters ($\\sigma$, $\\gamma$) for the generalized Pareto distribution and of the scale parameter ($\\sigma_e$) for the exponential model for the Italian \\Istat{} data as a function of threshold, with number of threshold exceedances ($n_u$), $p$-value for the likelihood ratio test of $\\gamma=0$ and probability under the fitted model that the endpoint is infinite, corresponding to $\\Pr(\\gamma>0)$.",
                          align = paste0(rep("r",ncol(table.res)+1),collapse = ""))
  setwd(table_dir)
  addtorow <- list()
  addtorow$pos <- list(0)
  addtorow$command <-  paste0("threshold", paste0(" & $", 105:111, "$", collapse = " "), "\\\\\n", collapse = " ")
  print(xtab1, file = "Table1a.tex", 
        size = "footnotesize", 
        floating.environment = "table*",  
        booktabs = TRUE, caption.placement = "top",
        add.to.row = addtorow, include.colnames = FALSE,
        sanitize.colnames.function = identity, 
        sanitize.text.function = identity, 
        sanitize.rownames.function = identity)
  setwd(code_dir)
}


#########################################################################################
#########################################################################################
# Fig 3: Cumulative hazard (Nelson-Aalen) 
# Survival object
italcent_surv <- Surv(time = slow, time2 = italcent$numdays - u, event = !italcent$rightcens)
# Fit Nelson--Alen estimator
NelAal <- survfit(italcent_surv ~ 1, stype = 2, ctype = 1)
# Fit Kaplan-Meier estimator
KM <- survfit(italcent_surv ~ 1)
# Extract cumulative hazard and its variance, plus time points
cumhaz0 <- NelAal$cumhaz
cumhazvar0 <- NelAal$std.chaz^2
time0 <- (NelAal$time + u)/365.25

fitsurv <- summary(survfit(italcent_surv ~ 1))
h.sort.of <- fitsurv$n.event / fitsurv$n.risk
cumhaz <- rep(0, length.out = length(h.sort.of))
cumhazvar <- rep(0, length.out = length(h.sort.of))
for(i in 1:length(h.sort.of)){
  cumhaz[i] <- sum(h.sort.of[1:i])
  cumhazvar[i] <- sum(h.sort.of[1:i]^2)
}



if(figures){
  setwd(fig_dir)
  fig <- "Fig3.tex" 
  tikz(fig, width = dwidth ,height = dheight, standAlone = TRUE)
}

par(mar = c(4,5.5,0.1,0.1))
# Plot of the Nelson-Aalen estimator
plot(time0, cumhaz0, 
     xlab='Survival time (in years)', 
     ylab='Conditional cumulative\n hazard above 105',
     yaxt="n", ylim=range(cumhaz0),
     bty = "l", type = "s", lwd = 2,
     panel.first = {
       abline(a = -time0[id0]*0.69+cumhaz0[id0], b = 0.69,
                           col = 4, lwd = 2)
      })
axis(side=2, at=0:7, labels=105:112, yaxs = "i")
id0 <- min(which(time0 > 108))

bands <- cumhazbands(time = time0, cumhaz = cumhaz0, type = "ptwise", 
                     transfo= "a", sigma.sq = cumhazvar0, confLevel = 0.95, n = length(xord))
lines(bands$time, bands$lower, col = "grey50", lty = 2, lwd = 2)
lines(bands$time, bands$upper, col = "grey50", lty = 2, lwd = 2)
# Expected warning about data(hw.k05)
bands <- cumhazbands(time = time0, cumhaz = cumhaz0, type = "band", 
                     transfo= "a", sigma.sq = cumhazvar0, confLevel = 0.95, n = length(xord))
lines(bands$time, cummax(bands$lower), col = "grey50", lty = 3, lwd = 2)
lines(bands$time, bands$upper, col = "grey50", lty = 3, lwd = 2)

if(figures){
  dev.off()
  system(command = paste0("lualatex ",fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}

#########################################################################################
#########################################################################################
# Fig 6: Exponential quantile-quantile plot with parametric bootstrap

ind108 <- which((italcent$numdays - u) > 365*3)
slow108 <-  pmax(0, slow[ind108]-365*3)
numd108 <- italcent$numdays[ind108] - u - 365*3
uplim108 <- c2 - italcent$birth[ind108]
cens108 <- italcent$rightcens[ind108]
italcent_surv_108 <- Surv(time = slow108, 
                          time2 = numd108, 
                          event = !italcent$rightcens[ind108])

# Compute maximum likelihood estimator
param_exp108  <- exp_mle_lt_rc(dat = numd108, slow = slow108, rightcens = italcent$rightcens[ind108])

B <- 10000L
set.seed(1234)
bootsamp <- matrix(0, ncol = length(numd108), nrow = B)
for(i in 1:ncol(bootsamp)){
  if(cens108[i]){
    bootsamp[,i] <- numd108[i]
  } else{
    # Due to censoring, ind108ividuals that are not censored are doubly truncated
    bootsamp[,i] <- qexp(pexp(slow108[i], rate = 1/param_exp108) + runif(B)*(pexp(uplim108[i], rate = 1/param_exp108) - pexp(slow108[i], rate = 1/param_exp108)), rate = 1/param_exp108)
  }
}

qqpts <- qqptsb <- matrix(0, nrow = B, ncol = sum(!cens108))
thetab <- rep(0, length.out = B)
for(j in 1:B){
  thetab[j] <- exp_mle_lt_rc(dat = bootsamp[j,],
                     rightcens = cens108, 
                     slow = slow108)
  
  #Fit a Kaplan-Meier for left-truncated right-censored data
  KM.b <- survfit(Surv(time = slow108, time2 = bootsamp[j,], event = !cens108)~1)
  #Extract the timing and make into a empirical distribution function (weighted)
  cecdf <- mev:::.wecdf(KM.b$time, diff(c(0, 1-exp(-KM.b$cumhaz))))
  # Evaluate on a regular grid to get confidence interval
  qqpts[j,] <- sort(qexp(cecdf(seq(0, 2800, length = sum(!cens108))), rate = 1/thetab[j]))
  # Evaluate at simulated data
  qqptsb[j,] <- sort(qexp(cecdf(bootsamp[j,][!cens108]), rate = 1/thetab[j]))
}

KM <- survfit(Surv(time = slow108, time2 = numd108, event = !cens108)~1)
cecdf <- mev:::.wecdf(KM$time, diff(c(0, 1-KM$surv)))


if(figures){
  setwd(fig_dir)
  fig <- "Fig6.tex"
  tikz(fig, width = dwidth, height = dheight, standAlone = TRUE)
}
par(mar = c(4,4,0.4,0.4), mfrow = c(1,2), cex = 1, bty = "l")
plot(x = NA,
     pch = 20, bty = "l", xlab = "Theoretical quantiles", ylab = "Observed quantiles",
     panel.first = {abline(a=0, b=1)}, yaxs = "i", xaxs = "i",
     ylim = 108+c(0, 3000)/365.25, xlim = 108+c(0, 3000)/365.25)
for(i in 1:100){
  lines(y = sort(bootsamp[i,!cens108])/365.25+108, qqptsb[j,]/365.25+108, 
        lty = 2, col = scales::alpha(colour = "black", 0.2))
}
points(y = sort(numd108[!cens108])/365.25+108, sort(qexp(cecdf(numd108[!cens108]), rate = 1/param_exp108))/365.25+108, pch = 1, cex = 0.8)
points(y = sort(numd108[!cens108])/365.25+108, sort(qexp(cecdf(numd108[!cens108]), rate = 1/param_exp108))/365.25+108, pch = 20, cex = 0.2, col = "white")
plot(y = sort(numd108[!cens108])/365.25+108, sort(qexp(cecdf(numd108[!cens108]), rate = 1/param_exp108))/365.25+108,
     pch = 20, bty = "l", xlab = "Theoretical quantiles", ylab = "Observed quantiles",
     panel.first = {abline(a=0, b=1)}, yaxs = "i", xaxs = "i",
     ylim = 108+c(0, 3000)/365.25, xlim = 108+c(0, 3000)/365.25, cex = 0.8)
env <- boot::envelope(mat = qqpts) 
ss <- seq(0, 2800, length = sum(!cens108))/365.25+108
lines(y = ss, env$point[1,]/365.25+108, lty = 2, lwd =1.5)
lines(y = ss, env$point[2,]/365.25+108, lty = 2, lwd =1.5)
lines(y = ss, env$overall[1,]/365.25+108, lty = 3, lwd =1.5)
lines(y = ss, env$overall[2,]/365.25+108, lty = 3, lwd =1.5)

if(figures){
  dev.off()
  system(command = paste0("lualatex ",fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}


#########################################################################################
#########################################################################################
# Power study - comparison of men versus women (SI, Section D)

m108 <- which(italcent$gender[ind108] == "male")
w108 <- which(italcent$gender[ind108] != "male")
gend108 <- italcent$gender[ind108]   

r108 <- 1/exp_mle_lt_rc(dat = numd108, rightcens = cens108, slow = slow108)

lrt_fun108a <- function(dat){
  datt <- as.integer(pmin(uplim108 - u, dat))
  rcens108 <- I((uplim108 -u) < dat)
  sigmae <- sum(datt - slow108)/sum(!rcens108)
  sigmaem <- sum(datt[m108] - slow108[m108])/sum(!rcens108[m108])
  sigmaew <- sum(datt[w108] - slow108[w108])/sum(!rcens108[w108])
  llexp <- gpd_cens(par = sigmae, dat = datt, rightcens = rcens108, slow = slow108, expo = TRUE)
  llexpm <- gpd_cens(par = sigmaem, dat = datt[m108], rightcens = rcens108[m108], slow = slow108[m108], expo = TRUE)
  llexpw <- gpd_cens(par = sigmaew, dat = datt[w108], rightcens = rcens108[w108], slow = slow108[w108], expo = TRUE)
  2*(llexp - (llexpm + llexpw))
}


lrt_fun108b <- function(dat){
  ## Same as above, but without right-censoring
  datt <- dat
  rcens108 <- rep(FALSE, length(dat))
  sigmae <- sum(datt - slow108)/sum(!rcens108)
  sigmaem <- sum(datt[m108] - slow108[m108])/sum(!rcens108[m108])
  sigmaew <- sum(datt[w108] - slow108[w108])/sum(!rcens108[w108])
  llexp <- gpd_cens(par = sigmae, dat = datt, rightcens = rcens108, slow = slow108, expo = TRUE)
  llexpm <- gpd_cens(par = sigmaem, dat = datt[m108], rightcens = rcens108[m108], slow = slow108[m108], expo = TRUE)
  llexpw <- gpd_cens(par = sigmaew, dat = datt[w108], rightcens = rcens108[w108], slow = slow108[w108], expo = TRUE)
  2*(llexp - (llexpm + llexpw))
}
# Power study
# For each lambda in 1, 1.1, ...
# 1) simulate from exponential with rate \mu/\lambda (women) and \mu \lambda (men) above 108
# 2) fit GP by maximum likelihood for the pooled sample and compute likelihood for men/women separately
# 3) preliminary check: the asymptotic null distribution is likely adequate
# so bootstrap p-value under lambda=1 should match that,
# otherwise our bootstrap generating process is not good enough... 
## (we have a p-value of 0.897 based on exponential at 108 chisq(1)
## and 0.899 for GP chisq(2)
## 4) compute power (but store LRT values)
B <- 10000L
set.seed(1234)
lam_vals <- seq(1, 1.5, length = 50L)
lrt_stata <- matrix(0, ncol = length(lam_vals), nrow = B)
lrt_statb <- matrix(0, ncol = length(lam_vals), nrow = B)
for(lam_ind in 1:length(lam_vals)){
  #Matrix to store the bootstrap sample
  bootsamp <- matrix(0, ncol = length(numd108), nrow = B)
  for(i in 1:ncol(bootsamp)){
    if(gend108[i]=="male"){
      rate <- r108*lam_vals[lam_ind]
    } else{
      rate <- r108/lam_vals[lam_ind]
    }
    # Do not censor individuals
    bootsamp[,i] <- ceiling(qexp(pexp(slow108[i], rate = rate) + runif(B)*(1 - pexp(slow108[i], rate = rate)), rate = rate))
  }
  # For each bootstrap sample
  # compute MLE for men, women and pooled sample
  lrt_stata[, lam_ind] <- apply(bootsamp, 1, lrt_fun108a)
  lrt_statb[, lam_ind] <- apply(bootsamp, 1, lrt_fun108b)
}

par(mar = c(4,4,1,1), pch = 20, bty = "l")
plot(colMeans(lrt_stata > qchisq(0.95,1)) ~ lam_vals,
     ylab = "power", xlab = expression(lambda))
# save(lam_vals, lrt_stat,lrt_statb, file = "power.RData")
y <- lam_vals[6:40]
x <- colMeans(lrt_stata[,6:40] > qchisq(0.95,1))
m <- mgcv::gam(y ~ s(x))
plot(m)
plot(y~x)
lines(fitted(m) ~ x)
powextra <- mgcv::predict.gam(m, newdata = data.frame(x = c(0.2,0.4,0.5)))
round((powextra/r108 - 1/(r108*powextra)),0)
par(mar = c(4,4,1,1), bty = "l")
plot(NULL, NULL,      ylab = "power", xlab = "$\\lambda$",yaxs = "i",xaxs = "i",ylim = c(0,1), xlim = c(1,2.25),
     panel.first = {abline(h = seq(0.1, 1, by = 0.1), col = scales::alpha(1, 0.05), lty = 1, lwd = 0.25)})
gam1 <- mgcv::gam(colMeans(lrt_stata > qchisq(0.95,1)) ~ s(lam_vals^2, k = -1, bs = "tp"))
lines(lam_vals^2, fitted(gam1), lwd = 2)
gam2 <- mgcv::gam(colMeans(lrt_statb > qchisq(0.95,1)) ~ s(lam_vals^2, k = -1, bs = "tp"))
lines(lam_vals^2, fitted(gam2), lwd = 2, lty = 3,  col = "grey30")
legend("bottomright", legend = c("censoring","no censoring"), bty = "n",
       col = c(1, "grey30"), lwd = 2, lty = c(1,3))
setwd(code_dir)


#########################################################################################
#########################################################################################
# Figure 2: parameter stability plots for Istat data

# Profile likelihood for shape of GP and scale of exponential distribution
# used to compute the parameter stability plots
prof_xizero <- rep(0, length.out = length(thresh))
confint_xi <- matrix(0, ncol = 3, nrow = length(thresh))
confint_exp <- matrix(0, ncol = 3, nrow = length(thresh))
for(i in 1:length(thresh)){
  datu <- dat - thresh[i]
  ind <- which(datu > 0)
  confint_xi[i,] <- prof_gpd_cens_xi_confint(dat = (datu[ind])/365.25,
                                             rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25)
  confint_exp[i,] <- prof_exp_cens(dat = (datu[ind])/365.25,
                                   rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25)
  prof_xizero[i] <- nllhxizero(dat = (datu[ind])/365.25, rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25)
  
}

# Compute the probability that xi >= 0 based on profile and likelihood root
probgamma0 <- pnorm(sign(param_gpd[,'shape'])*sqrt(param_gpd[,'deviance'] - prof_xizero))

if(figures){
  fig <- "Fig2.tex"
  setwd(fig_dir)
  tikz(fig, width = dwidth, height = 0.8*dheight, standAlone = TRUE)
}
par(mar = c(4, 4, 0.1, 0.1), mfrow = c(1,2))
plot(x = thresh/365.25+105, confint_xi[,1], type= "p", panel.first = abline(h = 0, col = "grey"), 
     pch = 20, col = 1, bty = "l", ylab = "$\\gamma$", xlab = "threshold (in years)",ylim = c(-0.3,1.1))
for(i in 1:length(thresh)){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_xi[i,2], y1 = confint_xi[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
plot(x = thresh/365.25+105, confint_exp[,1], type= "p", 
     panel.first = abline(h=1.45142127, col = "grey"), pch = 20, col = 1,bty = "l", ylab = "$\\sigma_e$", xlab = "threshold (in years)",ylim = range(confint_exp))
for(i in 1:length(thresh)){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_exp[i,2], y1 = confint_exp[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}

if(figures){
  dev.off()
  system(command = paste0("lualatex ", fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}




#########################################################################################
#########################################################################################
# Fig. 4 Parameter stability plot, by cohort 1896-1905 vs 1906-1910

confint_xi1 <- matrix(0, ncol = 3, nrow = length(thresh))
confint_exp1 <- matrix(0, ncol = 3, nrow = length(thresh))
confint_xi2 <- matrix(0, ncol = 3, nrow = 8)
confint_exp2 <- matrix(0, ncol = 3, nrow = 8)
for(i in 1:length(thresh)){
  datu <- dat - thresh[i]
  indi <- intersect(which(datu > 0), which(italcent$birth < as_date("1906-01-01"))) 
  confint_xi1[i,] <- prof_gpd_cens_xi_confint(dat = (datu[indi])/365.25,
                                              rightcens = rightcens[indi], slow = (pmax(0, slow[indi]-thresh[i]))/365.25)
  confint_exp1[i,] <- prof_exp_cens(dat = (datu[indi])/365.25,
                                    rightcens = rightcens[indi], slow = (pmax(0, slow[indi]-thresh[i]))/365.25)
  indi <- intersect(which(datu > 0), which(italcent$birth >= as_date("1906-01-01"))) 
  if(length(indi) > 20 && i < 9){
    confint_xi2[i,] <- prof_gpd_cens_xi_confint(dat = (datu[indi])/365.25,
                                                rightcens = rightcens[indi], slow = (pmax(0, slow[indi]-thresh[i]))/365.25)
    confint_exp2[i,] <- prof_exp_cens(dat = (datu[indi])/365.25,
                                      rightcens = rightcens[indi], slow = (pmax(0, slow[indi]-thresh[i]))/365.25)
  }
}

if(figures){
  fig <- "Fig4.tex"
  setwd(fig_dir)
  tikz(fig, width = dwidth, height = 1.6*dheight, standAlone = TRUE)
}
par(mar = c(4.5, 4.5, 0.1, 0.1), mfrow = c(2,2))
plot(x = thresh/365.25+105, confint_xi1[,1], type= "p", panel.first = abline(h = 0, col = "grey"), 
     pch = 19, col = 1, bty = "l", ylab = "$\\gamma$", xlab = "threshold (in years)",ylim =  c(-0.5,1))
for(i in 1:length(thresh)){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_xi1[i,2], y1 = confint_xi1[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
plot(x = thresh/365.25+105, confint_exp1[,1], type= "p", 
     panel.first = abline(h=1.45142127, col = "grey"), pch = 19, col = 1,
     bty = "l", ylab = "$\\sigma_e$", xlab = "threshold (in years)",ylim = range(confint_exp1))
for(i in 1:length(thresh)){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_exp1[i,2], y1 = confint_exp1[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}

plot(x = thresh[1:8]/365.25+105, confint_xi2[,1], type= "p", 
     panel.first = abline(h = 0, col = "grey"), 
     pch = 19, col = 1, bty = "l", ylab = "$\\gamma$", xlab = "threshold (in years)",ylim = c(-0.5,1))
for(i in 1:8){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_xi2[i,2], y1 = confint_xi2[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
plot(x = thresh[1:8]/365.25+105, confint_exp2[1:8,1], type= "p", 
     panel.first = abline(h=1.45142127, col = "grey"), pch = 19, col = 1,
     bty = "l", ylab = "$\\sigma_e$", xlab = "threshold (in years)",ylim = range(confint_exp2))
for(i in 1:8){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_exp2[i,2], y1 = confint_exp2[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
dev.off()
system(command = paste0("lualatex ", fig_dir, "/", fig,"; rm *.aux; rm *.log"))
setwd(code_dir)




#########################################################################################
#########################################################################################
# Fig. 7 Local hazard plots by birth cohort

# Local hazard using all data (not shown)
localhazard <- prof_gpd_hazard_confint(dat = dat, rightcens = rightcens, slow = slow, thresh = param_gpd[1,1]+0:6)

ind_l1 <-  which(italcent$birth < as_date("1906-01-01"))
localhazard1 <- prof_gpd_hazard_confint(dat = dat[ind_l1], rightcens = rightcens[ind_l1], slow = slow[ind_l1], thresh = 105:111)
ind_l2 <-  which(italcent$birth >= as_date("1906-01-01"))
localhazard2 <- prof_gpd_hazard_confint(dat = dat[ind_l2], rightcens = rightcens[ind_l2], slow = slow[ind_l2], thresh = 105:109)

if(figures){
  fig <- "Fig7.tex"
  setwd(fig_dir)
  tikz(fig, width = dwidth, height = dheight, standAlone = TRUE)
}
# Consider two thresholds for which xi is negative
yran <- range(c(localhazard1[,1:3], localhazard2[,1:3]))
par(mar = c(4.5, 4.5, 0.1, 0.1), mfrow = c(1,2))
plot(105:111, localhazard1[,1],  ylim = yran, panel.first = abline(h = 0.7, col = "grey"), 
     type= "p", pch = 19, bty = 'l', ylab = "local hazard", xlab = "age (in years)")
for(i in 1:nrow(localhazard1)){
  arrows(x0 = (105:111)[i], y0 = localhazard1[i,2], y1 = localhazard1[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
plot(105:109, localhazard2[,1],  ylim = yran, 
     type= "p", pch = 19, bty = 'l', ylab = "local hazard", 
     panel.first = abline(h = 0.7, col = "grey"), xlab = "age (in years)")
for(i in 1:nrow(localhazard2)){
  arrows(x0 = (105:109)[i], y0 = localhazard2[i,2], y1 = localhazard2[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}

if(figures){
  dev.off()
  system(command = paste0("lualatex ", fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}

#########################################################################################
#########################################################################################
# Testing difference between sex - perform likelihood ratio tests (SI, Section D)

# Generalized Pareto model - men only
param_gpd_men <- matrix(0, nrow = 8, ncol = 9L)
colnames(param_gpd_men) <- c("loc","scale","shape","deviance", "scale.stderror", "shape.stderror","endpoint", "endpoint.stderror", "nu")
for(i in 1:8){ # not enough males (12) beyond this stage to estimate parameters of the GP
  datu <- dat - thresh[i]
  ind <- as.logical(I(datu > 0) * I(italcent$gender=="male"))
  ind <- which(ind)
  vals <- optim(par = c(1.2,-0.1), fn = gpd_cens, method = "N", dat = (datu[ind])/365.25,
                rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25, expo = FALSE, hessian = TRUE)
  endpt_stderror <- sqrt(solve(obs.infomat(c(- vals$par[1]/vals$par[2], vals$par[2]), dat = (datu[ind])/365.25,rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25))[1,1])
  param_gpd_men[i,] <- c(loc = (u + thresh[i])/365.25, vals$par, -2*vals$value, sqrt(diag(solve(vals$hessian))),
                         ifelse(vals$par[2] >0, Inf, (u + thresh[i])/365.25 - vals$par[1]/vals$par[2]), endpt_stderror,  length(ind))
}

# Generalized Pareto model - women only
param_gpd_women <- matrix(0, nrow = 8, ncol = 9L)
colnames(param_gpd_women) <- c("loc","scale","shape","deviance", "scale.stderror", "shape.stderror","endpoint", "endpoint.stderror", "nu")
for(i in 1:8){
  datu <- dat - thresh[i]
  ind <- as.logical(I(datu > 0) * I(italcent$gender=="female"))
  vals <- optim(par = c(5,-0.01), fn = gpd_cens, method = "N", dat = (datu[ind])/365.25,
                rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25, expo = FALSE, hessian = TRUE)
  endpt_stderror <- sqrt(solve(obs.infomat(c(- vals$par[1]/vals$par[2], vals$par[2]), dat = (datu[ind])/365.25,rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25))[1,1])
  param_gpd_women[i,] <- c(loc = (u + thresh[i])/365.25, vals$par, -2*vals$value, sqrt(diag(solve(vals$hessian))),
                           ifelse(vals$par[2] >0, Inf, (u + thresh[i])/365.25 - vals$par[1]/vals$par[2]), endpt_stderror,  length(ind))
}
#P-values for the test: GP for both men and women, with different sigma_i and xi_i vs common scale+shape
1-pchisq((param_gpd_women[1:8,"deviance"] + param_gpd_men[1:8,"deviance"]) - param_gpd[1:8,"deviance"], df=2)

# Repeat the test, this times with exponential distribution
param_exp_women <- matrix(0, nrow = 9, ncol = 4L)
colnames(param_exp_women) <- c("loc","scale","deviance", "scale.stderror")
for(i in 1:9){
  datu <- dat - thresh[i]
  ind <- ind <- as.logical(I(datu > 0) * I(italcent$gender=="female"))
  vals <- optim(par = c(500), fn = gpd_cens, method = "Brent", lower = 0, upper = 800, dat = (datu[ind])/365.25,
                rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25, expo = TRUE, hessian = TRUE)
  param_exp_women[i,] <- c(u + thresh[i], vals$par, -2*vals$value, sqrt(solve(vals$hessian)))
}

param_exp_men <- matrix(0, nrow = 9, ncol = 4L)
colnames(param_exp_men) <- c("loc","scale","deviance", "scale.stderror")
for(i in 1:9){
  datu <- dat - thresh[i]
  ind <- as.logical(I(datu > 0) * I(italcent$gender=="male"))
  vals <- optim(par = c(500), fn = gpd_cens, method = "Brent", lower = 0, upper = 800, dat = (datu[ind])/365.25,
                rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25, expo = TRUE, hessian = TRUE)
  param_exp_men[i,] <- c(u + thresh[i], vals$par, -2*vals$value, sqrt(solve(vals$hessian)))
}
# P-values
1-pchisq((param_exp_women[1:9,"deviance"] + param_exp_men[1:9,"deviance"]) - param_exp[1:9,"deviance"], df=1)

#########################################################################################
#########################################################################################
# Bayesian analysis using HMC with vague priors (SI, Section F)

# Only run the code if the file with the results is not provided
if("Table3_draws.RData" %in% list.files()){
  run <- FALSE
}
# transform data to yearly scale
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
thresholdseq <- seq(from = 0, by = 365.25, length = 7L)

# Test now model with a different threshold
if(run){
for(j in 1:length(thresholdseq)){
  datuy <- (italcent$numdays - u - thresholdseq[j])/365.25
  slowy <- pmax(0,  as.numeric(pmax(0, c1 - xcal)) - thresholdseq[j])/365.25
  ind <- which(datuy > 0)
  assign(x = paste0("fit", j), value = stan(file = 'GPD_cens.stan',
                                              data = list(n= length(datuy[ind]), x = datuy[ind], rightcens = rightcens[ind], slow = slowy[ind], xmax = max(datuy[ind])),
                                              warmup = 1000,
                                              iter = 26000,
                                              thin = 5,
                                              chains = 5,
                                              init = c(1.5,0.01)))
  }
  save(list = ls(pattern="fit*"), file = "Table3_draws.RData")
  summary(fit5)
  traceplot(fit5)
  pairs(fit5)
} else{
  load("Table3_draws.RData")
}

interv <- matrix(0, nrow = length(thresholdseq), ncol = 5)
pgammapos <- rep(0, length.out = length(thresholdseq))
for(j in 1:length(thresholdseq)){
  # tau is endpoint
  tau <- (u + thresholdseq[j])/365.25 + ifelse(extract(get(paste0("fit",j)), "xi")$xi < 0,
                                          -extract(get(paste0("fit",j)), "sigma")$sigma/extract(get(paste0("fit",j)), "xi")$xi, Inf)
  interv[j,] <- quantile(tau, c(0.025, 0.05, 0.5, 0.95, 0.975))
  pgammapos[j] <- mean(extract(get(paste0("fit",j)), "xi")$xi >= 0)
}

medi <- format(round(interv[,3], 1), nsmall=1)
lowcr <- format(round(interv[,2], 1), nsmall=1)
uppcr <- format(round(interv[,4], 1), nsmall=1)
medi[which(medi == "  Inf")] <- "\\hphantom{00}\\infty\\hphantom{.}"
lowcr[which(lowcr == "  Inf")] <- "\\hphantom{00}\\infty\\hphantom{.}"
uppcr[which(uppcr == "  Inf")] <- "\\hphantom{00}\\infty\\hphantom{.}"
tab <- cbind(paste0("$", round((u+thresholdseq)/365.25,1),"$"), paste0("$", medi, "$ ($", lowcr, ", ", uppcr ,"$)"), paste0("$", round(pgammapos, 2),"$"))
colnames(tab) <- c("$u$", "endpoint $\\tau$ (90\\% CI)", "$\\Pr(\\gamma \\geq 0)$")
xtab <- xtable(tab, caption = c("Posterior median (90\\% credible intervals) for the upper limit to lifetime $\\tau$ (in years) and posterior probability of infinite endpoint, both as a function of the threshold $u$."),
               label = "tab_italcent_credint", align = "llrl")
if(tables){
  setwd(table_dir)
  print(xtab, 
        right = TRUE, 
        booktabs = TRUE, 
        include.rownames=FALSE, 
        sanitize.text.function = identity,
        sanitize.colnames.function = function(x){paste0("\\multicolumn{1}{c}{",x,"}")},
        sanitize.rownames.function = identity, 
        table.placement = "t!",
        caption.placement = "top",
        file = "Table3.tex")
  setwd(code_dir)
}

#########################################################################################
#########################################################################################
# Fig. 8 Local GP fit using splines with bootstrap analysis of variability
# Code by Anthony Davison
# Exceedances, with d=1 indicating full observation, d=0 indicating censoring
#  x = event time, s = start of observation time, t= censoring/death time, all in days
high <- as.numeric(c2 - xcal)
data <- data.frame(x=dat, t=high, s=slow, d=1-rightcens, yob=lubridate::year(italcent$birth), gender=italcent$gender)


# Use ordinary polynomial splines, with corresponding X matrix made with function make.X
make.X <- function(len=15, n.knots=5, knot.spacing=365, noise=0){
  # X matrix for discrete GPD
  x <- c(1:(len*365))/365
  one <- rep(1,len*365)
  X <- cbind(one,x)
  if (n.knots>0) for (i in 1:n.knots)
  { knot <- round( i*knot.spacing+noise*rnorm(1) )
  z <- pmax(x[knot]-x,0)
  X <- cbind(X, z^3) }
  list(X=X, nrows=nrow(X), ncols=ncol(X))
}

# X matrix for exponential model 
make.X.exp <- function(X){ 
  Y <- X 
  Y$X <- Y$X[,-2]
  Y$ncols <- Y$ncols - 1
  Y$X <- matrix(Y$X, nrow=X$nrows)
  Y
}

make.rh <- function(lam, X) list( rh=(X$X %*% lam), x=c(1:nrow(X$X)) ) #
# make discrete hazard function for given parameters and X matrix
make.hazard <- function(th, X){  
  rh <- make.rh(th, X)
  haz <- 1/pmax(rh$rh,0)
  Haz <- cumsum(haz)/365
  surv <- exp(-Haz)
  surv[is.nan(surv)] <- 0
  dens <- haz*surv
  list(x=rh$x/365, dens=dens, surv=surv, y=dens/surv)
}

nlogL.rh <- function(th, X, d){
  m <- make.hazard(th, X)
  denom <- m$surv[1+d$s]
  numer <- d$d*m$dens[d$x] + (1-d$d)*m$surv[d$x]
  ts <- Inf
  if (all(!is.nan(c(numer,denom)))) ts <- -sum(log(numer/denom))
  ts
}



boot.stat <- function(data, i, X){
  d <- data[i,]
  fit.rh <- optim(init, nlogL.rh, hessian=T, method="BFGS", X=X, d=d)
  c(fit.rh$par,fit.rh$value,fit.rh$convergence)
}

boot.stat1 <- function(data, i, n.knots=5, knot.spacing=365, noise=0){
  d <- data[i,]  # bootstrap dataset
  
  X <- make.X(n.knots=n.knots, knot.spacing=knot.spacing, noise=noise) # matrix X with random knots
  X.exp <- make.X.exp(X)  # cut-down matrix for exponential fit
  init <- c(1.6,rep(0,n.knots)) # initial values for optimisation
  
  fit.rh.exp <- optim(init, nlogL.rh, method="BFGS", X=X.exp, d=d) # exp optimisation
  h.exp <- make.hazard(fit.rh.exp$par, X.exp)
  
  init2 <- fit.rh.exp$par
  init2 <- c(init2[1],0,init2[-1])
  fit.rh <- optim(init2, nlogL.rh, method="Nelder-Mead", X=X, d=d) # GP optimisation
  fit.rh <- optim(fit.rh$par, nlogL.rh, method="BFGS", X=X, d=d) # GP optimisation
  h <- make.hazard(fit.rh$par, X)
  c(fit.rh$par,fit.rh$value,fit.rh$convergence, fit.rh.exp$par,fit.rh.exp$value,fit.rh.exp$convergence, length(h$x), h$x, h$y, h.exp$y)
}

# Warning: computationally intensive! Reduce number of bootstrap replicates to speed up calculations
# the paper used R=5000
system.time( rh.boot <- boot(data=data, statistic=boot.stat1, R=500L, 
                             n.knots=5, knot.spacing=365, noise=60))
# indices for plotting output
n.knots <- 5
N <- rh.boot$t0[2*n.knots+8]
ind.x <- 2*n.knots+8 + c(1:N)
ind.h <- 2*n.knots+8 + N + c(1:N)
ind.h.exp <- 2*n.knots+8 + 2*N + c(1:N)

# average locations of knots (in years) 
knots <- c(1:n.knots)

# make envelopes from bootstrap output
# sometimes there are NaNs in the output, 
# the next lines remove those bootstrap replicates
haz.out <- rh.boot$t[, ind.h]

notnumber <- rep(0, nrow(haz.out))
for (i in 1:nrow(haz.out)) notnumber[i] <- sum(is.nan(haz.out[i, ]))
which(notnumber>0)  # which lines are dud?

haz.out <- rh.boot$t[-which(notnumber>0), ind.h]
haz.out.exp <- rh.boot$t[-which(notnumber>0), ind.h.exp]


haz.env <- envelope( mat= haz.out) 
haz.env.exp <- envelope(mat = haz.out.exp)
# plotting
data <- data.frame(x = italcent$numdays - 105*365.25,
                   d =  I(italcent$rightcens))
ind.x <- 2*n.knots+8 + c(1:N)
ind.h <- 2*n.knots+8 + N + c(1:N)
ind.h.exp <- 2*n.knots+8 + 2*N + c(1:N)
subseq <- seq(1, by = 25, length.out = 219)
ind.x <- ind.x[subseq]
ind.h <- ind.h[subseq]
ind.h.exp <- ind.h.exp[subseq]
if(figures){
  setwd(fig_dir)
  fig <- "Fig8.tex"
  tikz(fig, width = dwidth, height = dheight, standAlone = TRUE)
}
par(mfrow=c(1, 2), mar = c(4, 4, 0.1, 0.1), bty = "l")
plot(rh.boot$t0[ind.x]+105, rh.boot$t0[ind.h], type="n", ylim=c(0, 1.5), 
     ylab="hazard function (1/year)", xlab="years", yaxs = "i",
     panel.first = {
       abline(v=knots + 105, lty=2, col="grey"); 
       rug(105 + jitter(data$x[as.logical(I(data$x >365.25*3)*I(data$d==1))])/365.25, ticksize=-0.02, col="red");
       rug(105 + jitter(data$x[as.logical(I(data$x >365.25*3)*I(data$d==0))])/365.25, ticksize=0.02);
     } )

for (r in 1:100) lines(rh.boot$t[r, ind.x] + 105, 
                       rh.boot$t[r, ind.h], col="grey", lwd=0.5)
lines(rh.boot$t0[ind.x] + 105, rh.boot$t0[ind.h], lwd=2)


lines(rh.boot$t0[ind.x]+105, haz.env$point[1,subseq ], lty=3,lwd = 2)
lines(rh.boot$t0[ind.x]+105, haz.env$point[2, subseq], lty=3,lwd = 2)
lines(rh.boot$t0[ind.x]+105, haz.env$overall[1, subseq], lty=2,lwd = 2)
lines(rh.boot$t0[ind.x]+105, haz.env$overall[2, subseq], lty=2,lwd = 2)

plot(rh.boot$t0[ind.x] + 105, rh.boot$t0[ind.h.exp], type="n", ylim=c(0, 1.5),
     ylab="hazard function (1/year)", xlab="years", yaxs = "i",
     panel.first = {
       abline(v = 105 + knots, lty=2, col="grey"); 
       rug(105 + jitter(data$x[as.logical(I(data$x >365.25*3)*I(data$d==1))])/365.25, ticksize=-0.02, col="red");
       rug(105 + jitter(data$x[as.logical(I(data$x >365.25*3)*I(data$d==0))])/365.25, ticksize=0.02) })

for (r in 1:100) lines(rh.boot$t[r, ind.x]+105, rh.boot$t[r, ind.h.exp], col="grey", lwd=0.5)
lines(rh.boot$t0[ind.x] + 105, rh.boot$t0[ind.h.exp], lwd=2)
lines(rh.boot$t0[ind.x] + 105, haz.env.exp$point[1, subseq], lty=3,lwd=2)
lines(rh.boot$t0[ind.x] + 105, haz.env.exp$point[2, subseq], lty=3,lwd=2)
lines(rh.boot$t0[ind.x] + 105, haz.env.exp$overall[1, subseq], lty=2,lwd=2)
lines(rh.boot$t0[ind.x] + 105, haz.env.exp$overall[2, subseq], lty=2,lwd=2)

if(figures){
  dev.off()
  setwd(code_dir)
}

#########################################################################################
#########################################################################################
# Analysis of "France 2019", i.e., French semi-supercentenarian from IDL,
# Data extracted October 2019
# Table 1
load("francent.rda")
xcal <- francent$birth + u
# Calendar time for sampling frame
c1a <- lubridate::dmy("01-01-1978")
c1b <- lubridate::dmy("01-01-1987")
c2 <- lubridate::dmy("31-12-2017")

# Lower truncation level, zero if individual reached 105 between c1 and c2
francent$slow <- ifelse(francent$numdays > 110*365.24, 
                        as.numeric(pmax(0, c1b - xcal)),
                        as.numeric(pmax(0, c1a - xcal)))
francent$supp <- as.numeric(pmax(0, c2 - xcal))


# Sequence of thresholds, from 105 years (in days), increments of half years
thresh <- seq(from = 0, by = 365.25/2, length = 15L)

# Fit a generalized Pareto model
param_gpd <- matrix(0, nrow = length(thresh), ncol = 10L)
colnames(param_gpd) <- c("loc","scale","shape","deviance", "scale.stderror", "shape.stderror","endpoint", "nu", "pval","pgeqzero")
param_exp <- matrix(0, nrow = length(thresh), ncol = 5L)
colnames(param_exp) <- c("loc","scale","deviance", "scale.stderror", "nu")

dat <- francent$numdays - u
slow <- francent$slow
supp <- francent$supp
for(i in 1:length(thresh)){
  #For each threshold, compute new threshold exceedances
  datu <- dat - thresh[i]
  # Keep only exceedances
  ind <- which(datu > 0)
  vals <- optim(par = fit.gpd(datu/365.25)$est, fn = gpd_dtrunc, 
                method = "N", dat = (datu[ind])/365.25,
                supp = pmax(0,(supp[ind]-thresh[i]))/365.25, 
                slow = (pmax(0, slow[ind]-thresh[i]))/365.25, expo = FALSE, 
                control = list(parscale = c(1,0.01), 
                               reltol = 1e-10, maxit= 1e5),
                hessian = TRUE)
  vals2 <- optim(par = fit.gpd(datu/365.25)$est[1], fn = gpd_dtrunc, 
                 method = "Brent", dat = (datu[ind])/365.25,
                 supp = pmax(0,(supp[ind]-thresh[i]))/365.25, 
                 slow = (pmax(0, slow[ind]-thresh[i]))/365.25, expo = TRUE,
                 lower = 0.05, upper = 200,
                 control = list(reltol = 1e-10, maxit= 1e5),
                 hessian = TRUE)
  param_exp[i,] <- c(loc = (u + thresh[i])/365.25, vals2$par, -2*vals2$value, sqrt(diag(solve(vals2$hessian))),
                     length(ind))
  
  param_gpd[i,] <- c(loc = (u + thresh[i])/365.25, vals$par, -2*vals$value, sqrt(diag(solve(vals$hessian))),
                     ifelse(vals$par[2] >0, Inf, (u + thresh[i])/365.25 - vals$par[1]/vals$par[2]),  
                     length(ind), pchisq(-2*vals$value - param_exp[i,3], 1, lower.tail = FALSE),
                     pnorm(sign(vals$par[2])*sqrt(-2*vals$value - param_exp[i,'deviance'])))
}

  # Table of parameters for French data (Table 1)
  table.res <- matrix(nrow = 6, ncol = 7)
  ind <- seq(1L, 13L, by = 2L)
  table.res[1,] <- paste0("$", round(param_gpd[ind,"nu"],0),"$")
  table.res[2,] <- paste0("$", round(param_gpd[ind, "scale"],2), "\\; (",round(param_gpd[ind, "scale.stderror"],2),")$")
  table.res[3,] <- paste0("$", round(param_gpd[ind, "shape"],2), "\\; (",round(param_gpd[ind, "shape.stderror"],2),")$")
  table.res[4,] <- paste0("$", round(param_exp[ind, "scale"],2), "\\; (",round(param_exp[ind, "scale.stderror"],2),")$")
  table.res[5,] <- paste0("$", round(param_gpd[ind,9],2),"$")
  table.res[6,] <- paste0("$", round(param_gpd[ind,10],2),"$")
  rownames(table.res) <- c("$n_u$", "$\\sigma$", "$\\gamma$", "$\\sigma_e$", "$p$-value","$p_{\\infty}$")
  xtab1 <- xtable::xtable(table.res, 
                          caption = "Estimates (standard errors) of scale  and shape parameters ($\\sigma$, $\\gamma$) for the generalized Pareto distribution and of the scale parameter ($\\sigma_e$) for the exponential model for the French \\IDL{} data as a function of threshold, with number of threshold exceedances ($n_u$), $p$-value for the likelihood ratio test of $\\gamma=0$ and probability under the fitted model that the endpoint is infinite, corresponding to $\\Pr(\\gamma>0)$.",
                          align = paste0(rep("r",ncol(table.res)+1),collapse = ""))
  addtorow <- list()
  addtorow$pos <- list(0)
  addtorow$command <-  paste0("threshold", paste0(" & $", 105:111, "$", collapse = " "), "\\\\\n", collapse = " ")
  
  if(tables){
  setwd(table_dir)
  print(xtab1, file = "Table1b.tex", 
        size = "footnotesize", 
        floating.environment = "table*",   
        booktabs = TRUE, caption.placement = "top",
        add.to.row = addtorow, include.colnames = FALSE,
        sanitize.colnames.function = identity, 
        sanitize.text.function = identity, 
        sanitize.rownames.function = identity)
  setwd(code_dir)
}


mint <- francent$gender == "male"
wint <- francent$gender == "female"

fjit <- fitgpddtrunc(dat = dat, thresh = thresh[1:11], supp = supp, slow = slow)
fmen <- fitgpddtrunc(dat = dat[mint], thresh = thresh[1:11], supp = supp[mint], slow = slow[mint])
fwom <- fitgpddtrunc(dat = dat[wint], thresh = thresh[1:11], supp = supp[wint], slow = slow[wint])

oyr <- seq(1,11, by =2)
pvals <- pchisq((fwom[oyr,4] + fmen[oyr,4])-fjit[oyr,4], df = 2, lower.tail = FALSE)
pvalt <- rbind(paste0("$",round(100*pvals, digits= 2), "$"))
colnames(pvalt) <- paste0( "$", as.vector(round(fjit[oyr,1],1)), "$")

((fmen[7,1] + fmen[7,2]/(1-fmen[7,3])) - (fwom[7,1] + fwom[7,2]/(1-fwom[7,3])) )*365.25


ejit <- fitgpddtrunc(dat = dat, thresh = thresh[1:11], supp = supp, slow = slow, expo = TRUE)
emen <- fitgpddtrunc(dat = dat[mint], thresh = thresh[1:11], supp = supp[mint], slow = slow[mint], expo = TRUE)
ewom <- fitgpddtrunc(dat = dat[wint], thresh = thresh[1:11], supp = supp[wint], slow = slow[wint], expo = TRUE)
pvalt <- rbind(pvalt, paste0("$",round(100*pchisq((ewom[oyr,3] + emen[oyr,3])-ejit[oyr,3], df = 1, lower.tail = FALSE), 2),"$"))
# All differences are significative, except 110
rownames(pvalt) <- c("generalized Pareto", "exponential")


#########################################################################################
#########################################################################################
# Fig 5 Threshold stability plots for France 2019 data

dat <- francent$numdays - u
slow <- francent$slow
supp <- francent$supp


confint_xi <- matrix(0, ncol = 3, nrow = length(thresh))
confint_exp <- matrix(0, ncol = 3, nrow = length(thresh))
for(i in 1:length(thresh)){
  datu <- dat - thresh[i]
  ind <- which(datu > 0)
  confint_xi[i,] <- prof_gpd_dtrunc_xi_confint(dat = (datu[ind])/365.25,
                                               supp = (supp[ind]-thresh[i])/365.25, slow = (pmax(0, slow[ind]-thresh[i]))/365.25)
  confint_exp[i,] <- prof_exp_dtrunc(dat = (datu[ind])/365.25,
                                     supp = (supp[ind]-thresh[i])/365.25, slow = (pmax(0, slow[ind]-thresh[i]))/365.25)
  
}

if(figures){
  fig <- "Fig5.tex"
  setwd(fig_dir)
  tikz(fig, width = dwidth, height = 0.8*dheight, standAlone = TRUE)
}
# Consider two thresholds for which xi is negative
par(mar = c(4, 4, 0.1, 0.1), mfrow = c(1,2))
plot(x = thresh/365.25+105, confint_xi[,1], 
     type= "p", 
     panel.first = abline(h = 0, col = "grey"), 
     pch = 20, col = 1, bty = "l", 
     ylab = "$\\gamma$", 
     xlab = "threshold (in years)",
     ylim = c(-0.3,1.1))
for(i in 1:length(thresh)){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_xi[i,2], y1 = confint_xi[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
plot(x = thresh/365.25+105, confint_exp[,1], type= "p", 
     panel.first = abline(h = 1.42, col = "grey"),
     pch = 20, col = 1,bty = "l", ylab = "$\\sigma_e$", xlab = "threshold (in years)",ylim = range(confint_exp))
for(i in 1:length(thresh)){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_exp[i,2], y1 = confint_exp[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
if(figures){
  dev.off()
  system(command = paste0("pdflatex ",fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}

