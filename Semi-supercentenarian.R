rm(list=ls())
# Three files are used, italcent.rda, idl2016.rda and francent.rda
# The first can be bought for a small fee from the
# National Institute of Statistics by registering at the Contact Center
# (https://contact.istat.it) and mentioning the Semisupercentenarian
# Survey and Marco Marsili as contact person.
# 
# The other two datsets can be downloaded by registering on 
# http://www.supercentenarians.org/, and/or by forwarding proof 
# of registration to Leo Belzile
# Load libraries and packages used in the analyses
# 
# Figure numbers refer to arXiv version of the paper
library(boot)  
library(evd)
library(mev)
library(xts)
library(lubridate)
library(survival)
library(xtable)
library(tikzDevice)
options(tikzLatexPackages =
          c("\\usepackage{tikz}\n",
            "\\usepackage[active,tightpage,psfixbb]{preview}\n",
            "\\usepackage{amsmath}",
            "\\PreviewEnvironment{pgfpicture}\n",
            "\\setlength\\PreviewBorder{0pt}\n",
            "\\usepackage{fourier}\n",
            "\\DeclareMathAlphabet{\\mathdis}{OT1}{pag}{m}{n}\n"
          )
)
# The following package is not available on the CRAN
# devtools::install_github("OpenIntroStat/openintro-r-package", subdir = "OIsurv")
library(OIsurv)
source("Semi-supercentenarian_fn.R")
dwidth <- 4
dheight <- 2.5
set.seed(1234)
#Set directory when sourcing file
# setwd(utils::getSrcDirectory()[1])
#Set directory from Rstudio to current document
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# TODO create figure and tables directories
code_dir <- getwd()
fig_dir <- paste0(getwd(), "/figure")
if(!dir.exists(fig_dir)){
  dir.create("figure")
}
table_dir <- paste0(getwd(), "/tables")
if(!dir.exists(table_dir)){
  dir.create("tables")
}

# To save Figures and Tables (.tex format), change the following to TRUE
figures <- tables <- TRUE
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
     bty = "l", pch = 20, cex = 0.8, ylab = "age at death (in years)", xaxt = 'n',
     xlab = "year of death", yaxs = "i", ylim = c(105, 117), xlim = c(14244, 16825), xaxs = "i")
axis(side=1,at=seq(14244, 16825, by = 365.25), labels = 2009:2016, tick = TRUE)
points(x = (italcent$birth + italcent$numdays)[menunc], 
       y = italcent$numdays[menunc]/365.25,  pch = 4, cex = 0.8,
       col = scales::alpha("red", alpha = 0.5),lwd = 2)
rug(italcent$numdays[mencens]/365.25, side = 4, line =  1, col = scales::alpha("red", 0.5), ticksize = 0.05)
rug(italcent$numdays[womcens]/365.25, side = 4, line =  2, col = scales::alpha("black", 0.5), ticksize = 0.05)
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
     pch = 20, bty = "l", xlab = "theoretical quantiles", ylab = "observed quantiles",
     panel.first = {abline(a=0, b=1)}, yaxs = "i", xaxs = "i",
     ylim = 108+c(0, 3000)/365.25, xlim = 108+c(0, 3000)/365.25)
for(i in 1:100){
  lines(y = sort(bootsamp[i,!cens108])/365.25+108, qqptsb[j,]/365.25+108, 
        lty = 2, col = scales::alpha(colour = "black", 0.2))
}
points(y = sort(numd108[!cens108])/365.25+108, sort(qexp(cecdf(numd108[!cens108]), rate = 1/param_exp108))/365.25+108, pch = 1, cex = 0.8)
points(y = sort(numd108[!cens108])/365.25+108, sort(qexp(cecdf(numd108[!cens108]), rate = 1/param_exp108))/365.25+108, pch = 20, cex = 0.2, col = "white")
plot(y = sort(numd108[!cens108])/365.25+108, sort(qexp(cecdf(numd108[!cens108]), rate = 1/param_exp108))/365.25+108,
     pch = 20, bty = "l", xlab = "theoretical quantiles", ylab = "observed quantiles",
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
  system(command = paste0("pdflatex ",fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}


#########################################################################################
#########################################################################################
# Power study - comparison of men versus women (SIM, Section A.5)

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
     panel.first = {abline(h = seq(0.1, 1, by = 0.1), col = scales::alpha("black", 0.05), lty = 1, lwd = 0.25)})
gam1 <- mgcv::gam(colMeans(lrt_stata > qchisq(0.95,1)) ~ s(lam_vals^2, k = -1, bs = "tp"))
lines(lam_vals^2, fitted(gam1), lwd = 2)
gam2 <- mgcv::gam(colMeans(lrt_statb > qchisq(0.95,1)) ~ s(lam_vals^2, k = -1, bs = "tp"))
lines(lam_vals^2, fitted(gam2), lwd = 2, lty = 3,  col = "grey30")
legend("bottomright", legend = c("censoring","no censoring"), bty = "n",
       col = c(1, "grey30"), lwd = 2, lty = c(1,3))
setwd(code_dir)


#########################################################################################
#########################################################################################
# Figure 2a: parameter stability plots for Istat data

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
  fig <- "Fig2a.tex"
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
     panel.first = abline(h=1.45142127, col = "grey"), pch = 20, col = 1,
     bty = "l", ylab = "$\\sigma_e$", xlab = "threshold (in years)",
     ylim = c(0.7,2.1))
for(i in 1:length(thresh)){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_exp[i,2], y1 = confint_exp[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}

if(figures){
  dev.off()
  system(command = paste0("pdflatex ", fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}




#########################################################################################
#########################################################################################
# Fig. 5 Parameter stability plot, by cohort 1896-1905 vs 1906-1910

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
  fig <- "Fig5.tex"
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
     bty = "l", ylab = "$\\sigma_e$", xlab = "threshold (in years)",
     ylim = c(0.7,2.1))
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
     bty = "l", ylab = "$\\sigma_e$", xlab = "threshold (in years)",
     ylim = c(0.7, 2.1))
for(i in 1:8){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_exp2[i,2], y1 = confint_exp2[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
dev.off()
system(command = paste0("pdflatex ", fig_dir, "/", fig,"; rm *.aux; rm *.log"))
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
yran <- range(c(0, c(localhazard1[,1:3], localhazard2[,1:3])))
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
  system(command = paste0("pdflatex ", fig_dir, "/", fig,"; rm *.aux; rm *.log"))
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
# Fig. 7 Local GP fit using splines with bootstrap analysis of variability
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
system.time( rh.boot <- boot(data=data, statistic=boot.stat1, R=5000L, 
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
  fig <- "Fig7.tex"
  tikz(fig, width = dwidth, height = dheight, standAlone = TRUE)
}
par(mfrow=c(1, 2), mar = c(4, 4, 0.1, 0.4), bty = "l")
plot(rh.boot$t0[ind.x]+105, rh.boot$t0[ind.h], type="n", ylim = c(0, 1.5), 
     ylab="hazard function (1/year)", xlab="age (in years)", yaxs = "i",
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
     ylab="hazard function (1/year)", xlab="age (in years)", yaxs = "i",
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
  system(command = paste0("pdflatex ",fig_dir, "/", fig,"; rm *.aux; rm *.log"))
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
# Fig 2 (bottom) Threshold stability plots for France 2019 data

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
  fig <- "Fig2b.tex"
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
     panel.first = abline(h = 1.45142127, col = "grey"),
     pch = 20, col = 1,bty = "l", ylab = "$\\sigma_e$", xlab = "threshold (in years)",
     ylim = c(0.7,2.1))
for(i in 1:length(thresh)){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_exp[i,2], y1 = confint_exp[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
if(figures){
  dev.off()
  system(command = paste0("pdflatex ",fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}

#########################################################################################
#########################################################################################
# Power plots (Figure 3)
# Power of GP for H0: xi>=0 versus Ha:xi < 0 (finite endpoint)
# Three database: IDL, France (2019) and IStat
# We condition for all three on entry points
# Right-censor for Istat data, else simulate doubly truncated data

# Compute both directed likelihood root and Wald statistics and compare power

# Power analysis for one-sided test with Ha: gamma>=0 versus H0: gamma < 0
# and for two-sided tests with Ha gamma=0 versus H0: gamma != 0 (not shown in text)
# 
# 

set.seed(20200618)
B <- 1e4
#Set directory from Rstudio to current document
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
#setwd("/home/lbelzile/Documents/Dropbox/SuperCentenarian/code")
source("Semi-supercentenarian_fn.R")
figures <- FALSE
fig_dir <- paste0(getwd(), "/figure")
library(progress)

u108 <- 39447L #108 years
xis <- seq(-0.25, 0.25, length.out = 101L)
nxis <- length(xis)
save <- FALSE # Uncomment to run the power study

load("italcent.rda")
# Calendar date at which individual reaches 108 years
italcent$xcal <- italcent$birth + u108
# Calendar time for sampling frame
c1 <- lubridate::dmy("01-01-2009")
c2 <- lubridate::dmy("01-01-2016")

# Lower truncation level, zero if individual reached 108 between c1 and c2
italcent$slow <- as.numeric(pmax(0, c1 - italcent$xcal))
# Exceedances
datu <- italcent$numdays - u108
# Fit a generalized Pareto model
ind <- which(datu > 0)
uplim <- as.numeric(c2 - italcent$birth)
# Keep only exceedances above 108
uplim <- (uplim[ind] - u108)/365.25
rcens <- italcent$rightcens[ind]
datu <- datu[ind]/365.25
slow <- italcent$slow[ind]/365.25

itex <- data.frame(datu = datu, slow = slow, rcens = rcens, uplim = uplim)
profile_italcent <- t(sapply(1:nxis, function(i){
  prof_gpd_cens_xi(xi = xis[i], 
                   dat = itex$datu, 
                   rightcens = itex$rcens, 
                   slow = itex$slow)}))[,3:2]

load("francent.rda")
# Calendar date at which individual reaches 105 years
xcal <- francent$birth + u108
# Calendar time for sampling frame
c1a <- lubridate::dmy("01-01-1978")
c1b <- lubridate::dmy("01-01-1987")
c2 <- lubridate::dmy("31-12-2017")

# Lower truncation level, zero if individual reached 105 between c1 and c2
francent$slow <- ifelse(francent$numdays > 110*365.24, 
                        as.numeric(pmax(0, c1b - xcal)),
                        as.numeric(pmax(0, c1a - xcal)))
francent$supp <- as.numeric(pmax(0, c2 - xcal))    
datu <- francent$numdays - u108
ind <- which(datu > 0)
slow <- francent$slow[ind]/365.25
supp <- francent$supp[ind]/365.25
datu <- (francent$numdays - u108)[ind]/365.25

frex <- data.frame(datu = datu, slow = slow, supp = supp)
profile_francent <- prof_gpd_dtrunc_xi(xi = xis,
                                       dat = frex$datu, 
                                       supp = frex$supp, 
                                       slow = frex$slow)$param



load("IDL2016.rda")
idl2016 <- idl2016[!idl2016$countrydeath == "FRA",]
u110 <- min(idl2016$numdays)-1 #some people are 2-3 days from 110...
# Calendar date at which individual reaches 105 years
datu <- idl2016$numdays - u110
slow <- pmax(0, (idl2016$slow - u110)/365.25)
supp <- (idl2016$supp - u110)/365.25
datu <- datu/365.25

idlex <- data.frame(datu = datu, slow = slow, supp = supp)
# Why is lower truncation level sometimes higher than the excess lifetime (in days)?

bootsamp <- matrix(0, ncol = length(datu), nrow = B)
prof_xi_idl <- prof_gpd_dtrunc_xi(xi = xis,
                                  dat = idlex$datu, 
                                  supp = idlex$supp, 
                                  slow = idlex$slow)$param
profile_idl <- prof_xi_idl
if(save){
  set.seed(20200824)
  power_italcent <- array(0, dim = c(B, nxis, 2))
  power_idl <- array(0, dim = c(B, nxis, 2))
  power_francent <- array(0, dim = c(B, nxis, 2))
  shapes <- array(0, dim = c(B, 3, nxis, 2))
  bootsampital <- matrix(0, ncol = length(itex$datu), nrow = B)
  bootsampfr <- matrix(0, ncol = length(frex$datu), nrow = B)
  bootsampidl <- matrix(0, ncol = length(idlex$datu), nrow = B)
  rightcensb <- matrix(FALSE, ncol = length(itex$datu), nrow = B)
  
  
  for(i in 1:nxis){
    #Generate new data
    for(j in 1:ncol(bootsampfr)){
      bootsampfr[,j] <- revddtrunc(n = B, param = profile_francent[i,], 
                                   lower = frex$slow[j], 
                                   upper = frex$supp[j])
    }
    for(j in 1:ncol(bootsampidl)){
      bootsampidl[,j] <- revddtrunc(n = B, param = profile_idl[i,], 
                                    lower = frex$slow[j], 
                                    upper = frex$supp[j])
    }
    for(j in 1:ncol(bootsampital)){
      bootsampital[,j] <- evd::qgpd(evd::pgpd(itex$slow[j], loc = 0, scale = profile_italcent[i,1], shape = profile_italcent[i,2]) + runif(B)*
                                      (1-evd::pgpd(itex$slow[j], loc = 0, scale = profile_italcent[i,1], shape = profile_italcent[i,2])), 
                                    loc = 0, scale = profile_italcent[i,1], shape = profile_italcent[i,2])
      rightcensb[,j] <- bootsampital[,j] > itex$uplim[j]
      bootsampital[rightcensb[,j],j] <- itex$uplim[j]
    }
    for(b in 1:B){
      #Find maximum likelihood estimates
     
      mleboot_it <- try(optim(par =  profile_italcent[i,],
                              fn = gpd_cens,
                              method = "N",
                              dat = bootsampital[b,],
                              rightcens = rightcensb[b,],
                              slow = itex$slow,
                              expo = FALSE,
                              control = list(fnscale = c(1.5,0.01), reltol = 1e-12, maxit = 1e5),
                              hessian = TRUE))
      mleboot_fr <- try(optim(par =  profile_francent[i,],
                              fn = gpd_dtrunc,
                              method = "N",
                              dat = bootsampfr[b,],
                              supp = frex$supp,
                              slow = frex$slow,
                              expo = FALSE,
                              control = list(fnscale = c(1.5,0.01), reltol = 1e-12, maxit= 1e5),
                              hessian = TRUE))
      mleboot_idl <- try(optim(par =  profile_idl[i,],
                               fn = gpd_dtrunc,
                               method = "N",
                               dat = bootsampidl[b,],
                               slow = idlex$slow,
                               supp = idlex$supp,
                               expo = FALSE,
                               control = list(fnscale = c(1.5,0.01), reltol = 1e-12, maxit= 1e5),
                               hessian = TRUE))
      
      shapes[b,1,i,] <- c(mleboot_it$par[2], sqrt(diag(solve(mleboot_it$hessian))[2]))
      shapes[b,2,i,] <- c(mleboot_fr$par[2], sqrt(diag(solve(mleboot_fr$hessian))[2]))
      shapes[b,3,i,] <- c(mleboot_idl$par[2], sqrt(diag(solve(mleboot_idl$hessian))[2]))
    
      mbootfr <- mean(bootsampfr[b,])
      mbootidl <- mean(bootsampidl[b,])
      pexpboot_idl <- try(optim(par = mbootidl, 
                                fn = gpd_dtrunc, 
                                method = "Brent", 
                                lower = 0.01*mbootidl, 
                                upper = 3*mbootidl, 
                                dat = bootsampidl[b,],
                                supp = idlex$supp, 
                                slow = idlex$slow, 
                                expo = TRUE))
      pexpboot_fr <- try(optim(par = mbootfr, 
                               fn = gpd_dtrunc, 
                               method = "Brent", 
                               lower = 0.01*mbootfr, 
                               upper = 3*mbootfr, 
                               dat = bootsampfr[b,],
                               supp = frex$supp, 
                               slow = frex$slow, 
                               expo = TRUE))
      rate_it <- exp_mle_lt_rc(dat =  bootsampital[b,], 
                               rightcens = rightcensb[b,], 
                               slow = itex$slow)
      pexpboot_it <- list(par = rate_it, 
                          value = gpd_cens(par = c(rate_it, 0), 
                                           dat =  bootsampital[b,], 
                                           rightcens = rightcensb[b,], 
                                           slow = itex$slow))
      
      if(is.character(mleboot_it)){
        power_italcent[b,i,] <- rep(NA, 2)
      } else{
        #Compute directed  profile likelihood ratio root test statistic
        power_italcent[b, i, 1] <- sign(mleboot_it$par[2])*
          sqrt(2*(pexpboot_it$value - mleboot_it$value))
        #Compute Wald-test
        if(mleboot_it$par[2] > -0.5){
          power_italcent[b,i,2] <- mleboot_it$par[2]/sqrt(diag(solve(mleboot_it$hessian))[2])
        } else{ 
          power_italcent[b,i,2] <- NA
        }
      }
      
      if(is.character(mleboot_fr) || is.character(pexpboot_fr)){
        power_francent[b,i,] <- rep(NA, 2)
      } else{
        #Compute directed profile likelihood ratio root test statistic
        power_francent[b, i, 1] <- sign(mleboot_fr$par[2])*
          sqrt(2*(pexpboot_fr$value - mleboot_fr$value))
        #Compute Wald-test
        if(mleboot_fr$par[2] > -0.5){
          power_francent[b,i,2] <- mleboot_fr$par[2]/sqrt(diag(solve(mleboot_fr$hessian))[2])
        } else{ 
          power_francent[b,i,2] <- NA
        }
      }
      
      if(is.character(mleboot_idl) || is.character(pexpboot_idl)){
        power_idl[b,i,] <- rep(NA, 2)
      } else{
        #Compute directed profile likelihood ratio root test statistic
        power_idl[b, i, 1] <- sign(mleboot_idl$par[2])*
          sqrt(-2*(mleboot_idl$value - pexpboot_idl$value))
        #Compute Wald-test
        if(mleboot_idl$par[2] > -0.5){
          power_idl[b,i,2] <- mleboot_idl$par[2]/sqrt(diag(solve(mleboot_idl$hessian))[2])
        } else{ 
          power_idl[b,i,2] <- NA
        }
      }
    }
  }
  save(xis, power_italcent, power_idl, power_francent, shapes, file = "power.RData")
} else{
  load("power.RData") 
}

# Power for the endpoint

set.seed(20200618)
B <- 1e4
#Set directory from Rstudio to current document
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
#setwd("/home/lbelzile/Documents/Dropbox/SuperCentenarian/code")
source("Semi-supercentenarian_fn.R")
figures <- TRUE
fig_dir <- paste0(getwd(), "/figure")
library(progress)

u108 <- 39447L #108 years
endpoints <- c(seq(116.2, 125, by = 0.2), seq(126, 150, by = 1))
nendpoints <- length(endpoints)
save <- FALSE

load("italcent.rda")
# Calendar date at which individual reaches 108 years
italcent$xcal <- italcent$birth + u108
# Calendar time for sampling frame
c1 <- lubridate::dmy("01-01-2009")
c2 <- lubridate::dmy("01-01-2016")

# Lower truncation level, zero if individual reached 108 between c1 and c2
italcent$slow <- as.numeric(pmax(0, c1 - italcent$xcal))
maxital <- max(italcent$numdays/365.25)
# Exceedances
datu <- italcent$numdays - u108
# Fit a generalized Pareto model
ind <- which(datu > 0)
uplim <- as.numeric(c2 - italcent$birth)
# Keep only exceedances above 108
uplim <- (uplim[ind] - u108)/365.25
rcens <- italcent$rightcens[ind]
datu <- datu[ind]/365.25
slow <- italcent$slow[ind]/365.25

itex <- data.frame(datu = datu, slow = slow, rcens = rcens, uplim = uplim)
profile_italcent <- t(sapply(1:nendpoints, function(i){
  prof_gpd_cens_endpoint(endpoint = (endpoints[i] - 108), 
                         dat = itex$datu, 
                         rightcens = itex$rcens, 
                         slow = itex$slow)}))[,2:3]

load("francent.rda")
# Calendar date at which individual reaches 105 years
xcal <- francent$birth + u108
# Calendar time for sampling frame
c1a <- lubridate::dmy("01-01-1978")
c1b <- lubridate::dmy("01-01-1987")
c2 <- lubridate::dmy("31-12-2017")

# Lower truncation level, zero if individual reached 105 between c1 and c2
francent$slow <- ifelse(francent$numdays > 110*365.24, 
                        as.numeric(pmax(0, c1b - xcal)),
                        as.numeric(pmax(0, c1a - xcal)))
francent$supp <- as.numeric(pmax(0, c2 - xcal))    
datu <- francent$numdays - u108
ind <- which(datu > 0)
slow <- francent$slow[ind]/365.25
supp <- francent$supp[ind]/365.25
datu <- (francent$numdays - u108)[ind]/365.25
maxfran <- max(francent$numdays/365.25)
frex <- data.frame(datu = datu, slow = slow, supp = supp)
profile_francent <- t(sapply(1:nendpoints, function(i){
  prof_gpd_dtrunc_endpoint(endpoint = (endpoints[i] - 108), 
                           dat = frex$datu, slow = frex$slow, supp = frex$supp)}))[,2:3]


load("IDL2016.rda")
idl2016 <- idl2016[!idl2016$countrydeath == "FRA",]
u110 <- min(idl2016$numdays)-1 #some people are 2-3 days from 110...
# Calendar date at which individual reaches 105 years
datu <- idl2016$numdays - u110
slow <- pmax(0, (idl2016$slow - u110)/365.25)
supp <- (idl2016$supp - u110)/365.25
datu <- datu/365.25
maxidl <- max(idl2016$numdays/365.25)
idlex <- data.frame(datu = datu, slow = slow, supp = supp)
# Why is lower truncation level sometimes higher than the excess lifetime (in days)?
profile_idl <- t(sapply(1:nendpoints, function(i){
  prof_gpd_dtrunc_endpoint(endpoint = (endpoints[i] - 110), 
                           dat = idlex$datu, slow = idlex$slow, supp = idlex$supp)}))[,2:3]


if(save){
  set.seed(20200824)
  power_ep_combo <- array(NA, dim = c(B, nendpoints))
  power_ep_francent <- array(NA, dim = c(B, nendpoints))
  power_ep_idl <- array(NA, dim = c(B, nendpoints))
  power_ep_italcent <- array(NA, dim = c(B, nendpoints))
  
  
  bootsampital <- matrix(0, ncol = length(itex$datu), nrow = B)
  bootsampfr <- matrix(0, ncol = length(frex$datu), nrow = B)
  bootsampidl <- matrix(0, ncol = length(idlex$datu), nrow = B)
  rightcensb <- matrix(FALSE, ncol = length(itex$datu), nrow = B)
  
  combllep <- function(par,  maxfr, maxidl,  maxital, datidl, datfr, datital, rcensit){
    endpoint <- par[1]
    shape <- par[2:4]
    scale <- -shape*(endpoint - c(108,108,110))
    par_it <- c(scale[1], shape[1])
    par_fr <- c(scale[2], shape[2])
    par_idl <- c(scale[3], shape[3])
    ub <- endpoint - c(108,108,110)
    gpd_dtrunc(par_fr, dat = datfr, slow = frex$slow, supp = frex$supp, expo = FALSE) +
      gpd_dtrunc(par_idl, dat = datidl, slow = idlex$slow, supp = idlex$supp, expo = FALSE) +  
      gpd_cens(par_it, dat = datital, rightcens = rcensit, slow = itex$slow, expo = FALSE)
  }
  combllep_scale <- function(par,  maxfr, maxidl,  maxital, datidl, datfr, datital, rcensit){
    endpoint <- par[1]
    scale <- par[2:4]
    shape <- -scale/(endpoint - c(108,108,110))
    par_it <- c(scale[1], shape[1])
    par_fr <- c(scale[2], shape[2])
    par_idl <- c(scale[3], shape[3])
    ub <- endpoint - c(108,108,110)
    gpd_dtrunc(par_fr, dat = datfr, slow = frex$slow, supp = frex$supp, expo = FALSE) +
      gpd_dtrunc(par_idl, dat = datidl, slow = idlex$slow, supp = idlex$supp, expo = FALSE) +  
      gpd_cens(par_it, dat = datital, rightcens = rcensit, slow = itex$slow, expo = FALSE)
  }
  ineq <- function(par, ...){
    args <- list(...)
    c(par[1]-pmax(args$maxital, args$maxfr, args$maxidl), par[2:4] + 1, -par[2:4])}
  ineq_scale <- function(par, ...){
    args <- list(...)
    endpoint <- par[1]
    scale <- par[2:4]
    shape <- -scale/(endpoint - c(108,108,110))
    c(endpoint-pmax(args$maxital, args$maxfr, args$maxidl), shape + 1, -shape)}
  
  for(i in 1:nendpoints){
    #Generate new data
    if(endpoints[i] > maxfran){
      for(j in 1:ncol(bootsampfr)){
        bootsampfr[,j] <- revddtrunc(n = B, param = profile_francent[i,], 
                                     lower = frex$slow[j], 
                                     upper = frex$supp[j])
      }
    }
    if(endpoints[i] > maxidl){
      for(j in 1:ncol(bootsampidl)){
        bootsampidl[,j] <- revddtrunc(n = B, param = profile_idl[i,], 
                                      lower = frex$slow[j], 
                                      upper = frex$supp[j])
      }
    }
    if(endpoints[i] > maxital){
      for(j in 1:ncol(bootsampital)){
        bootsampital[,j] <- evd::qgpd(evd::pgpd(itex$slow[j], loc = 0, scale = profile_italcent[i,1], shape = profile_italcent[i,2]) + runif(B)*
                                        (1-evd::pgpd(itex$slow[j], loc = 0, scale = profile_italcent[i,1], shape = profile_italcent[i,2])), 
                                      loc = 0, scale = profile_italcent[i,1], shape = profile_italcent[i,2])
        rightcensb[,j] <- bootsampital[,j] > itex$uplim[j]
        bootsampital[rightcensb[,j],j] <- itex$uplim[j]
      }
    }
    for(b in 1:B){
       #Find maximum likelihood estimates
      if(endpoints[i] > maxfran){
        maxital = max(bootsampital[b,])
        maxfr = max(bootsampfr[b,])
        maxidl = max(bootsampidl[b,])
          mleboot_combo <- try(alabama::auglag(par = c(pmax(maxital, maxfr, maxidl, endpoints[i]) + 0.5, profile_italcent[i,1], profile_francent[i,1], profile_idl[i,1]),
                                             fn = combllep_scale,
                                             hin = ineq_scale,
                                             maxital = maxital,
                                             maxfr = maxfr,
                                             maxidl = maxidl,
                                             datital = bootsampital[b,],
                                             datfr = bootsampfr[b,],
                                             datidl = bootsampidl[b,],
                                             rcensit = rightcensb[b,],
                                             control.outer=list(method="BFGS", trace=0),
                                             control.optim = list(parscale = c(endpoints[i], profile_italcent[i,1], profile_francent[i,1], profile_idl[i,1]))))
      }
      if(endpoints[i] > maxital){
        mleboot_it <- try(optim(par =  profile_italcent[i,],
                                fn = gpd_cens,
                                method = "N",
                                dat = bootsampital[b,],
                                rightcens = rightcensb[b,],
                                slow = itex$slow,
                                expo = FALSE,
                                control = list(fnscale = c(1.5,0.01), reltol = 1e-12, maxit = 1e5),
                                hessian = TRUE))
        rate_it <- exp_mle_lt_rc(dat =  bootsampital[b,], 
                                 rightcens = rightcensb[b,], 
                                 slow = itex$slow)
        pexpboot_it <- list(par = rate_it, 
                            value = gpd_cens(par = c(rate_it, 0), 
                                             dat =  bootsampital[b,], 
                                             rightcens = rightcensb[b,], 
                                             slow = itex$slow))
        if(!is.character(mleboot_it)){
          #Compute directed profile likelihood ratio root test statistic
          power_ep_italcent[b, i] <- sign(mleboot_it$par[2])*sqrt(2*(pexpboot_it$value - mleboot_it$value))
        }
        
      }
      if(endpoints[i] > maxfran){
        mleboot_fr <- try(optim(par =  profile_francent[i,],
                                fn = gpd_dtrunc,
                                method = "N",
                                dat = bootsampfr[b,],
                                supp = frex$supp,
                                slow = frex$slow,
                                expo = FALSE,
                                control = list(fnscale = c(1.5,0.01), reltol = 1e-12, maxit= 1e5),
                                hessian = TRUE))
        mbootfr <- mean(bootsampfr[b,])
        pexpboot_fr <- try(optim(par = mbootfr, 
                                 fn = gpd_dtrunc, 
                                 method = "Brent", 
                                 lower = 0.01*mbootfr, 
                                 upper = 3*mbootfr, 
                                 dat = bootsampfr[b,],
                                 supp = frex$supp, 
                                 slow = frex$slow, 
                                 expo = TRUE))
        if(!(is.character(mleboot_fr) || is.character(pexpboot_fr))){
          #Compute directed profile likelihood ratio root test statistic
          power_ep_francent[b, i] <- sign(mleboot_fr$par[2])*sqrt(2*(pexpboot_fr$value - mleboot_fr$value))
        }
      }
      if(endpoints[i] > maxidl){
        mleboot_idl <- try(optim(par =  profile_idl[i,],
                                 fn = gpd_dtrunc,
                                 method = "N",
                                 dat = bootsampidl[b,],
                                 slow = idlex$slow,
                                 supp = idlex$supp,
                                 expo = FALSE,
                                 control = list(fnscale = c(1.5,0.01), reltol = 1e-12, maxit= 1e5),
                                 hessian = TRUE))
        mbootidl <- mean(bootsampidl[b,])
        pexpboot_idl <- try(optim(par = mbootidl, 
                                  fn = gpd_dtrunc, 
                                  method = "Brent", 
                                  lower = 0.01*mbootidl, 
                                  upper = 3*mbootidl, 
                                  dat = bootsampidl[b,],
                                  supp = idlex$supp, 
                                  slow = idlex$slow, 
                                  expo = TRUE))
        if(!(is.character(mleboot_idl) || is.character(pexpboot_idl))){
          power_ep_idl[b, i] <- sign(mleboot_idl$par[2])*sqrt(-2*(mleboot_idl$value - pexpboot_idl$value))
        }
      }
      if(endpoints[i] > maxfran){
        if(!any(is.character(mleboot_combo), 
                is.character(pexpboot_idl), 
                is.character(pexpboot_fr))){
          power_ep_combo[b, i] <- -sqrt(pmax(0,-2*(mleboot_combo$value - (pexpboot_fr$value + pexpboot_it$value + pexpboot_idl$value))))
        }
      } 
    }
    }
  save(endpoints, power_ep_italcent, power_ep_idl, power_ep_combo, power_ep_francent, file = "power_endpoint.RData")
} else{
  load("power_endpoint.RData") 
}
powerIDL_ep <- colMeans(power_ep_idl < qnorm(0.05), na.rm = TRUE)
powerIstat_ep <- colMeans(power_ep_italcent < qnorm(0.05), na.rm = TRUE)
powerFrance_ep <- colMeans(power_ep_francent < qnorm(0.05), na.rm = TRUE)
powerCombo_ep <- colMeans(power_ep_combo < qnorm(0.05), na.rm = TRUE)


# Smoothing of the power curve to get rid of the artifacts of Monte-Carlo variability
library(cobs)
library(ggplot2)
library(poorman)
library(viridis)
library(patchwork)
library(tikzDevice)
nknots <- 40L
endpts <- seq(115, 180, by = 0.1)
p1IT <- is.na(powerIDL_ep)
powersmoothIstat <- predict(cobs::cobs(
  x = endpoints[-p1IT], 
  y = powerIstat_ep[-p1IT],
  pointwise = cbind(0, maxital, 1),   
  constraint = "convex", 
  nknots = nknots), z = endpts)[,"fit"]
powersmoothIstat[endpts < maxital] <- 1

p1FR <- is.na(powerFrance_ep)
powersmoothFrance <- predict(
  cobs::cobs(x = endpoints[-p1FR], 
             y = powerFrance_ep[-p1FR], 
             pointwise = cbind(0, maxfran, 1),
             constraint = "convex",
             nknots = nknots),
  z = endpts)[,"fit"]
powersmoothFrance[endpts < maxfran] <- 1
plot(endpts, powersmoothFrance, type = "l")
points(endpoints, powerFrance_ep)

p1IDL <- is.na(powerIDL_ep)
powersmoothIDL <- predict(cobs::cobs(
  x = endpoints[-p1IDL], 
  y = powerIDL_ep[-p1IDL], 
  pointwise = cbind(0, maxidl, 1),
  constraint = "convex", nknots = nknots),
  z = endpts)[,"fit"]
powersmoothIDL[endpts < maxidl] <- 1


powersmoothCombo <- predict(
  cobs::cobs(x = endpoints[-p1FR], 
             y = powerCombo_ep[-p1FR], 
             pointwise = cbind(0, maxfran, 1),
             constraint = "convex", 
             nknots = nknots),
  z = endpts)[,"fit"]
powersmoothCombo[endpts < maxfran] <- 1



power_endp <- data.frame(
  endpoint = rep(endpts, length.out = 4*length(endpts)),
  data = factor(rep(c("Istat","France","IDL2016","combined"), each = length(endpts))),
  power = c(powersmoothIstat, powersmoothFrance, powersmoothIDL, powersmoothCombo)
)

power_df_ep <- data.frame(
  endpoint = rep(endpoints, length.out = 4*length(endpoints)),
  data = factor(rep(c("Istat","France","IDL2016","combined"), each = length(endpoints))),
  power = c(powerIstat_ep, powerFrance_ep, powerIDL_ep, powerCombo_ep)
)

lifetimes <- data.frame(
  death = c(italcent$numdays[italcent$numdays>365.25*115],
            francent$numdays[francent$numdays>365.25*115],
            idl2016$numdays[idl2016$numdays>365.25*115]
  )/365.25,
  data = factor(c(rep("Istat", sum(italcent$numdays>365.25*115)),
                  rep("France", sum(francent$numdays>365.25*115)),
                  rep("IDL2016", sum(idl2016$numdays>365.25*115)))))

load("power.RData") 
nxis <- length(xis)
power_italcent[,xis == 0,]

# Because of the size distortion, we approximate the null distribution of the Wald
# statistic, LRT, etc. using the bootstrap values
critIstat <- apply(power_italcent[,xis == 0,], 2, quantile, c(0.025,0.05, 0.975))
# powerIstat <- colMeans(cbind(abs(power_italcent[,,2]) > qnorm(0.975),
#                              abs(power_italcent[,,1]) > qnorm(0.975),
#                              power_italcent[,,2] < qnorm(0.05),
#                              power_italcent[,,1] < qnorm(0.05)), na.rm = TRUE)
powerIstat <- colMeans(cbind((power_italcent[,,2] > critIstat[3,2])&(power_italcent[,,2] < critIstat[1,2]),
                             (power_italcent[,,1] > critIstat[3,1])&(power_italcent[,,1] < critIstat[1,1]),
                             power_italcent[,,2] < critIstat[2,2],
                             power_italcent[,,1] < critIstat[2,1]), na.rm = TRUE)
apply(power_italcent, c(2,3), function(x){sum(is.na(x))})
critFran <- apply(power_francent[,xis == 0,], 2, quantile, c(0.025,0.05, 0.975))
# powerfrancent <- colMeans(cbind(abs(power_francent[,,2]) > qnorm(0.975),
#                                 abs(power_francent[,,1]) > qnorm(0.975),
#                                 power_francent[,,2] < qnorm(0.05),
#                                 power_francent[,,1] < qnorm(0.05)), na.rm = TRUE)
powerfrancent <- colMeans(cbind((power_francent[,,2] > critFran[3,2])&(power_francent[,,2] < critFran[1,2]),
                                (power_francent[,,1] > critFran[3,1])&(power_francent[,,1] < critFran[1,1]),
                                power_francent[,,2] < critFran[2,2],
                                power_francent[,,1] < critFran[2,1]), na.rm = TRUE)
apply(power_francent, c(2,3), function(x){sum(is.na(x))})
critIDL <- apply(power_idl[,xis == 0,], 2, quantile, c(0.025,0.05, 0.975))
# powerIDL <- colMeans(cbind(abs(power_idl[,,2]) > qnorm(0.975),
#                            abs(power_idl[,,1]) > qnorm(0.975),
#                            power_idl[,,2] < qnorm(0.05),
#                            power_idl[,,1] < qnorm(0.05)), na.rm = TRUE)
apply(power_idl, c(2,3), function(x){sum(is.na(x))})
powerIDL <- colMeans(cbind((power_idl[,,2] > critIDL[3,2])&(power_idl[,,2] < critIDL[1,2]),
                           (power_idl[,,1] > critIDL[3,1])&(power_idl[,,1] < critIDL[1,1]),
                           power_idl[,,2] < critIDL[2,2],
                           power_idl[,,1] < critIDL[2,1]), na.rm = TRUE)

powerAnyWald <- 1-(1-powerIDL)*(1-powerfrancent)*(1-powerIstat)

power_df <- data.frame(
  shape = rep(xis, length.out = 16*nxis),
  test = factor(rep(rep(c("Wald","directed likelihood root"), each = nxis), length.out = 16*nxis)),
  data = factor(rep(c("Istat","France","IDL2016","combined"), each = 4*nxis)),
  hypothesis = factor(rep(rep(c("two-sided","one-sided"), each = 2*nxis), length.out = 16*nxis)),
  power = c(powerIstat, powerfrancent, powerIDL, powerCombo))

power_any <- data.frame(shape = rep(xis, length.out = 4*nxis),
                        test = factor(rep(rep(c("Wald","directed likelihood root"), each = nxis), length.out = 2*nxis)),
                        hypothesis = factor(rep(c("two-sided","one-sided"), each = 2*nxis)),
                        power = powerAny)
power_anyW <- power_any %>% filter((hypothesis == "one-sided")&(test == "Wald"))

g1 <- power_df %>% filter((hypothesis == "one-sided")&(test == "directed likelihood root")) %>%
  ggplot(aes(x = shape, y = power, col = data)) + 
  geom_hline(yintercept = 0.05, alpha = 0.5) + 
  geom_smooth(method = "gam", se = FALSE,
              formula = y ~ s(x, bs = "tp", fx = TRUE, k=25), show.legend = FALSE) +
  theme_classic() + 
  theme(panel.grid.major = element_line())  + #seq(0,1, by = 0.1)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.01),
                     breaks = seq(0,1, by= 0.25),
                     labels = c("$0$","$0.25$","$0.5$", "$0.75$", "$1$"))  + 
  scale_x_continuous(breaks = seq(-0.25,0, by = 0.05),
                     limits = c(-0.25,0),
                     labels = paste0("$",seq(-0.25,0, by = 0.05),"$")) + 
  xlab("$\\gamma$")
Vcols <- viridis(4)[c(2:4,1)]
g1b <- power_df %>% filter((hypothesis == "one-sided")&(test == "Wald")&(data != "combined")) %>%
  ggplot(aes(x = shape, y = power)) + 
  geom_hline(yintercept = 0.05, alpha = 0.5) + 
  geom_smooth(mapping = aes(col = data), method = "gam", se = FALSE,
              formula = y ~ s(x, bs = "tp", fx = TRUE, k=25), show.legend = FALSE) +
  geom_smooth(data = power_anyW, aes(x = shape, y = power), color = "black", lty = 2, method = "gam", se = FALSE,
              formula = y ~ s(x, bs = "tp", fx = TRUE, k=25), show.legend = FALSE) +
  scale_color_manual(values = Vcols) + 
  theme_classic() + 
  theme(panel.grid.major = element_line(), 
        legend.position = "bottom")  + #seq(0,1, by = 0.1)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.01),
                     breaks = seq(0,1, by= 0.25),
                     labels = c("$0$","$0.25$","$0.5$", "$0.75$", "$1$"))  + 
  scale_x_continuous(breaks = seq(-0.25,0, by = 0.05),
                     limits = c(-0.25,0),
                     labels = paste0("$",seq(-0.25,0, by = 0.05),"$")) + 
  xlab("$\\gamma$")

g2 <- ggplot(data = power_endp, aes(x = endpoint, y = power, col = data)) + 
  geom_hline(yintercept = 0.05, alpha = 0.5) + 
  geom_rug(data=lifetimes, mapping = aes(x = death, y = NULL, col = data)) + 
  geom_line(data = power_endp, 
            aes(x = endpoint, 
                y = power, col = data),
            size = 1) + 
  scale_color_viridis_d() + 
  theme_classic() + 
  theme(legend.position = "bottom",
        panel.grid.major = element_line())  + 
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.01),
                     breaks = seq(0,1, by= 0.25),
                     labels = c("$0$","$0.25$","$0.5$", "$0.75$", "$1$"))  + 
  scale_x_continuous(breaks = seq(120, 180, by = 10),
                     limits = c(min(endpoints),181),
                     expand = c(0,0),
                     labels = paste0("$",seq(120, 180, by = 10),"$")) + 
  xlab("upper limit to human lifespan") 


if(figures){
  fig <- "Fig9.tex"
  setwd(fig_dir)
  tikz(fig, width = 8, height = 4, standAlone = TRUE)
}
g2 + g1b + plot_layout(guides = 'collect') & theme(legend.position="bottom")
if(figures){
  dev.off()
  system(command = paste0("pdflatex ",fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd("..")
}


# Compute the power at different values
round(rbind(
  powersmoothIstat[which(endpts %in% c(125, 130, 135))],
  powersmoothFrance[which(endpts %in% c(125, 130, 135))],
  powersmoothIDL[which(endpts %in% c(125, 130, 135))],
  powersmoothCombo[which(endpts %in% c(125, 130, 135))]
),2)

#########################################################################################
#########################################################################################
# Extended generalized Pareto, Gompertz model (Tables 3 and 4)
source("Semi-supercentenarian_fn.R")
# This database must be purchased from Istat
load("italcent.rda")

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

# Sequence of thresholds, from 105 years (in days), increments of one year
thresh <- seq(from = 0, to = 2191.5, by = 365.25)


# Exponential model
param_exp <- matrix(0, nrow = length(thresh), ncol = 4L)
colnames(param_exp) <- c("loc","scale","deviance", "scale.stderror")
for(i in 1:length(thresh)){
  datu <- dat - thresh[i]
  ind <- which(datu > 0)
  vals <- list(par = exp_mle_lt_rc(dat = (datu[ind])/365.25, rightcens = rightcens[ind], slow =  (pmax(0, slow[ind]-thresh[i]))/365.25))
  vals$value <- gpd_cens(par = c(vals$par, 0), dat = (datu[ind])/365.25, rightcens = rightcens[ind], slow =  (pmax(0, slow[ind]-thresh[i]))/365.25)
  vals$hessian <- numDeriv::hessian(func = function(par){gpd_cens(par = par, dat = (datu[ind])/365.25, rightcens = rightcens[ind], slow =  (pmax(0, slow[ind]-thresh[i]))/365.25, expo = TRUE)}, x = vals$par)
  # vals <- optim(par = c(500), fn = gpd_cens, method = "Brent", lower = 0, upper = 800, dat = (datu[ind])/365.25,
  # rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25, expo = TRUE, hessian = TRUE)
  param_exp[i,] <- c(round((u + thresh[i])/365.25,2), vals$par, -2*vals$value, sqrt(solve(vals$hessian)))
}



# Fit an extended model
param_exggp <- matrix(0, nrow = length(thresh), ncol = 9L)
colnames(param_exggp) <- c("thresh","scale","shape1", "shape2","deviance", "scale.stderror", "shape1.stderror", "shape2.stderror", "nu")
for(i in 1:length(thresh)){
  #For each threshold, compute new threshold exceedances
  datu <- dat - thresh[i]
  # Keep only exceedances
  ind <- which(datu > 0)
  if(i == 1){
    start <- c(1.6,0.1,-0.05)
  } else{
    start <- param_exggp[i-1,2:4]
  }
  objfun <- function(par){
    exggp_cens(par = par, dat = (datu[ind])/365.25, 
               rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25)
  }
  # Compute maximum likelihood estimates numerically
  #   # Compute maximum likelihood estimates numerically
  test <- optim(par = start,
                fn = objfun, 
                method = "BFGS")
  vals <- nlm(p = test$par, 
              f = objfun,
              iterlim = 1e8,
              stepmax = 1, gradtol = 1e-8, hessian = TRUE, 
              typsize = start, 
              fscale = param_exp[i,"deviance"]/2)
  test <- nloptr::slsqp(x0 = vals$estimate, fn = objfun, control = list(xtol_rel = 1e-12))
  test$gradient <- numDeriv::grad(func = objfun, x =test$par)
  cat(round((u + thresh[i])/365.25,2),":", test$gradient,"\n")
  std <- try(sqrt(diag(solve(numDeriv::hessian(func = objfun, x = vals$estimate)))))
  if(is.character(std)){
    std <- rep(NA, 3)
  }
  param_exggp[i,] <- c(round((u + thresh[i])/365.25,2), 
                       test$par, -2*test$value, std, length(ind))
}

# Fit a generalized Pareto model
param_gpd <- matrix(0, nrow = length(thresh), ncol = 7L)
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
  # Bundle estimates in a vector
  param_gpd[i,] <- c(loc = (u + thresh[i])/365.25, vals$par, -2*vals$value, sqrt(diag(solve(vals$hessian))), length(ind))
}
colnames(param_gpd) <- c("loc","scale","shape","deviance", "scale.stderror", "shape.stderror", "nu")

# Fit Gompertz
param_gomp <- matrix(0, nrow = length(thresh), ncol = 7L)
colnames(param_gomp) <- c("thresh","scale","shape","deviance", "scale.stderror", "shape.stderror", "nu")
for(i in 1:length(thresh)){
  #For each threshold, compute new threshold exceedances
  datu <- dat - thresh[i]
  # Keep only exceedances
  ind <- which(datu > 0)
  if(i == 1){
    start <- c(1.6,0.1)
  } else{
    start <- param_gomp[i-1,2:3]
  }
  # Compute maximum likelihood estimates numerically
  #   # Compute maximum likelihood estimates numerically
  test <- optim(par = start, 
                fn = function(par){
                  exggp_cens(par = c(par,0), dat = (datu[ind])/365.25, 
                             rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25)}, method = "N")
  vals <- nlm(p = start, 
              f = function(par){
                exggp_cens(par = c(par,0), dat = (datu[ind])/365.25, 
                           rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25)},
              iterlim = 1e8,
              stepmax = 1, gradtol = 1e-8, hessian = TRUE, 
              typsize = start, 
              fscale = param_exp[i,"deviance"]/2)
  print(round((u + thresh[i])/365.25,2)); print(vals$gradient)
  # Bundle estimates in a vector
  std <- try(sqrt(diag(solve(numDeriv::hessian(func = function(par){exggp_cens(par = c(par,0), dat = (datu[ind])/365.25, rightcens = rightcens[ind], slow = (pmax(0, slow[ind]-thresh[i]))/365.25)}, x = vals$estimate)))))
  if(is.character(std)){
    std <- rep(NA, 2)
  }
  if(vals$code %in% 1:3){
    param_gomp[i,] <- c(round((u + thresh[i])/365.25,2), 
                        ifelse(rep(vals$minimum <test$value,2), vals$estimate, test$par),
                        -2*ifelse(vals$minimum <test$value, vals$minimum, test$value), 
                        std,
                        length(ind))
  }
}


# LRT weibull versus exponential (shape = 1)
param_weibull <- cbind(param_weibull, pvalexp = pchisq(param_weibull[,"deviance"] - param_exp[,"deviance"], 1, lower.tail = FALSE))

# Extended model versus exponential
round(pchisq(param_exggp[,"deviance"] - param_exp[,"deviance"], 2, lower.tail = FALSE),3)
round(pchisq(param_exggp[,"deviance"] - param_gpd[,"deviance"], 1, lower.tail = FALSE)/2,3)
round(pchisq(param_exggp[,"deviance"] - param_gomp[,"deviance"], 1, lower.tail = FALSE),3)

round(pchisq(param_gomp[,"deviance"] - param_exp[,"deviance"], 1, lower.tail = FALSE),3)
round(pchisq(param_gpd[,"deviance"] - param_exp[,"deviance"], 1, lower.tail = FALSE),3)

if(tables){
  tab_dev <- matrix(NA, ncol = 7, nrow = 3)
  ii <- 1:7
  tab_dev[1,] <- paste0("$", sprintf(param_exggp[ii,"deviance"] - param_gomp[ii,"deviance"], f = "%1.2f"), "$")
  tab_dev[2,] <- paste0("$", sprintf(param_exggp[ii,"deviance"] - param_gpd[ii,"deviance"], f = "%1.2f"), "$")
  tab_dev[3,] <- paste0("$", sprintf(param_exggp[ii,"deviance"] - param_exp[ii,"deviance"], f = "%1.2f"), "$")
  
  rownames(tab_dev) <- c("Gompertz", "gen. Pareto","exponential")
  xtab1 <- xtable::xtable(tab_dev, 
                          caption = "Likelihood ratio statistic (deviance) for the extended generalized Pareto for comparisons with the Gompertz, generalized Pareto and exponential sub-models, for different thresholds.",
                          align = paste0(rep("r",ncol(tab_dev)+1),collapse = ""))
  table_dir <- paste0(getwd(), "/tables")
  setwd(table_dir)
  addtorow <- list()
  addtorow$pos <- list(0)
  addtorow$command <-  paste0("threshold", paste0(" & $", 105:110, "$", collapse = " "), "\\\\\n", collapse = " ")
  print(xtab1, file = "Table4.tex", 
        size = "footnotesize", 
        floating.environment = "table*",  
        booktabs = TRUE, caption.placement = "top",
        add.to.row = addtorow, include.colnames = FALSE,
        sanitize.colnames.function = identity, 
        sanitize.text.function = identity, 
        sanitize.rownames.function = identity)
  setwd(code_dir)
}
# Gompertz model has a = shape*scale^(-shape) and b = shape-1 provided shape > 1
cbind(param_weibull[,"shape"]*param_weibull[,"scale"]^(-param_weibull[,"shape"]),
      param_weibull[,"shape"]-1)

mev::smith.penult(family = "weibull", method = "pot", u=thresh[-1]/365.25, 
                  shape = param_weibull[1,"shape"], scale = param_weibull[1,"scale"] )

mev::smith.penult(family = "weibull", method = "pot", u=thresh[-1]/365.25, 
                  shape = param_weibull[7,"shape"], scale = param_weibull[7,"scale"] )

# Compute the bootstrap p-value for the test of exponentiality for Gompertz with beta=0
# Repeatedly sample from exponential and compute the likelihood ratio statistic
# the sampling mimics that of the power study
if(save){
  B <- 1e4L
  set.seed(1245)
  slow <- as.numeric(pmax(0, c1 - xcal))
  boot_lrt_gomp_vs_exp <- matrix(0, nrow = B, ncol = length(thresh))
  for(j in 1:length(thresh)){
    #For each threshold, compute new threshold exceedances
    datu <- (dat - thresh[j])/365.25
    ind <- which(datu > 0)
    datu <- datu[ind]
    uplim <- as.numeric(c2 - italcent$birth)
    # Keep only exceedances above 108
    uplim <- (uplim[ind] - 365.25*105 - thresh[j])/365.25
    slowmu <- pmax(0,(slow[ind] - thresh[j])/365.25)
    # Keep only exceedances
    
    bootsamp <- matrix(0, nrow = length(ind), ncol = B)
    bootrcens <- matrix(FALSE, nrow = length(ind), ncol = B)
    for(i in 1:length(ind)){
      bootsamp[i,] <- qexp(pexp(slowmu[i], rate = 1/param_exp[j,"scale"]) + runif(B)*
                             (1-pexp(slowmu[i], rate = 1/param_exp[j,"scale"])), 
                           rate = 1/param_exp[j,"scale"])
      bootrcens[i,] <- bootsamp[i,] > uplim[i]
      bootsamp[i,bootrcens[i,]] <- uplim[i]
    }
    start <- param_gomp[j,2:3]
    for(b in 2:B){
      # Fit the Gompertz model on the bth bootstrap sample
      test <- optim(par = start, 
                    fn = function(par){
                      exggp_cens(par = c(par,0), dat = bootsamp[,b], 
                                 rightcens = bootrcens[,b], slow = slowmu)}, method = "N")
      vals <- nlm(p = test$par, 
                  f = function(par){
                    exggp_cens(par = c(par,0), dat = bootsamp[,b], 
                               rightcens = bootrcens[,b], slow = slowmu)},
                  iterlim = 1e8,
                  stepmax = 1, gradtol = 1e-8, hessian = TRUE, 
                  typsize = start, 
                  fscale = param_exp[j,"deviance"]/2)
      exp_mle <- exp_mle_lt_rc(dat = bootsamp[,b], rightcens = bootrcens[,b], slow = slowmu)
      boot_lrt_gomp_vs_exp[b, j] <- -2*(exggp_cens(par = c(vals$estimate,0), dat = bootsamp[,b], rightcens = bootrcens[,b], slow = slowmu) - 
                                          gpd_cens(par = exp_mle, dat = bootsamp[,b], rightcens = bootrcens[,b], slow = slowmu, expo = TRUE))
    }
  }
  boot_lrt_gomp_vs_exp[1,] <- pmax(0, param_gomp[,"deviance"] - param_exp[,"deviance"])
  pvalboot <- apply(pmax(boot_lrt_gomp_vs_exp,0), 2, function(x){mean(x>= x[1])})
  save(file = "pvalboot_Istat_gompvsexp.Rdata", pvalboot, thresh, boot_lrt_gomp_vs_exp)
  
} else{
  load(file = "pvalboot_Istat_gompvsexp.Rdata")
}

gomp_par <- data.frame(thresh = paste0(105+thresh/365.25), 
                       param_gomp[,1:3])
table.res <- matrix(nrow = 5, ncol = 7)
ind <- 1:length(thresh)
table.res[1,] <- paste0("$", round(param_gpd[ind,"nu"],0),"$")
table.res[2,] <- paste0("$", round(param_gomp[ind, "scale"],2), ifelse(is.na(param_gomp[ind, "scale.stderror"]), paste0("\\; (",round(param_exp[ind, "scale.stderror"],1),")$"), paste0("\\; (",round(param_gomp[ind, "scale.stderror"],2),")$")))
table.res[3,] <- paste0("$", round(param_gomp[ind, "shape"],2), ifelse(is.na(param_gomp[ind, "shape.stderror"]), "$", paste0("\\; (",round(param_gomp[ind, "shape.stderror"],2),")$")))
# table.res[4,] <- paste0("$", round(param_exp[ind, "scale"],2), "\\; (",round(param_exp[ind, "scale.stderror"],1),")$")
table.res[4,] <- paste0("$", round(0.5*(param_gomp[ind,"deviance"] <= param_exp[ind,"deviance"]) + pchisq(param_gomp[ind,"deviance"]- param_exp[ind,"deviance"],1, lower.tail = FALSE)/2, 2),"$")
table.res[5,] <- paste0("$", round(pvalboot,2), "$")
rownames(table.res) <- c("$n_u$", "$\\sigma$", "$\\beta$", "$p$-value (asymptotic)", "$p$-value (bootstrap)")
xtab2 <- xtable::xtable(table.res, 
                        caption = "Estimates (standard errors) of Gompertz parameters ($\\beta$, $\\sigma$) for the Italian \\Istat{} data as a function of threshold, with number of threshold exceedances ($n_u$), $p$-value for the likelihood ratio test of $\\beta=0$. Estimates reported as zero for $\\beta$ are smaller than $10^{-7}$.",
                        align = paste0(rep("r",ncol(table.res)+1),collapse = ""))
addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <-  paste0("threshold", paste0(" & $", 105:111, "$", collapse = " "), "\\\\\n", collapse = " ")
print(xtab2, file = "tables/Table3.tex", 
      size = "footnotesize", 
      floating.environment = "table*",  
      booktabs = TRUE, caption.placement = "top",
      add.to.row = addtorow, include.colnames = FALSE,
      sanitize.colnames.function = identity, 
      sanitize.text.function = identity, 
      sanitize.rownames.function = identity)

