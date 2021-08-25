########################################################################
##   This script must be called from within Semi-supercentenarian.R   ##
########################################################################

########################################################################
#  Fig 6: Exponential quantile-quantile plot with parametric bootstrap and cumulative hazard plot
########################################################################

# Calendar time for sampling frame
# # 105 years threshold (criterion for inclusion in dataset)
u108 <- 39447L
load("italcent.rda")
c1 <- lubridate::dmy("01-01-2009")
c2 <- lubridate::dmy("01-01-2016")
ind108 <- which((italcent$numdays - u108) > 0)
slow108 <-  as.numeric(pmax(c1 - italcent$birth[ind108] - u108, 0))
numd108 <- italcent$numdays[ind108] - u108
uplim108 <- c2 - italcent$birth[ind108]
cens108 <- italcent$rightcens[ind108]
italcent_surv_108 <- kaplan.meier(time = numd108,
                                  lefttrunc = slow108,
                                  rightcens = cens108)
# Compute maximum likelihood estimator
param_exp108  <- 
  exp_mle_lt_rc(dat = numd108, 
                slow = slow108, 
                rightcens = italcent$rightcens[ind108])

B <- 10000L
set.seed(1234)
# We condition on the censoring / truncation values
bootsamp <- matrix(0, ncol = length(numd108), nrow = B)
for(i in 1:ncol(bootsamp)){
  if(cens108[i]){
    # the censored records are left untouched.
    bootsamp[,i] <- numd108[i]
  } else{
    # Due to censoring, individuals that are not censored are doubly truncated
    bootsamp[,i] <- qexp(pexp(slow108[i], rate = 1/param_exp108) +
                           runif(B)*(pexp(uplim108[i], rate = 1/param_exp108) - pexp(slow108[i], rate = 1/param_exp108)), rate = 1/param_exp108)
  }
}

qqpts <- qqptsb <- matrix(0, nrow = B, ncol = sum(!cens108))
thetab <- rep(0, length.out = B)
for(j in 1:B){
  thetab[j] <- exp_mle_lt_rc(dat = bootsamp[j,],
                             rightcens = cens108, 
                             slow = slow108)
  
  #Fit a Kaplan-Meier for left-truncated right-censored data
  KM.b <- kaplan.meier(time = bootsamp[j,],
                       lefttrunc = slow108,
                       rightcens = cens108)
  
  #Extract the timing and make into a empirical distribution function (weighted)
  cecdf <- mev:::.wecdf(KM.b$time, 
                        -diff(c(KM.b$surv, 0)))
  # Evaluate on a regular grid to get confidence intervals
  qqpts[j,] <- sort(qexp(cecdf(seq(0, 2800, length = sum(!cens108))),
                         rate = 1/thetab[j]))
  # Evaluate at simulated data
  qqptsb[j,] <- sort(qexp(cecdf(bootsamp[j,][!cens108]), 
                          rate = 1/thetab[j]))
}
surv_obj <-  kaplan.meier(time = numd108/365.25,
                          lefttrunc = slow108/365.25, 
                          rightcens = cens108)

#Extract the timing and make into a empirical distribution function (weighted)
cecdf <- mev:::.wecdf(surv_obj$time, 
                      -diff(c(surv_obj$surv, 0)))
# Plotting positions (x-axis)
xpos <- (cecdf(numd108[!cens108]/365.25) - cecdf(slow108[!cens108]/365.25))/(1-cecdf(slow108[!cens108]/365.25))
Fa <- pexp(slow108[!cens108],rate = 1/param_exp108)
txpos <- qexp(Fa + (1-Fa)*xpos, rate = 1/param_exp108)/365.25 + 108
env <- boot::envelope(mat = qqpts) 

# Fit nonparametric cumulative hazard via Nelson--Aalen
npfit <- nelson.aalen(time = numd108/365.25,
                      lefttrunc = slow108/365.25, 
                      rightcens = cens108)

# Log transformed bands
cumhaz.log.confint <- npfit$confband

# Fit exponential model
exp_fit <- longevity::prof_exp_scale(
         time = numd108/365.25,
         event = as.integer(!cens108),
         thresh = 0,
         ltrunc = slow108/365.25,
         type = "right")

# The cumulative hazard of the Generalized Pareto distribution
# is H(t) = log(1+t*xi/sigma)/xi

# Compute 95% pointwise profile confidence intervals
t_grid <- seq(0.01, 7, length.out = 101)
cumhazgp <- matrix(nrow = length(t_grid), ncol = 3)
for(i in seq_along(t_grid)){
  gp_fit <- try(cumhaz_gp(x = t_grid[i], 
                          time = numd108/365.25,
                          event = as.integer(!cens108),
                          thresh = 0,
                          ltrunc = slow108/365.25,
                          type = "right"))
  if(!is.character(gp_fit)){
    cumhazgp[i,] <- c(gp_fit$par, gp_fit$confint)
  }
}



if(figures){
  setwd(fig_dir)
  fig <- "Fig6.tex"
  tikz(fig,
       width = 1.42*dwidth, 
       height = 1.2*dheight, 
       standAlone = TRUE)
}
par(mar = c(4,4,0.4,0.4), 
     mfrow = c(1,2), 
     cex = 1, 
     bty = "l")
plot(y = numd108[!cens108]/365.25+108, 
     x = txpos,
     pch = 20, 
     bty = "l", 
     xlab = "", 
     ylab = "",
     panel.first = {abline(a=0, b=1)}, 
     yaxs = "i", 
     xaxs = "i",
     ylim = 108+c(0, 3100)/365.25, 
     xlim = 108+c(0, 3100)/365.25, 
     cex = 0.8)
mtext(text = "theoretical quantiles",
      side = 1,
      line = 2)
mtext(text = "empirical quantiles",
      side = 2,
      line = 2)

ss <- seq(0, 2800, length = sum(!cens108))/365.25+108
lines(y = ss, 
      x = env$point[1,]/365.25+108, 
      lty = 2, 
      lwd = 1.5)
lines(y = ss, 
      x = env$point[2,]/365.25+108, 
      lty = 2, 
      lwd = 1.5)
lines(y = ss, 
      x = env$overall[1,]/365.25+108, 
      lty = 3, 
      lwd = 1.5)
lines(y = ss, 
      x = env$overall[2,]/365.25+108, 
      lty = 3, 
      lwd = 1.5)


# Plot the cumulative hazard function along with fitted model
plot(x = c(npfit$time, 1e10),
     y = c(npfit$cumhaz, max(npfit$cumhaz)),
     bty = "l",
     xaxs = "i",
     yaxs = "i",
     xaxt = "n",
     lwd = 2,
     xlim = c(0,7),
     ylim = c(0,7),
     yaxt = 'n',
     type = "s",
     xlab = "",
     ylab = "")
axis(side = 1,
     at = 0:7,
     labels = c("$0$","","$2$","","$4$","","$6$",""))
axis(side = 2,
     at = 0:7,
     labels = c("$0$","","$2$","","$4$","","$6$",""))
mtext(text = "excess lifetime above 108",
      side = 1,
      line = 2)
mtext(text = "conditional cumulative hazard",
      side = 2,
      line = 2)
lines(x = npfit$time, 
      y = cumhaz.log.confint[,1],
      lty = 2,
      type = "s")
lines(x = npfit$time, 
      y = cumhaz.log.confint[,2],
      lty = 2,
      type = "s")
# Add exponential pointwise confidence intervals
# Because there is a single parameter and the cumulative
# hazard is linear in the parameter, we need only compute
# the profile likelihood and appeal to invariance of maximum
# likelihood estimators for the curves
polygon(x = c(0,7,7), 
        y = c(0, 7/exp_fit[2], 7/exp_fit[3]),
        col = scales::alpha("blue", alpha = 0.1),
        border = NA)
abline(a = 0, b = 1/exp_fit[1], col = 4)
abline(a = 0, b = 1/exp_fit[2], col = 4, lty = 2)
abline(a = 0, b = 1/exp_fit[3], col = 4, lty = 2)

# Generalized Pareto 
polygon(x = c(t_grid, rev(t_grid)), 
        y = c(cumhazgp[,2], rev(cumhazgp[,3])),
        col = scales::alpha("red", alpha = 0.1),
        border = NA)
lines(t_grid, cumhazgp[,1], col = 2) #point estimates
lines(t_grid, cumhazgp[,2], col = 2, lty = 2) #lower confint
lines(t_grid, cumhazgp[,3], col = 2, lty = 2) #upper confint

if(figures){
  dev.off()
  system(command = paste0("pdflatex ",fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}
