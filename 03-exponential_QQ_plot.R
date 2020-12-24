########################################################################
##   This script must be called from within Semi-supercentenarian.R   ##
########################################################################

########################################################################
#  Fig 6: Exponential quantile-quantile plot with parametric bootstrap
########################################################################

load("italcent.rda")
# 105 years threshold (criterion for inclusion in dataset
u <- 38351L
# Calendar date at which individual reaches 105 years
xcal <- italcent$birth + u
# Calendar time for sampling frame
c1 <- lubridate::dmy("01-01-2009")
c2 <- lubridate::dmy("01-01-2016")

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
