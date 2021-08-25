########################################################################
##   This script must be called from within Semi-supercentenarian.R   ##
########################################################################

# Graphical goodness-of-fit diagnostic plots
# for the French data
load("francent.rda")
u108 <- 39447L
xcal <- francent$birth + u108
# Calendar time for (common) sampling frame
c1 <- lubridate::dmy("01-01-1987")
c2 <- lubridate::dmy("31-12-2017")

# Lower truncation level, zero if individual reached 108 between c1 and c2
surv108 <- francent$numdays > u108
etime108 <- francent$numdays[surv108] - u108
ltrunc108 <- as.numeric(pmax(0, c1 - xcal[surv108]))
rtrunc108 <- as.numeric(c2 - xcal[surv108])



# Fit the nonparametric maximum likelihood estimator

# Parametric models for threshold exceedances above 108
# Fit exponential model
exp_fit_fr <- longevity::prof_exp_scale(
  time = etime108/365.25,
  thresh = 0,
  ltrunc = ltrunc108/365.25,
  rtrunc = rtrunc108/365.25)

# Compute pointwise 95% confidence intervals 
# for the cumulative hazard
t_grid <- seq(0.01, 15, length.out = 101)
cumhaz_gp_fr <- matrix(nrow = length(t_grid), ncol = 3)
for(i in seq_along(t_grid)){
  cumhaz_fit_gp_fr <- try(cumhaz_gp(x = t_grid[i],
                          time = etime108/365.25,
                          thresh = 0,
                          ltrunc = ltrunc108/365.25,
                          rtrunc = rtrunc108/365.25))
   if(!is.character(cumhaz_fit_gp_fr)){
    cumhaz_gp_fr[i,] <- c(cumhaz_fit_gp_fr$par, 
                          cumhaz_fit_gp_fr$confint)
  }
}
# Fix outliers due to lack of convergence of software routine
cumhaz_gp_fr[,2] <- 
  fitted(cobs::cobs(x = t_grid, 
                     y = cumhaz_gp_fr[,2], 
                     constraint = "increase"))
cumhaz_gp_fr[,3] <- 
  fitted(cobs::cobs(x = t_grid, 
                    y = cumhaz_gp_fr[,3], 
                    constraint = "increase"))
# Fit the exponential and generalized 
# Pareto distributions to exceedances
#  above 108
fit_exp_fr <- fit_elife(time = etime108/365.25,
                     ltrunc = ltrunc108/365.25,
                     rtrunc = rtrunc108/365.25,
                     family = "exp",
                     export = TRUE)
fit_gp_fr <-
  fit_elife(time = etime108/365.25,
            ltrunc = ltrunc108/365.25,
            rtrunc = rtrunc108/365.25,
            family = "gp",
            export = TRUE)



# Nonparametric maximum likelihood estimation
# of the survival function and estimation of the
# cumulative hazard with equal precision bands
npsurv_fr <- longevity::np_elife(
  time = etime108,
  ltrunc = ltrunc108,
  rtrunc = rtrunc108,
  vcov = TRUE)

utimes_fr <- npsurv_fr$xval/365.25
# Convert the object to cumulative hazard
cumhaz_fr <- -log(npsurv_fr$surv)
# Convert the covariance matrix of jump to 
# cumulative hazard using the delta-method
jacobian <- matrix(0,
                   nrow = length(npsurv_fr$surv),
                   ncol = length(npsurv_fr$surv))
for(i in seq_len(nrow(jacobian))){
  for(j in seq_len(i)){
    jacobian[i,j] <- 1/(1-sum(npsurv_fr$prob[1:i]))
  }
}
cumhaz_std_fr <- sqrt(diag(jacobian %*% npsurv_fr$vcov %*% t(jacobian)))
# Equal precision bounds over whole range (log transformed)
ep_critical <- function(d, cumhazstd, tL, tU){
  4*dnorm(d)/d+2*dnorm(d)*(d-1/d)*
    (log(cumhazstd[tU])-log(cumhazstd[tL])) - 0.05
}
# Critical value
dalpha <- uniroot(ep_critical,
                  cumhazstd = cumhaz_std_fr,
                  tL = 1,
                  tU = length(cumhaz_std_fr),
                  interval = c(2, 5))$root
# Compute log transformed bands
cumhaz_fr.log.confint <- cbind(
  cumhaz_fr*exp(-dalpha*cumhaz_std_fr/cumhaz_fr),
  cumhaz_fr*exp(dalpha*cumhaz_std_fr/cumhaz_fr))


# Repeat the plots, this time for the IDL data
# Load IDL data
load("IDL2021.rda")
# Exponential MLE with 95% confidence intervals
exp_fit_idl <- 
  with(idlex,
       longevity::prof_exp_scale(
         time = datu, # exceedances (day) above 110
         thresh = 0, # threshold already subtracted
         ltrunc = slow, #left truncation
         rtrunc = supp)) #right truncation)

fit_exp_idl <- with(idlex,
     longevity::fit_elife(time = datu,
                          ltrunc = slow, 
                          rtrunc = supp,
                          thresh = 0,
                          family = "exp",
                          export = TRUE))

fit_gp_idl <- 
  with(idlex,
      longevity::fit_elife(time = datu,
                           ltrunc = slow, 
                           rtrunc = supp,
                           thresh = 0,
                           family = "gp",
                           export = TRUE))


# Nonparametric maximum likelihood estimation
# of the survival function and estimation of the
# cumulative hazard with equal precision bands
npsurv_idl <- 
  with(idlex,
  longevity::np_elife(
  time = datu,
  ltrunc = slow,
  rtrunc = supp,
  thresh = 0,
  vcov = TRUE))


cumhaz_gp_idl <- matrix(nrow = length(t_grid), 
                        ncol = 3)
for(i in seq_along(t_grid)){
    cumhaz_fit_gp_idl <- 
      with(idlex,
           try(cumhaz_gp(x = t_grid[i],
                         time = datu,
                         thresh = 0,
                         ltrunc = slow,
                         rtrunc = supp)))
  if(!is.character(cumhaz_fit_gp_idl)){
    cumhaz_gp_idl[i,] <- c(cumhaz_fit_gp_idl$par, 
                          cumhaz_fit_gp_idl$confint)
  }
}
# Poor man fix: some outliers due to problems 
# with the numerical optimization routine
cumhaz_gp_idl[,2] <- fitted(cobs::cobs(constraint = "none",
                                       x = t_grid, 
                                       y = cumhaz_gp_idl[,2]))

utimes_idl <- npsurv_idl$xval
# Convert the object to cumulative hazard
cumhaz_idl <- -log(npsurv_idl$surv)
# Convert the covariance matrix of jump to 
# cumulative hazard using the delta-method
jacobian <- matrix(0,
                   nrow = length(npsurv_idl$surv),
                   ncol = length(npsurv_idl$surv))
for(i in seq_len(nrow(jacobian))){
  for(j in seq_len(i)){
    jacobian[i,j] <- 1/(1-sum(npsurv_idl$prob[1:i]))
  }
}
cumhaz_std_idl <- sqrt(diag(jacobian %*% npsurv_idl$vcov %*% t(jacobian)))
# Critical value
dalpha <- uniroot(ep_critical,
                  cumhazstd = cumhaz_std_idl,
                  tL = 1,
                  tU = length(cumhaz_std_idl),
                  interval = c(1, 7))$root
# Compute log transformed bands
cumhaz_idl.log.confint <- cbind(
  cumhaz_idl*exp(-dalpha*cumhaz_std_idl/cumhaz_idl),
  cumhaz_idl*exp(dalpha*cumhaz_std_idl/cumhaz_idl))



# Create plots
if(figures){
  setwd(fig_dir)
  fig <- "Fig10.tex"
  tikz(fig, 
       width = 1.42*dwidth, 
       height = 1.2*dheight,
       standAlone = TRUE)
}
par(mar = c(4,4,2,2),
    mfrow = c(1,2),
    cex = 1,
    bty = "l")
# Cumulative hazard
# repeat last element to create the step
plot(x = utimes_fr,
     y = c(cumhaz_fr, cumhaz_fr[length(cumhaz_fr)]),
     bty = "l",
     type = "s",
     lwd = 2,
     ylab = "",
     xlab = "",
     xaxs = "i",
     xlim = c(0, 14.5),
     xaxt = 'n',
     ylim = c(0, 14.5),
     yaxs = "i")
axis(side = 1,
     at = seq(0L, 16L, by = 2L),
     labels = paste("$",
                    seq(0L, 16L, by = 2L),
                    "$"))
mtext(text = "excess lifetime above 108",
      side = 1,
      line = 2)
mtext(text = "conditional cumulative hazard",
      side = 2,
      line = 2)
lines(x = utimes_fr,
      y = c(cumhaz_fr.log.confint[,1], 
            tail(cumhaz_fr.log.confint[,1],1)),
      lty = 2,
      type = "s")
lines(x = utimes_fr,
      y = c(cumhaz_fr.log.confint[,2],
            tail(cumhaz_fr.log.confint[,2],1)),
      lty = 2,
      type = "s")
polygon(x = c(0,16,16,0),
        y = c(0,16/exp_fit_fr[2], 16/exp_fit_fr[3],0),
        col = scales::alpha("blue", alpha = 0.1),
        border = NA)
abline(a = 0, b = 1/exp_fit_fr[1], col = 4)
abline(a = 0, b = 1/exp_fit_fr[2], col = 4, lty = 2)
abline(a = 0, b = 1/exp_fit_fr[3], col = 4, lty = 2)

# Generalized Pareto
polygon(x = c(t_grid, rev(t_grid)),
        y = c(cumhaz_gp_fr[,2], rev(cumhaz_gp_fr[,3])),
        col = scales::alpha("red", alpha = 0.1),
        border = NA)
lines(t_grid, cumhaz_gp_fr[,1], col = 2) #point estimates
lines(t_grid, cumhaz_gp_fr[,2], col = 2, lty = 2) #lower confint
lines(t_grid, cumhaz_gp_fr[,3], col = 2, lty = 2) #upper confint


plot(x = utimes_idl,
     y = c(cumhaz_idl, tail(cumhaz_idl, 1)),
     bty = "l",
     type = "s",
     lwd = 2,
     ylab = "",
     xlab = "",
     xaxs = "i",
     xlim = c(0, 8),
     xaxt = 'n',
     ylim = c(0, 8),
     yaxs = "i")
axis(side = 1,
     at = seq(0L, 16L, by = 2L),
     labels = paste("$",
                    seq(0L, 16L, by = 2L),
                    "$"))
mtext(text = "excess lifetime above 110",
      side = 1,
      line = 2)
mtext(text = "conditional cumulative hazard",
      side = 2,
      line = 2)
lines(x = utimes_idl,
      y = c(cumhaz_idl.log.confint[,1],
            tail(cumhaz_idl.log.confint[,1], 1)),
      lty = 2,
      type = "s")
lines(x = utimes_idl,
      y = c(cumhaz_idl.log.confint[,2],
            tail(cumhaz_idl.log.confint[,2], 1)),
      lty = 2,
      type = "s")
polygon(x = c(0,13,13,0),
        y = c(0,13/exp_fit_idl[2], 13/exp_fit_idl[3],0),
        col = scales::alpha("blue", alpha = 0.1),
        border = NA)
abline(a = 0, b = 1/exp_fit_idl[1], col = 4)
abline(a = 0, b = 1/exp_fit_idl[2], col = 4, lty = 2)
abline(a = 0, b = 1/exp_fit_idl[3], col = 4, lty = 2)

# Generalized Pareto
polygon(x = c(t_grid, rev(t_grid)),
        y = c(cumhaz_gp_idl[,2], rev(cumhaz_gp_idl[,3])),
        col = scales::alpha("red", alpha = 0.1),
        border = NA)
lines(t_grid, cumhaz_gp_idl[,1], col = 2) #point estimates
lines(t_grid, cumhaz_gp_idl[,2], col = 2, lty = 2) #lower confint
lines(t_grid, cumhaz_gp_idl[,3], col = 2, lty = 2) #upper confint

if(figures){
  dev.off()
  system(command = paste0("lualatex ",
                          fig_dir, 
                          "/", 
                          fig,
                          "; rm *.aux; rm *.log"))
  setwd(code_dir)
}
