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
exp_fit <- longevity::prof_exp_scale(
  time = etime108/365.25,
  thresh = 0,
  ltrunc = ltrunc108/365.25,
  rtrunc = rtrunc108/365.25)

# Compute pointwise 95% confidence intervals 
# for the cumulative hazard
t_grid <- seq(0.01, 15, length.out = 101)
cumhazgp <- matrix(nrow = length(t_grid), ncol = 3)
for(i in seq_along(t_grid)){
  gp_fit <- try(cumhaz_gp(x = t_grid[i],
                          time = etime108/365.25,
                          thresh = 0,
                          ltrunc = ltrunc108/365.25,
                          rtrunc = rtrunc108/365.25))
  if(!is.character(gp_fit)){
    cumhazgp[i,] <- c(gp_fit$par, gp_fit$confint)
  }
}


# Fit the exponential and generalized 
# Pareto distributions to exceedances
#  above 108
fit_exp <- fit_elife(time = etime108/365.25,
                     ltrunc = ltrunc108/365.25,
                     rtrunc = rtrunc108/365.25,
                     family = "exp",
                     export = TRUE)
fit_gp <-
  fit_elife(time = etime108/365.25,
            ltrunc = ltrunc108/365.25,
            rtrunc = rtrunc108/365.25,
            family = "gp",
            export = TRUE)

B <- 10000L
xgrid <- seq(0, 14, length.out = 201)
qqpts_gp <- matrix(0, 
                nrow = B, 
                ncol = length(xgrid))
set.seed(1234)
for(b in seq_len(B)){
  bootsamp <- with(fit_gp,
    longevity::r_dtrunc_elife(
      n = nexc,
      scale = par['scale'],
      shape = par['shape'],
      lower = ltrunc,
      upper = rtrunc,
      family = "gp"))
  boot_np <- longevity::npsurv(time = bootsamp,
                            ltrunc = fit_gp$ltrunc,
                            rtrunc = fit_gp$rtrunc)
  boot_fit_gp <- longevity::fit_elife(
    time = bootsamp,
    ltrunc = fit_gp$ltrunc,
    rtrunc = fit_gp$rtrunc,
    family = "gp")
  qqpts_gp[b,] <- 
    longevity::qelife(
      p = fit_gp$nexc/(fit_gp$nexc + 1L)* boot_np$cdf(xgrid),
      scale = boot_fit_gp$par['scale'], 
      shape = boot_fit_gp$par['shape'],
      family = "gp")
}

qqpts_exp <- matrix(0, 
                nrow = B, 
                ncol = length(xgrid))
set.seed(1234)
for(b in seq_len(B)){
  bootsamp <- with(fit_exp,
                   longevity::r_dtrunc_elife(
                     n = nexc,
                     scale = par['scale'],
                     lower = ltrunc,
                     upper = rtrunc,
                     family = "exp"))
  boot_np <- longevity::npsurv(time = bootsamp,
                               ltrunc = fit_gp$ltrunc,
                               rtrunc = fit_gp$rtrunc)
  boot_fit_exp <- longevity::fit_elife(
    time = bootsamp,
    ltrunc = fit_exp$ltrunc,
    rtrunc = fit_exp$rtrunc,
    family = "exp")
  qqpts_exp[b,] <- 
    longevity::qelife(
      p = fit_exp$nexc/(fit_exp$nexc + 1L)* boot_np$cdf(xgrid),
      scale = boot_fit_exp$par['scale'], 
      family = "exp")
}
# Compute simulation envelopes
envel_gp <- boot::envelope(mat = qqpts_gp)
envel_exp <- boot::envelope(mat = qqpts_exp)
# Use this function to calculate the plotting positions
plots_gp <- plot(fit_gp,
              plot.type = "ggplot",
              which.plot = c("pp","qq"),
              plot = FALSE)

plots_exp <- plot(fit_exp,
                  plot.type = "ggplot",
                  which.plot = c("pp","qq"),
                  plot = FALSE)


# Nonparametric maximum likelihood estimation
# of the survival function and estimation of the
# cumulative hazard with equal precision bands
npsurv_fr <- longevity::np_elife(
  time = etime108,
  ltrunc = ltrunc108,
  rtrunc = rtrunc108,
  vcov = TRUE)

utimes <- npsurv_fr$xval[-length(npsurv_fr$xval)]/365.25
# Convert the object to cumulative hazard
cumhaz <- -log(npsurv_fr$surv)
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
cumhaz.std <- sqrt(diag(jacobian %*% npsurv_fr$vcov %*% t(jacobian)))
# Equal precision bounds over whole range (log transformed)
ep_critical <- function(d, cumhazstd, tL, tU){
  4*dnorm(d)/d+2*dnorm(d)*(d-1/d)*
    (log(cumhazstd[tU])-log(cumhazstd[tL])) - 0.05
}
# Critical value
dalpha <- uniroot(ep_critical,
                  cumhazstd = cumhaz.std,
                  tL = 1,
                  tU = length(cumhaz.std),
                  interval = c(2, 5))$root
# Compute log transformed bands
cumhaz.log.confint <- cbind(
  cumhaz*exp(-dalpha*cumhaz.std/cumhaz),
  cumhaz*exp(dalpha*cumhaz.std/cumhaz))


# Create plots
if(figures){
  setwd(fig_dir)
  fig <- "Fig10.tex"
  tikz(fig, width = 1.2*dwidth,
       height = 1.2*dheight,
       standAlone = TRUE)
}

par(mar = c(4,4,0.4,0.6),
    mfrow = c(1,2),
    cex = 1,
    bty = "l")


# Base R QQ-plot
plot(x = plots_exp$qq$data[,"x"],
     y = plots_exp$qq$data[,"y"],
     xlim = c(0, 15),
     ylim = c(0, 15),
     xlab = "theoretical quantiles",
     ylab = "empirical quantiles",
     bty = "l",
     pch = 19,
     type = "n",
     axes = FALSE,
     yaxs = "i",
     xaxs = "i",
     panel.first = {
       abline(a = 0, b = 1, col = "grey10")})
axis(side = 1,
     at = seq(0L, 16L, by = 2L),
     labels = paste("$",
                    seq(108L, 124L, by = 2L),
                    "$"))
axis(side = 2,
     at = seq(0L, 16L, by = 2L),
     labels = paste("$",
                    seq(108L, 124L, by = 2L),
                    "$"))
lines(y = envel_exp$point[1,],
      x = xgrid,
      lty = 2)
lines(y = envel_exp$point[2,],
      x = xgrid,
      lty = 2)
lines(y = envel_exp$overall[1,],
      x = xgrid,
      lty = 3)
lines(y = envel_exp$overall[2,],
      x = xgrid,
      lty = 3)
points(x = plots_exp$qq$data[,"x"],
       y = plots_exp$qq$data[,"y"],
       pch = 20,
       col = 1)
# Base R QQ-plot
plot(x = plots_gp$qq$data[,"x"],
     y = plots_gp$qq$data[,"y"],
     xlim = c(0, 15),
     ylim = c(0, 15),
     xlab = "theoretical quantiles",
     ylab = "empirical quantiles",
     bty = "l",
     pch = 19,
     type = "n",
     axes = FALSE,
     yaxs = "i",
     xaxs = "i",
     panel.first = {
       abline(a = 0, b = 1, col = "grey10")})
axis(side = 1,
     at = seq(0L, 16L, by = 2L),
     labels = paste("$",
                    seq(108L, 124L, by = 2L),
                    "$"))
axis(side = 2,
     at = seq(0L, 16L, by = 2L),
     labels = paste("$",
                    seq(108L, 124L, by = 2L),
                    "$"))
lines(y = envel_gp$point[1,],
      x = xgrid,
      lty = 2)
lines(y = envel_gp$point[2,],
      x = xgrid,
      lty = 2)
lines(y = envel_gp$overall[1,],
      x = xgrid,
      lty = 3)
lines(y = envel_gp$overall[2,],
      x = xgrid,
      lty = 3)
# for(i in 1:100){
#  lines(x = boot_qq[i,],
#        y = fit_gp$time,
#        lty = 2,
#        col = scales::alpha(colour = "black", 0.05))
# }
# points(x = plots$qq$data[,"x"],
#        y = plots$qq$data[,"y"],
#        pch = 1,
#        col = "white",
#        cex = 1)
points(x = plots_gp$qq$data[,"x"],
       y = plots_gp$qq$data[,"y"],
       pch = 20,
       col = 1)


# lines(x = envel_gp$overall[2,order(fit_gp$time)],
#       y = sort(fit_gp$time),
#       lty = 3)
# lines(x = envel_gp$overall[1,order(fit_gp$time)],
#       y = sort(fit_gp$time),
#       lty = 3)
#       


if(figures){
  dev.off()
  system(command = paste0("lualatex ",fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}

# Load IDL data
load("idl2021.rda")


fit_exp_idl <- with(idlex,
     longevity::fit_elife(time = datu,
                          ltrunc = slow, 
                          rtrunc = supp,
                          thresh = 0,
                          family = "exp",
                          export = TRUE))




xgrid_idl <- seq(0, 10, length.out = 201)
qqpts_exp_idl <- matrix(0, 
                    nrow = B, 
                    ncol = length(xgrid_idl))
set.seed(1234)
for(b in seq_len(B)){
  bootsamp <- with(fit_exp_idl,
                   longevity::r_dtrunc_elife(
                     n = nexc,
                     scale = par['scale'],
                     lower = ltrunc,
                     upper = rtrunc,
                     family = "exp"))
  boot_np <- longevity::npsurv(time = bootsamp,
                               ltrunc = fit_exp_idl$ltrunc,
                               rtrunc = fit_exp_idl$rtrunc)
  boot_fit_exp <- longevity::fit_elife(
    time = bootsamp,
    ltrunc = fit_exp_idl$ltrunc,
    rtrunc = fit_exp_idl$rtrunc,
    family = "exp")
  qqpts_exp_idl[b,] <- 
    longevity::qelife(
      p = fit_exp_idl$nexc/(fit_exp_idl$nexc + 1L)* boot_np$cdf(xgrid_idl),
      scale = boot_fit_exp$par['scale'], 
      family = "exp")
}
envel_exp_idl <- boot::envelope(mat = qqpts_exp_idl)

plots_exp_idl <- plot(fit_exp_idl,
                 plot.type = "ggplot",
                 which.plot = c("pp","qq"),
                 plot = FALSE)
# Create plots
if(figures){
  setwd(fig_dir)
  fig <- "Fig11.tex"
  tikz(fig, width = 1.2*dwidth,
       height = 1.2*dheight,
       standAlone = TRUE)
}

par(mar = c(4,4,0.4,0.6),
    mfrow = c(1,2),
    cex = 1,
    bty = "l")



plot(x = plots_exp_idl$qq$data[,"x"],
     y = plots_exp_idl$qq$data[,"y"],
     xlim = c(0, 9.5),
     ylim = c(0, 9.5),
     xlab = "theoretical quantiles",
     ylab = "empirical quantiles",
     bty = "l",
     pch = 19,
     type = "n",
     axes = FALSE,
     yaxs = "i",
     xaxs = "i",
     panel.first = {
       abline(a = 0, b = 1, col = "grey10")})
axis(side = 1,
     at = seq(0L, 10L, by = 2L),
     labels = paste("$",
                    seq(108L, 118L, by = 2L),
                    "$"))
axis(side = 2,
     at = seq(0L, 10L, by = 2L),
     labels = paste("$",
                    seq(108L, 118L, by = 2L),
                    "$"))
lines(y = envel_exp_idl$point[1,],
      x = xgrid_idl,
      lty = 2)
lines(y = envel_exp_idl$point[2,],
      x = xgrid_idl,
      lty = 2)
lines(y = envel_exp_idl$overall[1,],
      x = xgrid_idl,
      lty = 3)
lines(y = envel_exp_idl$overall[2,],
      x = xgrid_idl,
      lty = 3)
points(x = plots_exp_idl$qq$data[,"x"],
       y = plots_exp_idl$qq$data[,"y"],
       pch = 20,
       col = 1)


# Plot of cumulative hazard
plot(x = utimes,
     y = cumhaz,
     bty = "l",
     type = "s",
     lwd = 2,
     ylab = "conditional cumulative hazard",
     xlab = "excess lifetime (years)",
     xaxs = "i",
     xlim = c(0, 8),
     xaxt = 'n',
     ylim = c(0, 8),
     yaxs = "i")
axis(side = 1,
     at = seq(0L, 16L, by = 2L),
     labels = paste("$",
                    seq(108L, 124L, by = 2L),
                    "$"))
lines(x = utimes,
      y = cumhaz.log.confint[,1],
      lty = 2,
      type = "s")
lines(x = utimes,
      y = cumhaz.log.confint[,2],
      lty = 2,
      type = "s")
polygon(x = c(0,13,13,0),
        y = c(0,13/exp_fit[2], 13/exp_fit[3],0),
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
  system(command = paste0("lualatex ",fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}
