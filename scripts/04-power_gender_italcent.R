########################################################################
##   This script must be called from within Semi-supercentenarian.R   ##
########################################################################

########################################################################
# Power study - comparison of men versus women (SM, Section E)         #
########################################################################
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

## The following figure is not included in the paper

par(mar = c(4,4,1,1), pch = 20, bty = "l")
plot(colMeans(lrt_stata > qchisq(0.95,1)) ~ lam_vals,
     ylab = "power", xlab = expression(lambda))
## Save the results to avoid repeating the analysis
# save(lam_vals, lrt_stat,lrt_statb, file = "power.RData")
 
# To find the lambda such that the power is X
# We fit a spline but flip x and y
# and predict at X = {0.2, 0.5, 0.8}

# First, remove points at which the power is too high/too low
x <- colMeans(lrt_stata > qchisq(0.95,1))
flat <- which(x > 0.985, x < 0.015)
x <- x[-flat]
y <- lam_vals[-flat]

# Fit a penalized spline
m <- mgcv::gam(y ~ s(x))
plot(y~x)
lines(fitted(m) ~ x)
pow_pred <- mgcv::predict.gam(m, newdata = data.frame(x = c(0.2,0.5,0.8)))
# Predicted lambda
pow_pred
# Difference in expected survival for men < women (in days)
# exponential distribution = mean, rounded to the nearest day
round((pow_pred/r108 - 1/(r108*pow_pred)),0)

## Additional plot not included in the paper
# Create plot as a function of lambda of the power
par(mar = c(5.5,4,1,1), bty = "l")
# Create empty frame
plot(NULL, NULL, ylab = "power", 
     xlab = "expected difference in survival time\n women versus men (in years)",
     yaxs = "i",
     xaxs = "i",
     ylim = c(0,1),
     xlim = c(0, 1.2),
     panel.first = {abline(h = seq(0.1, 1, by = 0.1), 
                           col = scales::alpha("black", 0.05), lty = 1, lwd = 0.25)})
gam1 <- mgcv::gam(colMeans(lrt_stata > qchisq(0.95,1)) ~ s(lam_vals^2, k = -1, bs = "tp"))
lines((lam_vals/r108 - 1/(r108*lam_vals))/365.25, fitted(gam1), lwd = 2)
gam2 <- mgcv::gam(colMeans(lrt_statb > qchisq(0.95,1)) ~ s(lam_vals^2, k = -1, bs = "tp"))
lines((lam_vals/r108 - 1/(r108*lam_vals))/365.25, fitted(gam2), lwd = 2, lty = 3,  col = "grey30")
legend("bottomright", legend = c("censoring","no censoring"), bty = "n",
       col = c(1, "grey30"), lwd = 2, lty = c(1,3))
dev.off()