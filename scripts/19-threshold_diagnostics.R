########################################################################
##   This script must be called from within Semi-supercentenarian.R   ##
########################################################################

# Custom tests for threshold selection
# Northrop and Coleman (2014) propose fitting
# a piecewise generalized Pareto with continuity
# constraints and performing a test of equality
# of shapes.
u105 <- 38350L
u108 <- 39447L
# Perform a likelihood ratio test
# comparing the piecewise generalized Pareto
# with the generalized Pareto distribution.


fit_gppiece <- with(francent,
     longevity::fit_elife(time = (numdays - u105)/365.25,
                          ltrunc = slow/365.25,
                          rtrunc = supp/365.25,
                          thresh = c(3, 5, 7), #108,110,112
                          family = "gppiece",
                          restart = TRUE)
     )

fit_gp <- with(francent,
               longevity::fit_elife(time = (numdays - u105)/365.25,
                                    ltrunc = slow/365.25,
                                    rtrunc = supp/365.25,
                                    thresh = 3,
                                    family = "gp",
                                    restart = TRUE)
)
fit_exp <- with(francent,
               longevity::fit_elife(time = (numdays - u105)/365.25,
                                    ltrunc = slow/365.25,
                                    rtrunc = supp/365.25,
                                    thresh = 3,
                                    family = "exp",
                                    restart = TRUE)
)
# Likelihood ratio test
anova(fit_gppiece, fit_gp)
anova(fit_gp, fit_exp)
anova(fit_gppiece, fit_exp)

## Threshold selection based on metric proposed in Varty et al.
## arXiv:2102.00884
## 
## Algorithm: simulate B bootstrap samples from model with MLE 
# 1) fit the parametric model
# 2) transform data to unit exponential
# 3) calculate sample quantile function and empirical distribution function
# 4) compute metrics
## Return the average metric for each threshold.

#' Transform interval truncated data to unit exponential scale
dtrunc_gp_to_exp <- function(time,
                             ltrunc,
                             rtrunc,
                             par){
  cdf <- function(x){
    longevity::pelife(x, 
                      scale = par[1],
                      shape = par[2],
                      family = "gp")
  }
  psunif <- (cdf(time) - cdf(ltrunc)) / (cdf(rtrunc) - cdf(ltrunc))
  return( -log(1 - psunif))
}

metric <- function(z, m = 1000L){
  pseq <- seq_len(m)/(m + 1)
  emp_cdf <- ecdf(z)(-log(1-pseq))
  quant_fn <- quantile(z, probs = pseq)
  # quantile, d = 1
  c(mean(abs(-log(1-pseq) - quant_fn)), #d(q,1)
    mean((-log(1-pseq) - quant_fn)^2), #d(q,2)
   mean(sqrt(sqrt(length(z))/(pseq*(1-pseq)))*abs(pseq - emp_cdf)), #d(p,1)
   mean(sqrt(sqrt(length(z))/(pseq*(1-pseq)))*(pseq - emp_cdf)^2) #d(p,2)
  )
}

th_seq <- seq(106, 112, length.out = 15L)
nt <- length(th_seq) # length of threshold
m <- 1000L # nb of points
B <- 1000L # nb of Monte Carlo replications
bootstrap <- array(NA, dim = c(B, nt, 4))
set.seed(1234)
for(i in seq_along(th_seq)){
  fit_gp <- with(francent,
                 longevity::fit_elife(time = (numdays - u105)/365.25,
                                      restart = TRUE,
                                      ltrunc = slow/365.25,
                                      rtrunc = supp/365.25,
                                      thresh = th_seq[i] - 105,
                                      family = "gp",
                                      export = TRUE))
for(b in seq_len(B)){
    boot_samp <- with(fit_gp, 
                      r_dtrunc_elife(n = nexc,
                                     scale = par[1],
                                     shape = par[2],
                                     lower = ltrunc,
                                     upper = rtrunc, 
                                     family = "gp"))
    boot_mle <- try(longevity::fit_elife(
      time = boot_samp,
      ltrunc = fit_gp$ltrunc,
      rtrunc = fit_gp$rtrunc,
      thresh = 0,
      family = "gp"))
    if(!is.character(boot_mle)){
    boot_samp_exp <- dtrunc_gp_to_exp(time = boot_samp,
                                      ltrunc = fit_gp$ltrunc,
                                      rtrunc = fit_gp$rtrunc,
                                      par = boot_mle$par)
    
  bootstrap[b,i,] <- metric(z = boot_samp_exp, m = m)
    }
  }  
}
#save(bootstrap, file = "bootstrap_varty.RData", version = 2)
method_means <- apply(bootstrap, 2:3, mean, na.rm = TRUE)
matplot(x = th_seq, 
        y = method_means, 
        type = "l", 
        col = 1,
        xlab = "threshold (years)", 
        ylab = "metric")
th_seq[apply(method_means, 2, which.min)]
# These results suggest picking the smallest possible threshold