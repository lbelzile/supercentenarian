#########################################################################################
# Extended generalized Pareto, Gompertz model (Tables 3 and 4)
#########################################################################################

# Turn this switch on to save the results of the simulation
save <- FALSE

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
if(tables){
gomp_par <- data.frame(thresh = paste0(105+thresh/365.25), 
                       param_gomp[,1:3])
table.res <- matrix(nrow = 4, ncol = 7)
ind <- 1:length(thresh)
table.res[1,] <- paste0("$", round(param_gpd[ind,"nu"],0),"$")
table.res[2,] <- paste0("$", round(param_gomp[ind, "scale"],2), ifelse(is.na(param_gomp[ind, "scale.stderror"]), paste0("\\; (",round(param_exp[ind, "scale.stderror"],1),")$"), paste0("\\; (",round(param_gomp[ind, "scale.stderror"],2),")$")))
table.res[3,] <- paste0("$", round(param_gomp[ind, "shape"],2), ifelse(is.na(param_gomp[ind, "shape.stderror"]), "$", paste0("\\; (",round(param_gomp[ind, "shape.stderror"],2),")$")))
# table.res[4,] <- paste0("$", round(param_exp[ind, "scale"],2), "\\; (",round(param_exp[ind, "scale.stderror"],1),")$")
# table.res[4,] <- paste0("$", round(0.5*(param_gomp[ind,"deviance"] <= param_exp[ind,"deviance"]) + pchisq(param_gomp[ind,"deviance"]- param_exp[ind,"deviance"],1, lower.tail = FALSE)/2, 2),"$")
table.res[4,] <- paste0("$", round(pvalboot,2), "$")
rownames(table.res) <- c("$n_u$", "$\\sigma$", "$\\beta$", "$p$-value")
xtab2 <- xtable::xtable(table.res, 
                        caption = "Estimates (standard errors) of Gompertz parameters ($\\beta$, $\\sigma$) for the Italian \\Istat{} data as a function of threshold, with number of threshold exceedances ($n_u$), $p$-value for the likelihood ratio test of $\\beta=0$. Estimates reported as zero for $\\beta$ are smaller than $10^{-7}$.",
                        align = paste0(rep("r",ncol(table.res)+1),collapse = ""))
addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <-  paste0("threshold", paste0(" & $", 105:111, "$", collapse = " "), "\\\\\n", collapse = " ")

print(xtab2, file = "tables/Table3a.tex", 
      size = "footnotesize", 
      floating.environment = "table*",  
      booktabs = TRUE, caption.placement = "top",
      add.to.row = addtorow, include.colnames = FALSE,
      sanitize.colnames.function = identity, 
      sanitize.text.function = identity, 
      sanitize.rownames.function = identity)
}
