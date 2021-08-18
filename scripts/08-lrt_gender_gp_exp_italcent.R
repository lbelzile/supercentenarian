########################################################################
##   This script must be called from within Semi-supercentenarian.R   ##
########################################################################

#########################################################################################
# Testing difference between sex - perform likelihood ratio tests (SI, Section D)
#########################################################################################

# Generalized Pareto model - men only
param_gpd_men <- matrix(0, nrow = 8, ncol = 9L)
colnames(param_gpd_men) <- c("loc","scale","shape","deviance", "scale.stderror", "shape.stderror","endpoint", "endpoint.stderror", "nu")
for(i in 1:8){ # not enough males (12) beyond this stage to estimate parameters of the GP
  datu <- dat - thresh[i]
  ind <- as.logical(I(datu > 0) * I(italcent$gender=="male"))
  ind <- which(ind)
  vals <- optim(par = c(1.2,-0.1), fn = gpd_cens, method = "N", dat = (datu[ind])/365.25,
                rightcens = rightcens[ind], slow = (pmax(0, itcalent$slow[ind]-thresh[i]))/365.25, expo = FALSE, hessian = TRUE)
  endpt_stderror <- sqrt(solve(obs.infomat(c(- vals$par[1]/vals$par[2], vals$par[2]), dat = (datu[ind])/365.25,rightcens = rightcens[ind], slow = (pmax(0, italcent$slow[ind]-thresh[i]))/365.25))[1,1])
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
                rightcens = rightcens[ind], slow = (pmax(0, italcent$slow[ind]-thresh[i]))/365.25, expo = FALSE, hessian = TRUE)
  endpt_stderror <- sqrt(solve(obs.infomat(c(- vals$par[1]/vals$par[2], vals$par[2]), dat = (datu[ind])/365.25,rightcens = rightcens[ind], slow = (pmax(0, italcent$slow[ind]-thresh[i]))/365.25))[1,1])
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
                rightcens = rightcens[ind], slow = (pmax(0, italcent$slow[ind]-thresh[i]))/365.25, expo = TRUE, hessian = TRUE)
  param_exp_women[i,] <- c(u + thresh[i], vals$par, -2*vals$value, sqrt(solve(vals$hessian)))
}

param_exp_men <- matrix(0, nrow = 9, ncol = 4L)
colnames(param_exp_men) <- c("loc","scale","deviance", "scale.stderror")
for(i in 1:9){
  datu <- dat - thresh[i]
  ind <- as.logical(I(datu > 0) * I(italcent$gender=="male"))
  vals <- optim(par = c(500), fn = gpd_cens, method = "Brent", lower = 0, upper = 800, dat = (datu[ind])/365.25,
                rightcens = rightcens[ind], slow = (pmax(0, italcent$slow[ind]-thresh[i]))/365.25, expo = TRUE, hessian = TRUE)
  param_exp_men[i,] <- c(u + thresh[i], vals$par, -2*vals$value, sqrt(solve(vals$hessian)))
}
# P-values
1-pchisq((param_exp_women[1:9,"deviance"] + param_exp_men[1:9,"deviance"]) - param_exp[1:9,"deviance"], df=1)
