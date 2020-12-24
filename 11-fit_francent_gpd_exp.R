########################################################################
##   This script must be called from within Semi-supercentenarian.R   ##
########################################################################

#########################################################################################
# Analysis of "France 2019"
# French semi-supercentenarian from IDL
#########################################################################################

# Table 1, bottom panel

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


######################################################
# Aside: compare with incorrect parameter estimates  #
#  (obtained without accounting for the truncation)  #
######################################################
wparams <- t(sapply(thresh/365.25, function(i){ 
  mev::fit.gpd(dat/365.25, threshold = i)$est
  }))
wparams - param_gpd[,2:3]
# Systematic different in scale parameters, not really in shape
# Compare endpoints
ifelse(wparams[,2] >=0, Inf, -wparams[,1]/wparams[,2]) - ifelse(param_gpd[,'shape'] >= 0, Inf, -param_gpd[,'scale']/param_gpd[,'shape'])
# Compare mean survival time
meansurv1 <- thresh / 365.25 + 105 + param_gpd[,'scale']/(1 - param_gpd[,'shape'])
meansurv2 <- thresh / 365.25 + 105 + wparams[,1]/(1 - wparams[,2])
# Mean survival time - difference in days
round((meansurv2 - meansurv1)*365.25,0)
# Median survival time - difference in days
round((param_gpd[,'scale']*(2^param_gpd[,'shape']-1)/param_gpd[,'shape'] - wparams[,1]*(2^wparams[,2]-1)/wparams[,2])*365.25,0)


# Likelihood ratio test for differences between gender
# male versus female

# split data by gender
mint <- francent$gender == "male"
wint <- francent$gender == "female"

# fit truncated gpd to all data, then separately per gender
# the objective value is twice the log-likelihood
fjit <- fitgpddtrunc(dat = dat, thresh = thresh[1:11], supp = supp, slow = slow)
fmen <- fitgpddtrunc(dat = dat[mint], thresh = thresh[1:11], supp = supp[mint], slow = slow[mint])
fwom <- fitgpddtrunc(dat = dat[wint], thresh = thresh[1:11], supp = supp[wint], slow = slow[wint])

# Build likelihood ratio statistic
oyr <- seq(1,11, by = 2)
# Compute p-value using asymptotic null (2 df)
pvals <- pchisq((fwom[oyr,4] + fmen[oyr,4])-fjit[oyr,4], df = 2, lower.tail = FALSE)
pvalt <- rbind(paste0("$",round(100*pvals, digits= 2), "$"))
colnames(pvalt) <- paste0( "$", as.vector(round(fjit[oyr,1],1)), "$")
# What is the expect difference in survival of men versus women?
# Compute the expected survival time (mean of GPD is mu + sigma / (1-xi))
((fmen[7,1] + fmen[7,2]/(1-fmen[7,3])) - (fwom[7,1] + fwom[7,2]/(1-fwom[7,3])) )*365.25

# Likelihood ratio test for the exponential distribution
ejit <- fitgpddtrunc(dat = dat, thresh = thresh[1:11], supp = supp, slow = slow, expo = TRUE)
emen <- fitgpddtrunc(dat = dat[mint], thresh = thresh[1:11], supp = supp[mint], slow = slow[mint], expo = TRUE)
ewom <- fitgpddtrunc(dat = dat[wint], thresh = thresh[1:11], supp = supp[wint], slow = slow[wint], expo = TRUE)
# Compute p-value and then bind them
pvalt <- rbind(pvalt, paste0("$",round(100*pchisq((ewom[oyr,3] + emen[oyr,3])-ejit[oyr,3], df = 1, lower.tail = FALSE), 2),"$"))
# All differences are significative, except 110
rownames(pvalt) <- c("generalized Pareto", "exponential")
pvalt
