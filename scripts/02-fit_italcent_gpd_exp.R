########################################################################
##   This script must be called from within Semi-supercentenarian.R   ##
########################################################################

#########################################################################
# Fit generalized Pareto and exponential distributions to 
# Istat data for multiple thresholds
########################################################################

# Exceedances
dat <- italcent$numdays - u
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
                rightcens = italcent$rightcens[ind], 
                slow = (pmax(0, italcent$slow[ind]-thresh[i]))/365.25, 
                expo = FALSE, 
                control = list(fnscale = c(1,0.1), reltol = 1e-12, maxit= 1e5),
                hessian = TRUE)
  # Compute standard errors
  endpt_stderror <- sqrt(solve(obs.infomat(c(- vals$par[1]/vals$par[2], vals$par[2]), 
                                           dat = (datu[ind])/365.25,
                                           rightcens = italcent$rightcens[ind], 
                                           slow = (pmax(0, italcent$slow[ind]-thresh[i]))/365.25))[1,1])
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
                rightcens =  italcent$rightcens[ind], slow = (pmax(0, italcent$slow[ind]-thresh[i]))/365.25, expo = TRUE, hessian = TRUE)
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

