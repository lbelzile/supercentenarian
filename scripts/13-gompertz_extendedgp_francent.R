########################################################################
##   This script must be called from within Semi-supercentenarian.R   ##
########################################################################

#########################################################################################
# Extended generalized Pareto, Gompertz model (Table 3, bottom panel)
#########################################################################################

# Turn this switch on/off to save/load the results of the simulation
save <- TRUE
dat <- francent$numdays - u
slow <- francent$slow
supp <- francent$supp
# Fit Gompertz model
param_gomp <- matrix(0, nrow = length(thresh), ncol = 7L)
colnames(param_gomp) <- c("thresh","scale","shape","deviance", "scale.stderror", "shape.stderror", "nu")
for(i in 1:length(thresh)){
  #For each threshold, compute new threshold exceedances
  datu <- dat - thresh[i]
  # Keep only exceedances
  ind <- which(datu > 0)
  datu <- datu[ind]/365.25
  supp <- (francent$supp[ind] - thresh[i])/365.25
  slow <- (pmax(0, francent$slow[ind] - thresh[i]))/365.25
  if(i == 1){
    start <- c(1.6,0.1)
  } else{
    start <- param_gomp[i-1,2:3]
  }
  # Compute maximum likelihood estimates numerically
  #   # Compute maximum likelihood estimates numerically
  test <- optim(par = start, 
                fn = function(par){
                  gomp_dtrunc(par = c(par,0), 
                              dat = datu, 
                              supp = supp, 
                              slow = slow)}, 
                method = "N")
  vals <- nlm(p = start, 
              f = function(par){
                gomp_dtrunc(par = c(par,0), 
                            dat = datu, 
                            supp = supp, 
                            slow = slow)},
              iterlim = 1e8,
              stepmax = 1, gradtol = 1e-8, hessian = TRUE, 
              typsize = start, 
              fscale = param_exp[i,"deviance"]/2)
  print(round((u + thresh[i])/365.25,2)); print(vals$gradient)
  # Bundle estimates in a vector
  std <- try(sqrt(diag(solve(numDeriv::hessian(func = function(par){
    gomp_dtrunc(par = par, dat = datu, 
                supp = supp, 
                slow = slow)}, x = vals$estimate)))))
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
param_gomp[param_gomp[,"shape"] < 1e-7, "scale.stderror"] <- param_exp[param_gomp[,"shape"] < 1e-7, "scale.stderror"]
# LRT test for exponential is non-regular - asymptotically half chi-square(1), half point mass at 0
# Compute the p-value using a parametric bootstrap.
if(save){
  B <- 1e4L
  set.seed(1245)
  boot_lrt_gomp_vs_exp <- matrix(0, nrow = B, ncol = length(thresh))
  for(j in 1:length(thresh)){
    #For each threshold, compute new threshold exceedances
    datu <- (francent$numdays - round(365.25*105) - thresh[j])/365.25
    ind <- which(datu > 0)
    datu <- datu[ind]/365.25
    supp <- (francent$supp[ind] - thresh[j])/365.25
    slow <- (pmax(0, francent$slow[ind] - thresh[j]))/365.25
    # Keep only exceedances
    bootsamp <- matrix(0, nrow = length(ind), ncol = B)
    for(i in 1:length(ind)){
      bootsamp[i,] <- qexp(pexp(slow[i], rate = 1/param_exp[j,"scale"]) + runif(B)*
                             (pexp(supp[i], rate = 1/param_exp[j,"scale"]) - pexp(slow[i], rate = 1/param_exp[j,"scale"])), 
                           rate = 1/param_exp[j,"scale"])
    }
    start <- param_gomp[j,2:3] + c(0,0.05)
    for(b in 2:B){
      # Fit the Gompertz model on the bth bootstrap sample
      test <- optim(par = start, 
                    fn = function(par){
                      gomp_dtrunc(par = par, dat = bootsamp[,b], 
                                  slow = slow, supp = supp)}, method = "N")
      vals <- nlm(p = test$par, 
                  f = function(par){
                    gomp_dtrunc(par = par, dat = bootsamp[,b], 
                                slow = slow, supp = supp)},
                  iterlim = 1e8,
                  stepmax = 1, gradtol = 1e-8, hessian = TRUE, 
                  typsize = start, 
                  fscale = param_exp[j,"deviance"]/2)
      mboot <- mean(bootsamp[,b], na.rm = TRUE)
      pexpboot <- nlm(f = function(x){
        gpd_dtrunc(par = x, dat = bootsamp[,b], supp = supp, slow = slow, expo = TRUE)
      }, p = mboot, typsize = mboot, gradtol = 1e-12)
      
      # try(optim(par = mboot, 
      #                    fn = gpd_dtrunc, 
      #                    method = "Brent", 
      #                    lower = 0.01*mboot, 
      #                    upper = 5*mboot, 
      #                    dat = bootsamp[,b],
      #                    supp = supp, 
      #                    slow = slow, 
      #                    expo = TRUE, 
      #                    hessian = FALSE))
      
      boot_lrt_gomp_vs_exp[b, j] <- 2*(pexpboot$minimum - vals$minimum)
    }
  }
  boot_lrt_gomp_vs_exp[1,] <- pmax(0, param_gomp[,"deviance"] - param_exp[,"deviance"])
  pvalboot <- apply(pmax(boot_lrt_gomp_vs_exp,0), 2, function(x){mean(x>=x[1])})
  save(file = "pvalboot_France_gompvsexp.Rdata", thresh, boot_lrt_gomp_vs_exp, pvalboot)
} else{
  load(file = "pvalboot_France_gompvsexp.Rdata")
}
if(tables){
  gomp_par <- data.frame(thresh = paste0(105+thresh/365.25), 
                         param_gomp[,1:3])
  table.res <- matrix(nrow = 5, ncol = 7)
  ind <- seq(1L, 13L, by = 2L)
  table.res[1,] <- paste0("$", round(param_gpd[ind,"nu"],0),"$")
  table.res[2,] <- paste0("$", round(param_gomp[ind, "scale"],2), ifelse(is.na(param_gomp[ind, "scale.stderror"]), paste0("\\; (",round(param_exp[ind, "scale.stderror"],1),")$"), paste0("\\; (",round(param_gomp[ind, "scale.stderror"],2),")$")))
  table.res[3,] <- paste0("$", round(param_gomp[ind, "shape"],2), ifelse(is.na(param_gomp[ind, "shape.stderror"]), "$", paste0("\\; (",round(param_gomp[ind, "shape.stderror"],2),")$")))
  table.res[5,] <- paste0("$", round(pvalboot[ind],2),"$")
  table.res[4,] <- paste0("$", round(I((param_exp[ind,"deviance"] + param_gomp[ind,"deviance"])< 0)*0.5 + pchisq(param_exp[ind,"deviance"] + param_gomp[ind,"deviance"],1, lower.tail = FALSE)/2, 2),"$")
  rownames(table.res) <- c("$n_u$", "$\\sigma$", "$\\beta$", "$p$-value (asymptotic)", "$p$-value (bootstrap)")
  xtab2 <- xtable::xtable(table.res, 
                          caption = "Estimates (standard errors) of Gompertz parameters ($\\beta$, $\\sigma$) for the Italian \\Istat{} data as a function of threshold, with number of threshold exceedances ($n_u$), $p$-value for the likelihood ratio test of $\\beta=0$. Estimates reported as zero for $\\beta$ are smaller than $10^{-7}$.",
                          align = paste0(rep("r",ncol(table.res)+1),collapse = ""))
  setwd(table_dir)
  addtorow <- list()
  addtorow$pos <- list(0)
  addtorow$command <-  paste0("threshold", paste0(" & $", 105:111, "$", collapse = " "), "\\\\\n", collapse = " ")
  print(xtab2, file = "Table3b.tex", 
        size = "footnotesize", 
        floating.environment = "table*",  
        booktabs = TRUE, caption.placement = "top",
        add.to.row = addtorow, include.colnames = FALSE,
        sanitize.colnames.function = identity, 
        sanitize.text.function = identity, 
        sanitize.rownames.function = identity)
  setwd(code_dir)
}