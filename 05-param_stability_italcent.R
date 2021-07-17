########################################################################
##   This script must be called from within Semi-supercentenarian.R   ##
########################################################################

#########################################################################################
# Figure 2a: parameter stability plots for Istat data
#########################################################################################
# Profile likelihood for shape of GP and scale of exponential distribution
# used to compute the parameter stability plots
prof_xizero <- rep(0, length.out = length(thresh))
confint_xi <- matrix(0, ncol = 3, nrow = length(thresh))
confint_exp <- matrix(0, ncol = 3, nrow = length(thresh))
for(i in 1:length(thresh)){
  datu <- dat - thresh[i]
  ind <- which(datu > 0)
  confint_xi[i,] <- prof_gpd_cens_xi_confint(dat = (datu[ind])/365.25,
                                             rightcens = rightcens[ind], slow = (pmax(0, italcent$slow[ind]-thresh[i]))/365.25)
  confint_exp[i,] <- prof_exp_cens(dat = (datu[ind])/365.25,
                                   rightcens = rightcens[ind], slow = (pmax(0, italcent$slow[ind]-thresh[i]))/365.25)
  prof_xizero[i] <- nllhxizero(dat = (datu[ind])/365.25, rightcens = rightcens[ind], slow = (pmax(0, italcent$slow[ind]-thresh[i]))/365.25)
  
}

# Compute the probability that xi >= 0 based on profile and likelihood root
probgamma0 <- pnorm(sign(param_gpd[,'shape'])*sqrt(param_gpd[,'deviance'] - prof_xizero))

if(figures){
  fig <- "Fig2a.tex"
  setwd(fig_dir)
  tikz(fig, width = dwidth, height = 0.8*dheight, standAlone = TRUE)
}
par(mar = c(4, 4, 0.1, 0.1), mfrow = c(1,2))
plot(x = thresh/365.25+105, confint_xi[,1], type= "p", panel.first = abline(h = 0, col = "grey"), 
     pch = 20, col = 1, bty = "l", ylab = "$\\gamma$", xlab = "threshold (in years)",ylim = c(-0.3,1.1))
for(i in 1:length(thresh)){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_xi[i,2], y1 = confint_xi[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
plot(x = thresh/365.25+105, confint_exp[,1], type= "p", 
     panel.first = abline(h=1.45142127, col = "grey"), pch = 20, col = 1,
     bty = "l", ylab = "$\\sigma_e$", xlab = "threshold (in years)",
     ylim = c(0.7,2.1))
for(i in 1:length(thresh)){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_exp[i,2], y1 = confint_exp[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}

if(figures){
  dev.off()
  system(command = paste0("pdflatex ", fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}
