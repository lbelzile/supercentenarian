########################################################################
##   This script must be called from within Semi-supercentenarian.R   ##
########################################################################

#########################################################################################
# Fig 2 (bottom) Threshold stability plots for France 2019 data
#########################################################################################

dat <- francent$numdays - u
slow <- francent$slow
supp <- francent$supp


confint_xi <- matrix(0, ncol = 3, nrow = length(thresh))
confint_exp <- matrix(0, ncol = 3, nrow = length(thresh))
for(i in 1:length(thresh)){
  datu <- dat - thresh[i]
  ind <- which(datu > 0)
  confint_xi[i,] <- prof_gpd_dtrunc_xi(dat = (datu[ind])/365.25,
                                       supp = (supp[ind]-thresh[i])/365.25, 
                                       slow = (pmax(0, slow[ind]-thresh[i]))/365.25,
                                       confint = TRUE)
  confint_exp[i,] <- prof_exp_dtrunc(dat = (datu[ind])/365.25,
                                     supp = (supp[ind]-thresh[i])/365.25, 
                                     slow = (pmax(0, slow[ind]-thresh[i]))/365.25)
  
}

if(figures){
  fig <- "Fig2b.tex"
  setwd(fig_dir)
  tikz(fig, width = dwidth, height = 0.8*dheight, standAlone = TRUE)
}
# Consider two thresholds for which xi is negative
par(mar = c(4, 4, 0.1, 0.1), mfrow = c(1,2))
plot(x = thresh/365.25+105, confint_xi[,1], 
     type= "p", 
     panel.first = abline(h = 0, col = "grey"), 
     pch = 20, col = 1, bty = "l", 
     ylab = "$\\gamma$", 
     xlab = "threshold (years)",
     ylim = c(-0.3,1.1))
for(i in 1:length(thresh)){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_xi[i,2], y1 = confint_xi[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
plot(x = thresh/365.25+105, confint_exp[,1], type= "p", 
     panel.first = abline(h = 1.45142127, col = "grey"),
     pch = 20, col = 1,
     bty = "l", 
     ylab = "$\\sigma_e$", 
     xlab = "threshold (years)",
     ylim = c(0.7,2.1))
for(i in 1:length(thresh)){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_exp[i,2], y1 = confint_exp[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
if(figures){
  dev.off()
  system(command = paste0("pdflatex ",fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}