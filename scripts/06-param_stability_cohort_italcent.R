########################################################################
##   This script must be called from within Semi-supercentenarian.R   ##
########################################################################

#########################################################################################
# Fig. 5 Parameter stability plot, by cohort 1896-1905 vs 1906-1910
#########################################################################################

confint_xi1 <- matrix(0, ncol = 3, nrow = length(thresh))
confint_exp1 <- matrix(0, ncol = 3, nrow = length(thresh))
confint_xi2 <- matrix(0, ncol = 3, nrow = 8)
confint_exp2 <- matrix(0, ncol = 3, nrow = 8)
for(i in 1:length(thresh)){
  datu <- dat - thresh[i]
  indi <- intersect(which(datu > 0), which(italcent$birth < as_date("1906-01-01"))) 
  confint_xi1[i,] <- prof_gpd_cens_xi_confint(dat = (datu[indi])/365.25,
                                              rightcens = rightcens[indi], slow = (pmax(0, italcent$slow[indi]-thresh[i]))/365.25)
  confint_exp1[i,] <- prof_exp_cens(dat = (datu[indi])/365.25,
                                    rightcens = rightcens[indi], slow = (pmax(0, italcent$slow[indi]-thresh[i]))/365.25)
  indi <- intersect(which(datu > 0), which(italcent$birth >= as_date("1906-01-01"))) 
  if(length(indi) > 20 && i < 9){
    confint_xi2[i,] <- prof_gpd_cens_xi_confint(dat = (datu[indi])/365.25,
                                                rightcens = rightcens[indi], slow = (pmax(0, italcent$slow[indi]-thresh[i]))/365.25)
    confint_exp2[i,] <- prof_exp_cens(dat = (datu[indi])/365.25,
                                      rightcens = rightcens[indi], slow = (pmax(0, italcent$slow[indi]-thresh[i]))/365.25)
  }
}

if(figures){
  fig <- "Fig5.tex"
  setwd(fig_dir)
  tikz(fig, width = dwidth, height = 1.6*dheight, standAlone = TRUE)
}
par(mar = c(4.5, 4.5, 0.1, 0.1), mfrow = c(2,2))
plot(x = thresh/365.25+105, confint_xi1[,1], type= "p", panel.first = abline(h = 0, col = "grey"), 
     pch = 19, col = 1, bty = "l", ylab = "$\\gamma$", xlab = "threshold (in years)",ylim =  c(-0.5,1))
for(i in 1:length(thresh)){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_xi1[i,2], y1 = confint_xi1[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
plot(x = thresh/365.25+105, confint_exp1[,1], type= "p", 
     panel.first = abline(h=1.45142127, col = "grey"), pch = 19, col = 1,
     bty = "l", ylab = "$\\sigma_e$", xlab = "threshold (in years)",
     ylim = c(0.7,2.1))
for(i in 1:length(thresh)){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_exp1[i,2], y1 = confint_exp1[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}

plot(x = thresh[1:8]/365.25+105, confint_xi2[,1], type= "p", 
     panel.first = abline(h = 0, col = "grey"), 
     pch = 19, col = 1, bty = "l", ylab = "$\\gamma$", xlab = "threshold (in years)",ylim = c(-0.5,1))
for(i in 1:8){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_xi2[i,2], y1 = confint_xi2[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
plot(x = thresh[1:8]/365.25+105, confint_exp2[1:8,1], type= "p", 
     panel.first = abline(h=1.45142127, col = "grey"), pch = 19, col = 1,
     bty = "l", ylab = "$\\sigma_e$", xlab = "threshold (in years)",
     ylim = c(0.7, 2.1))
for(i in 1:8){
  arrows(x0 = thresh[i]/365.25+105, y0 = confint_exp2[i,2], y1 = confint_exp2[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
dev.off()
system(command = paste0("pdflatex ", fig_dir, "/", fig,"; rm *.aux; rm *.log"))
setwd(code_dir)

