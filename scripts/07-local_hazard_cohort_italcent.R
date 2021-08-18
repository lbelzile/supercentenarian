########################################################################
##   This script must be called from within Semi-supercentenarian.R   ##
########################################################################

#########################################################################################
# Fig. 4 Local hazard plots by birth cohort
#########################################################################################

# Local hazard using all data (not shown)
localhazard <- prof_gpd_hazard_confint(dat = dat, rightcens = rightcens, slow = italcent$slow, thresh = param_gpd[1,1]+0:6)

ind_l1 <-  which(italcent$birth < as_date("1906-01-01"))
localhazard1 <- prof_gpd_hazard_confint(dat = dat[ind_l1], rightcens = rightcens[ind_l1], slow = italcent$slow[ind_l1], thresh = 105:111)
ind_l2 <-  which(italcent$birth >= as_date("1906-01-01"))
localhazard2 <- prof_gpd_hazard_confint(dat = dat[ind_l2], rightcens = rightcens[ind_l2], slow = italcent$slow[ind_l2], thresh = 105:109)
# Note: some optimization calls yield warnings, but these are unconsequential

if(figures){
  fig <- "Fig4.tex"
  setwd(fig_dir)
  tikz(fig, width = dwidth, height = dheight, standAlone = TRUE)
}
# Consider two thresholds for which xi is negative
yran <- range(c(0, c(localhazard1[,1:3], localhazard2[,1:3])))
par(mar = c(4.5, 4.5, 0.1, 0.1), mfrow = c(1,2))
plot(105:111, localhazard1[,1],  ylim = yran, panel.first = abline(h = 0.7, col = "grey"), 
     type= "p", pch = 19, bty = 'l', ylab = "local hazard", xlab = "age (in years)")
for(i in 1:nrow(localhazard1)){
  arrows(x0 = (105:111)[i], y0 = localhazard1[i,2], y1 = localhazard1[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}
plot(105:109, localhazard2[,1],  ylim = yran, 
     type= "p", pch = 19, bty = 'l', ylab = "local hazard", 
     panel.first = abline(h = 0.7, col = "grey"), xlab = "age (in years)")
for(i in 1:nrow(localhazard2)){
  arrows(x0 = (105:109)[i], y0 = localhazard2[i,2], y1 = localhazard2[i,3], 
         length = 0.02, angle = 90, code = 3, lwd = 2)
}

if(figures){
  dev.off()
  system(command = paste0("pdflatex ", fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}
