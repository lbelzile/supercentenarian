########################################################################
##   This script must be called from within Semi-supercentenarian.R   ##
########################################################################

#########################################################################################
# Figure 1: Lexis diagram of Istat data
#########################################################################################

# Split sample by sex and right-censoring status
menunc <- which(as.logical((as.numeric(italcent$gender)-1) * !italcent$rightcens))
mencens <- which(as.logical((as.numeric(italcent$gender)-1) * italcent$rightcens))
womunc <- which(as.logical(((as.numeric(italcent$gender)-1)+1) %%2 * !italcent$rightcens))
womcens <- which(as.logical(((as.numeric(italcent$gender)-1)+1) %%2 * italcent$rightcens))

if(figures){
  fig <- "Fig1.tex"
  setwd(fig_dir)
  tikz(fig, width = dwidth, height = dheight, standAlone = TRUE)
}
par(mar = c(4,4,1,3))
plot(x = (italcent$birth + italcent$numdays)[womunc], 
     y = italcent$numdays[womunc]/365.25,
     col = scales::alpha("black", alpha = 0.5),
     bty = "l", pch = 20, cex = 0.8, ylab = "age at death (in years)", xaxt = 'n',
     xlab = "year of death", yaxs = "i", ylim = c(105, 117), xlim = c(14244, 16825), xaxs = "i")
axis(side=1,at=seq(14244, 16825, by = 365.25), labels = 2009:2016, tick = TRUE)
points(x = (italcent$birth + italcent$numdays)[menunc], 
       y = italcent$numdays[menunc]/365.25,  pch = 4, cex = 0.8,
       col = scales::alpha("red", alpha = 0.5),lwd = 2)
rug(italcent$numdays[mencens]/365.25, side = 4, line =  1, col = scales::alpha("red", 0.5), ticksize = 0.05)
rug(italcent$numdays[womcens]/365.25, side = 4, line =  2, col = scales::alpha("black", 0.5), ticksize = 0.05)
# Number of death per year
deathcohorts_wom <- table(substr((italcent$birth + italcent$numdays)[womunc], 1, 4))
deathcohorts_men <- table(substr((italcent$birth + italcent$numdays)[menunc], 1, 4))
# Add counts to top of figure
text(labels = paste0(deathcohorts_wom, "/\\textcolor{red}{", deathcohorts_men,"}"), 
     x =as.Date(paste0(2009:2015, "-06-15")), y = 116, cex = 0.8)
if(figures){
  dev.off()
  system(command = paste0("pdflatex ", getwd(), "/", fig,"; rm *.aux; rm *.log"))
  setwd(code_dir)
}