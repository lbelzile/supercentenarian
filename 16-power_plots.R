library(cobs)
library(ggplot2)
library(viridis)
library(patchwork)

load("italcent.rda")
maxital <- max(italcent$numdays/365.25)
load("francent.rda")
maxfran <- max(francent$numdays/365.25)
load("IDL2016.rda")
idl2016 <- idl2016[!idl2016$countrydeath == "FRA",]
maxidl <- max(idl2016$numdays/365.25)

# Results from 14-power_shape.R output
load("power.RData") 
nxis <- length(xis)
power_italcent[,xis == 0,]

powerIstat <- colMeans(cbind(abs(power_italcent[,,2]) > qnorm(0.975),
                             abs(power_italcent[,,1]) > qnorm(0.975),
                             power_italcent[,,2] < qnorm(0.05),
                             power_italcent[,,1] < qnorm(0.05)), na.rm = TRUE)
powerfrancent <- colMeans(cbind(abs(power_francent[,,2]) > qnorm(0.975),
                                abs(power_francent[,,1]) > qnorm(0.975),
                                power_francent[,,2] < qnorm(0.05),
                                power_francent[,,1] < qnorm(0.05)), na.rm = TRUE)
powerIDL <- colMeans(cbind(abs(power_idl[,,2]) > qnorm(0.975),
                           abs(power_idl[,,1]) > qnorm(0.975),
                           power_idl[,,2] < qnorm(0.05),
                           power_idl[,,1] < qnorm(0.05)), na.rm = TRUE)
powerCombo <- colMeans(cbind(abs(power_combo[,,2]) > qnorm(0.975),
                             abs(power_combo[,,1]) > qnorm(0.975),
                             power_combo[,,2] < qnorm(0.05),
                             power_combo[,,1] < qnorm(0.05)), na.rm = TRUE)

power_df <- data.frame(
  shape = rep(xis, length.out = 16*nxis),
  test = factor(rep(rep(c("Wald","directed likelihood root"), each = nxis), length.out = 16*nxis)),
  data = factor(rep(c("Istat","France","IDL2016","combined"), each = 4*nxis)),
  hypothesis = factor(rep(rep(c("two-sided","one-sided"), each = 2*nxis), length.out = 16*nxis)),
  power = c(powerIstat, powerfrancent, powerIDL, powerCombo))


# Results from 15-power_endpoint.R
load("power_endpoint.RData") 

powerIDL_ep <- colMeans(power_ep_idl < qnorm(0.05), na.rm = TRUE)
powerIstat_ep <- colMeans(power_ep_italcent < qnorm(0.05), na.rm = TRUE)
powerFrance_ep <- colMeans(power_ep_francent < qnorm(0.05), na.rm = TRUE)
powerCombo_ep <- colMeans(power_ep_combo < qnorm(0.05), na.rm = TRUE)


# Smoothing of the power curve to get rid of the artifacts of Monte-Carlo variability

nknots <- 40L
endpts <- seq(115, max(endpoints), by = 0.1)
p1IT <- is.na(powerIDL_ep)
powersmoothIstat <- predict(cobs::cobs(
  x = endpoints[-p1IT], 
  y = powerIstat_ep[-p1IT],
  pointwise = cbind(0, maxital, 1),   
  constraint = "decrease", 
  nknots = nknots), z = endpts)[,"fit"]
powersmoothIstat[endpts < maxital] <- 1
plot(endpts, powersmoothIstat, type = "l")
points(endpoints, powerIstat_ep)

p1FR <- is.na(powerFrance_ep)
powersmoothFrance <- predict(
  cobs::cobs(x = endpoints[-p1FR], 
             y = powerFrance_ep[-p1FR], 
             pointwise = cbind(0, maxfran, 1),
             constraint = "convex",
             nknots = 25),
  z = endpts)[,"fit"]
powersmoothFrance[endpts < maxfran] <- 1
plot(endpts, powersmoothFrance, type = "l")
points(endpoints, powerFrance_ep)

p1IDL <- is.na(powerIDL_ep)
powersmoothIDL <- predict(cobs::cobs(
  x = endpoints[-p1IDL], 
  y = powerIDL_ep[-p1IDL], 
  pointwise = cbind(0, maxidl, 1),
  constraint = "decrease", nknots = nknots),
  z = endpts)[,"fit"]
powersmoothIDL[endpts < maxidl] <- 1
plot(endpts, powersmoothIDL, type = "l")
points(endpoints, powerIDL_ep)


powersmoothCombo <- predict(
  cobs::cobs(x = endpoints[-p1FR], 
             y = powerCombo_ep[-p1FR], 
             pointwise = cbind(0, maxfran, 1),
             constraint = "decrease", 
             nknots = nknots),
  z = endpts)[,"fit"]
powersmoothCombo[endpts < maxfran] <- 1
plot(endpts, powersmoothCombo, type = "l")
points(endpoints, powerCombo_ep)


power_endp <- data.frame(
  endpoint = rep(endpts, length.out = 4*length(endpts)),
  data = factor(rep(c("Istat","France","IDL2016","combined"), each = length(endpts))),
  power = c(powersmoothIstat, powersmoothFrance, powersmoothIDL, powersmoothCombo)
)

power_df_ep <- data.frame(
  endpoint = rep(endpoints, length.out = 4*length(endpoints)),
  data = factor(rep(c("Istat","France","IDL2016","combined"), each = length(endpoints))),
  power = c(powerIstat_ep, powerFrance_ep, powerIDL_ep, powerCombo_ep)
)





lifetimes <- data.frame(
  death = c(italcent$numdays[italcent$numdays>365.25*115],
            francent$numdays[francent$numdays>365.25*115],
            idl2016$numdays[idl2016$numdays>365.25*115]
  )/365.25,
  data = factor(c(rep("Istat", sum(italcent$numdays>365.25*115)),
                  rep("France", sum(francent$numdays>365.25*115)),
                  rep("IDL2016", sum(idl2016$numdays>365.25*115)))))


critIstat <- apply(power_italcent[,xis == 0,], 2, quantile, c(0.025,0.05, 0.975))
# powerIstat <- colMeans(cbind(abs(power_italcent[,,2]) > qnorm(0.975),
#                              abs(power_italcent[,,1]) > qnorm(0.975),
#                              power_italcent[,,2] < qnorm(0.05),
#                              power_italcent[,,1] < qnorm(0.05)), na.rm = TRUE)
powerIstat <- colMeans(cbind((power_italcent[,,2] > critIstat[3,2])&(power_italcent[,,2] < critIstat[1,2]),
                             (power_italcent[,,1] > critIstat[3,1])&(power_italcent[,,1] < critIstat[1,1]),
                             power_italcent[,,2] < critIstat[2,2],
                             power_italcent[,,1] < critIstat[2,1]), na.rm = TRUE)
apply(power_italcent, c(2,3), function(x){sum(is.na(x))})
critFran <- apply(power_francent[,xis == 0,], 2, quantile, c(0.025,0.05, 0.975))
# powerfrancent <- colMeans(cbind(abs(power_francent[,,2]) > qnorm(0.975),
#                                 abs(power_francent[,,1]) > qnorm(0.975),
#                                 power_francent[,,2] < qnorm(0.05),
#                                 power_francent[,,1] < qnorm(0.05)), na.rm = TRUE)
powerfrancent <- colMeans(cbind((power_francent[,,2] > critFran[3,2])&(power_francent[,,2] < critFran[1,2]),
                                (power_francent[,,1] > critFran[3,1])&(power_francent[,,1] < critFran[1,1]),
                                power_francent[,,2] < critFran[2,2],
                                power_francent[,,1] < critFran[2,1]), na.rm = TRUE)
apply(power_francent, c(2,3), function(x){sum(is.na(x))})
critIDL <- apply(power_idl[,xis == 0,], 2, quantile, c(0.025,0.05, 0.975))
# powerIDL <- colMeans(cbind(abs(power_idl[,,2]) > qnorm(0.975),
#                            abs(power_idl[,,1]) > qnorm(0.975),
#                            power_idl[,,2] < qnorm(0.05),
#                            power_idl[,,1] < qnorm(0.05)), na.rm = TRUE)
apply(power_idl, c(2,3), function(x){sum(is.na(x))})
powerIDL <- colMeans(cbind((power_idl[,,2] > critIDL[3,2])&(power_idl[,,2] < critIDL[1,2]),
                           (power_idl[,,1] > critIDL[3,1])&(power_idl[,,1] < critIDL[1,1]),
                           power_idl[,,2] < critIDL[2,2],
                           power_idl[,,1] < critIDL[2,1]), na.rm = TRUE)
apply(power_combo, c(2,3), function(x){sum(is.na(x))})
critCombo <- apply(power_combo[,xis == 0,], 2, quantile, c(0.025,0.05, 0.975), na.rm = TRUE)
# powerCombo <- colMeans(cbind(abs(power_combo[,,2]) > qnorm(0.975),
#                              abs(power_combo[,,1]) > qnorm(0.975),
#                              power_combo[,,2] < qnorm(0.05),
#                              power_combo[,,1] < qnorm(0.05)), na.rm = TRUE)
powerCombo <- colMeans(cbind((power_combo[,,2] > critCombo[3,2])&(power_combo[,,2] < critCombo[1,2]),
                             (power_combo[,,1] > critCombo[3,1])&(power_combo[,,1] < critCombo[1,1]),
                             power_combo[,,2] < critCombo[2,2],
                             power_combo[,,1] < critCombo[2,1]), na.rm = TRUE)
# For the shape parameter, forcing all three shape parameters to be zero

powerAny<- 1-(1-powerIDL)*(1-powerfrancent)*(1-powerIstat)

power_df <- data.frame(
  shape = rep(xis, length.out = 16*nxis),
  test = factor(rep(rep(c("Wald","directed likelihood root"), each = nxis), length.out = 16*nxis)),
  data = factor(rep(c("Istat","France","IDL2016","combined"), each = 4*nxis)),
  hypothesis = factor(rep(rep(c("two-sided","one-sided"), each = 2*nxis), length.out = 16*nxis)),
  power = c(powerIstat, powerfrancent, powerIDL, powerCombo))

power_any <- data.frame(shape = rep(xis, length.out = 4*nxis),
                        test = factor(rep(rep(c("Wald","directed likelihood root"), each = nxis), length.out = 2*nxis)),
                        hypothesis = factor(rep(c("two-sided","one-sided"), each = 2*nxis)),
                        power = powerAny)
power_anyW <- power_any %>% filter((hypothesis == "one-sided")&(test == "Wald"))
# power_anyE <- data.frame(endpoint = endpts, 
#                          power = 1-(1-powersmoothIstat)*(1-powersmoothFrance)*(1-powersmoothIDL))

g1 <- power_df %>% filter((hypothesis == "one-sided")&(test == "directed likelihood root")) %>%
  ggplot(aes(x = shape, y = power, col = data)) + 
  geom_hline(yintercept = 0.05, alpha = 0.5) + 
  geom_smooth(method = "gam", se = FALSE,
              formula = y ~ s(x, bs = "tp", fx = TRUE, k=25), show.legend = FALSE) +
  theme_classic() + 
  theme(panel.grid.major = element_line(),
        text = element_text(size = 16))  + #seq(0,1, by = 0.1)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.01),
                     breaks = seq(0,1, by= 0.25),
                     labels = c("$0$","$0.25$","$0.5$", "$0.75$", "$1$"))  + 
  scale_x_continuous(breaks = seq(-0.25,0, by = 0.05),
                     limits = c(-0.25,0),
                     labels = paste0("$",seq(-0.25,0, by = 0.05),"$")) + 
  xlab("$\\gamma$")
Vcols <- viridis(4)[c(2:4,1)]
g1b <- power_df %>% filter((hypothesis == "one-sided")&(test == "Wald")&(data != "combined")) %>%
  ggplot(aes(x = shape, y = power)) + 
  geom_hline(yintercept = 0.05, alpha = 0.5) + 
  geom_smooth(mapping = aes(col = data), method = "gam", se = FALSE,
              formula = y ~ s(x, bs = "tp", fx = TRUE, k=25), show.legend = FALSE) +
  geom_smooth(data = power_anyW, aes(x = shape, y = power), color = "black", lty = 2, method = "gam", se = FALSE,
              formula = y ~ s(x, bs = "tp", fx = TRUE, k=25), show.legend = FALSE) +
  scale_color_manual(values = Vcols) + 
  theme_classic() + 
  theme(panel.grid.major = element_line(), 
        legend.position = "bottom",
        text = element_text(size = 16))  + #seq(0,1, by = 0.1)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.01),
                     breaks = seq(0,1, by= 0.25),
                     labels = c("$0$","$0.25$","$0.5$", "$0.75$", "$1$"))  + 
  scale_x_continuous(breaks = seq(-0.25,0, by = 0.05),
                     limits = c(-0.25,0),
                     labels = paste0("$",seq(-0.25,0, by = 0.05),"$")) + 
  xlab("$\\gamma$")

g2 <- ggplot(data = power_endp, aes(x = endpoint, y = power, col = data)) + 
  geom_hline(yintercept = 0.05, alpha = 0.5) + 
  geom_rug(data=lifetimes, mapping = aes(x = death, y = NULL, col = data)) + 
  geom_line(data = power_endp, 
            aes(x = endpoint, 
                y = power, col = data),
            size = 1) + 
  scale_color_viridis_d() + 
  theme_classic() + 
  theme(legend.position = "bottom",
        panel.grid.major = element_line(),
        text = element_text(size = 16)) + 
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1.01),
                     breaks = seq(0,1, by= 0.25),
                     labels = c("$0$","$0.25$","$0.5$", "$0.75$", "$1$"))  + 
  scale_x_continuous(breaks = seq(120, 150, by = 10),
                     limits = c(min(endpoints),151),
                     expand = c(0,0),
                     labels = paste0("$",seq(120, 150, by = 10),"$")) + 
  xlab("human lifespan limit") 


if(figures){
  fig <- "Fig3.tex"
  setwd(fig_dir)
  tikz(fig, width = 8, height = 4, standAlone = TRUE)
}
g2 + g1b + plot_layout(guides = 'collect') & theme(legend.position="bottom")
if(figures){
  dev.off()
  system(command = paste0("pdflatex ",fig_dir, "/", fig,"; rm *.aux; rm *.log"))
  setwd("..")
}

# Compute the power at different values of the endpoint
# to add in the manuscript text
round(cbind(endpoints, 
        powerIstat_ep,
        powerIDL_ep,
        powerFrance_ep, 
        powerCombo_ep)[c(45,50,55),]
,2)
