
# Table 2: parameter estimates for exponential
# distribution above 108 (Istat, France) and 110 (IDL2016)
# with gender-specific estimates and sample sizes

u108 <- 39447L #108 years
load("italcent.rda")
# Calendar date at which individual reaches 108 years
italcent$xcal <- italcent$birth + u108
# Calendar time for sampling frame
c1 <- lubridate::dmy("01-01-2009")
c2 <- lubridate::dmy("01-01-2016")

# Lower truncation level, zero if individual reached 108 between c1 and c2
italcent$slow <- as.numeric(pmax(0, c1 - italcent$xcal))
# Exceedances
datu <- italcent$numdays - u108
# Fit a generalized Pareto model
ind <- which(datu > 0)
uplim <- as.numeric(c2 - italcent$birth)
# Keep only exceedances above 108
uplim <- (uplim[ind] - u108)/365.25
rcens <- italcent$rightcens[ind]
gender <- italcent$gender[ind]
datu <- datu[ind]/365.25
slow <- italcent$slow[ind]/365.25
mle_istat_exp <- c(
  exp_mle_lt_rc(dat = datu, 
                slow = slow, 
                rightcens = rcens),
  exp_mle_lt_rc(dat = datu[gender == "female"],
                slow = slow[gender == "female"],
                rightcens = rcens[gender == "female"]),
  exp_mle_lt_rc(dat = datu[gender == "male"],
                slow = slow[gender == "male"],
                rightcens = rcens[gender == "male"])
)
nsamp_istat <- c(sum(!rcens), 
                 sum(!rcens[gender == "female"]),
                 sum(!rcens[gender == "male"]))
stdev_istat_exp <- mle_istat_exp/sqrt(nsamp_istat)
lci_istat_exp <- mle_istat_exp - qnorm(0.975)*stdev_istat_exp
uci_istat_exp <- mle_istat_exp + qnorm(0.975)*stdev_istat_exp


load("francent.rda")
# Calendar date at which individual reaches 105 years
xcal <- francent$birth + u108
# Calendar time for sampling frame
c1a <- lubridate::dmy("01-01-1978")
c1b <- lubridate::dmy("01-01-1987")
c2 <- lubridate::dmy("31-12-2017")

# Lower truncation level, zero if individual reached 105 between c1 and c2
francent$slow <- ifelse(francent$numdays > 110*365.24, 
                        as.numeric(pmax(0, c1b - xcal)),
                        as.numeric(pmax(0, c1a - xcal)))
francent$supp <- as.numeric(pmax(0, c2 - xcal))    
datu <- francent$numdays - u108
ind <- which(datu > 0)
slow <- francent$slow[ind]/365.25
supp <- francent$supp[ind]/365.25
datu <- (francent$numdays - u108)[ind]/365.25
gender <- francent$gender[ind]
opt_francent_all <- optim(fn = gpd_dtrunc, 
                          par = 1.4, 
                          method = "Brent",
                          lower = 1,
                          upper = 3,
                          hessian = TRUE,
                          dat = datu, 
                          supp = supp, 
                          slow = slow, 
                          expo = TRUE)
opt_francent_female <- optim(fn = gpd_dtrunc, 
                          par = 1.4, 
                          method = "Brent",
                          lower = 1, 
                          upper = 3, 
                          hessian = TRUE,
                          dat = datu[gender == "female"], 
                          supp = supp[gender == "female"], 
                          slow = slow[gender == "female"],
                          expo = TRUE)
opt_francent_male <- optim(fn = gpd_dtrunc, 
                          par = 1.4, 
                          method = "Brent",
                          lower = 1, 
                          upper = 3, 
                          hessian = TRUE,
                          dat = datu[gender == "male"], 
                          supp = supp[gender == "male"], 
                          slow = slow[gender == "male"],
                          expo = TRUE)
mle_france_exp <- c(opt_francent_all$par,
                    opt_francent_female$par,
                    opt_francent_male$par)
nsamp_france <- c(length(gender), sum(gender == "female"), sum(gender == "male"))
stdev_france_exp <- sqrt(1/c(opt_francent_all$hessian,
                      opt_francent_female$hessian,
                      opt_francent_male$hessian))
lci_france_exp <- mle_france_exp - qnorm(0.975)*stdev_france_exp
uci_france_exp <- mle_france_exp + qnorm(0.975)*stdev_france_exp

load("IDL2016.rda")
# idl2016 <- idl2016[!idl2016$countrydeath == "FRA",]
u110 <- min(idl2016$numdays)-1 #some people are 2-3 days from 110...
# Calendar date at which individual reaches 110 years
datu <- idl2016$numdays - u110
slow <- pmax(0, (idl2016$slow - u110)/365.25)
supp <- (idl2016$supp - u110)/365.25
datu <- datu/365.25
gender <- factor(idl2016$sex, labels = c("female","male"))
opt_IDL_all <- optim(fn = gpd_dtrunc, 
                          par = 1.4, 
                          method = "Brent",
                          lower = 1,
                          upper = 3,
                          hessian = TRUE,
                          dat = datu, 
                          supp = supp, 
                          slow = slow, 
                          expo = TRUE)
opt_IDL_female <- optim(fn = gpd_dtrunc, 
                             par = 1.4, 
                             method = "Brent",
                             lower = 1, 
                             upper = 3, 
                             hessian = TRUE,
                             dat = datu[gender == "female"], 
                             supp = supp[gender == "female"], 
                             slow = slow[gender == "female"],
                             expo = TRUE)
opt_IDL_male <- optim(fn = gpd_dtrunc, 
                           par = 1.4, 
                           method = "Brent",
                           lower = 1, 
                           upper = 3, 
                           hessian = TRUE,
                           dat = datu[gender == "male"], 
                           supp = supp[gender == "male"], 
                           slow = slow[gender == "male"],
                           expo = TRUE)
mle_IDL_exp <- c(opt_IDL_all$par,
                    opt_IDL_female$par,
                    opt_IDL_male$par)
nsamp_IDL <- c(length(gender), sum(gender == "female"), sum(gender == "male"))
stdev_IDL_exp <- sqrt(1/c(opt_IDL_all$hessian,
                             opt_IDL_female$hessian,
                             opt_IDL_male$hessian))
lci_IDL_exp <- mle_IDL_exp - qnorm(0.975)*stdev_IDL_exp
uci_IDL_exp <- mle_IDL_exp + qnorm(0.975)*stdev_IDL_exp

tab <- cbind(paste0("$",nsamp_istat,"$"),
      paste0("$", sprintf(fmt = "%.2f", mle_istat_exp), " (",
             sprintf(fmt = "%.2f", lci_istat_exp), ", ",
             sprintf(fmt = "%.2f", uci_istat_exp), ")$"),
      paste0("$",nsamp_france,"$"),
      paste0("$", sprintf(fmt = "%.2f", mle_france_exp), " (",
             sprintf(fmt = "%.2f", lci_france_exp), ", ",
             sprintf(fmt = "%.2f", uci_france_exp), ")$"),
      paste0("$",nsamp_IDL,"$"),
      paste0("$", sprintf(fmt = "%.2f", mle_IDL_exp), " (",
             sprintf(fmt = "%.2f", lci_IDL_exp), ", ",
             sprintf(fmt = "%.2f", uci_IDL_exp), ")$"))
tab <- as.data.frame(tab[c(2,3,1),])
row.names(tab) <- c("women","men","all")
colnames(tab) <- rep(c("$n$","$\\sigma_e$ ($95$\\% CI)"), length.out = 6)
if(tables){
  
  print(xtable(tab),
        file = "Table2.tex",
        booktabs = TRUE, sanitize.text.function = identity)
}
