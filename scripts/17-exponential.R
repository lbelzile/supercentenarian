# Table 2: parameter estimates for exponential
# distribution above 108 (Istat, France) and 110 (IDL, rest)
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

ital_df <- data.frame(datu = datu, 
                      slow = slow, 
                      gender = gender,
                      rcens = rcens)
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
nsamp_istat <- c(length(gender), 
                 sum(gender == "female"),
                 sum(gender == "male"))
stdev_istat_exp <- mle_istat_exp/sqrt(nsamp_istat)
lci_istat_exp <- mle_istat_exp - qnorm(0.975)*stdev_istat_exp
uci_istat_exp <- mle_istat_exp + qnorm(0.975)*stdev_istat_exp


load("francent.rda")
# Calendar date at which individual reaches 108 years
xcal <- francent$birth + u108
# Calendar time for sampling frame
c1 <- lubridate::dmy("01-01-1987")
c2 <- lubridate::dmy("31-12-2017")

# Lower truncation level, zero if individual reached 105 between c1 and c2
francent$slow <- as.numeric(pmax(0, c1 - xcal))
francent$supp <- as.numeric(c2 - xcal)
datu <- francent$numdays - u108
ind <- which(datu > 0)
slow <- pmax(0, francent$slow[ind]/365.25)
supp <- francent$supp[ind]/365.25
datu <- (francent$numdays[ind] - u108)/365.25
gender <- francent$gender[ind]
fran_df <- data.frame(datu = datu, 
                      slow = slow, 
                      supp = supp,
                      gender = gender)
opt_francent_all <- optim(fn = gpd_dtrunc, 
                          par = 1.4, 
                          method = "Brent",
                          lower = 0.1,
                          upper = 6,
                          hessian = TRUE,
                          dat = datu, 
                          supp = supp, 
                          slow = slow, 
                          expo = TRUE)
opt_francent_female <- optim(fn = gpd_dtrunc, 
                          par = 1.4, 
                          method = "Brent",
                          lower = 0.1, 
                          upper = 6, 
                          hessian = TRUE,
                          dat = datu[gender == "female"], 
                          supp = supp[gender == "female"], 
                          slow = slow[gender == "female"],
                          expo = TRUE)
opt_francent_male <- optim(fn = gpd_dtrunc, 
                          par = 1.4, 
                          method = "Brent",
                          lower = 0.1, 
                          upper = 6, 
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

# Estimates for French women, excluding Jeanne Calment
opt_francent_female_mJC <- optim(fn = gpd_dtrunc, 
                             par = 1.4, 
                             method = "Brent",
                             lower = 0.1, 
                             upper = 6, 
                             hessian = TRUE,
                             dat = datu[gender == "female" & datu != max(datu)], 
                             supp = supp[gender == "female" & datu != max(datu)], 
                             slow = slow[gender == "female" & datu != max(datu)],
                             expo = TRUE)

# Compute likelihood ratio-based confidence intervals
# First, compute the log likelihood on a grid of values of sigma
psi_grid <- seq(0.5,3, by = 0.01)
confint_exp_francent <- 
  list(psi = psi_grid,
       pll = -sapply(psi_grid, function(scale){
         gpd_dtrunc(
           par = scale, 
       dat = datu, 
       supp = supp, 
       slow = slow, 
       expo = TRUE)}),
       psi.max = opt_francent_all$par,
       maxpll = -opt_francent_all$value
  )
confint_exp_francent_female <- 
  list(psi = psi_grid,
       pll = - sapply(psi_grid, function(scale){
         gpd_dtrunc(
           par = scale, 
           dat = datu[gender == "female"], 
           supp = supp[gender == "female"], 
           slow = slow[gender == "female"],
           expo = TRUE)}),
       psi.max = opt_francent_female$par,
       maxpll = -opt_francent_female$value
  )
confint_exp_francent_male <- 
  list(psi = seq(0.5,3, by = 0.01),
       pll = - sapply(seq(0.5,3, by = 0.01), function(scale){
         gpd_dtrunc(
           par = scale, 
           dat = datu[gender == "male"], 
           supp = supp[gender == "male"], 
           slow = slow[gender == "male"],
           expo = TRUE)}),
       psi.max = opt_francent_male$par,
       maxpll = -opt_francent_male$value
  )
# Compute the intervals by root search
conf_fr <- confint_int(object = confint_exp_francent, parm = "profile")
conf_fr_f <- confint_int(object = confint_exp_francent_female, parm = "profile")
conf_fr_m <- confint_int(object = confint_exp_francent_male, parm = "profile")

# Confidence interval for scale
round(conf_fr, 2)
# Confidence interval for hazard (reciprocal scale)
round(1/conf_fr, 2)
# Confidence interval for prob
round(exp(-1/conf_fr), 2) 
# Hazard ratio male/female
opt_francent_female$par/opt_francent_male$par

load("IDL2021.rda")
idl_df <- idlex


# load("IDL2016.rda")
# # idl2016 <- idl2016[!idl2016$countrydeath == "FRA",]
# u110 <- min(idl2016$numdays) #some people are 2-3 days from 110...
# # Calendar date at which individual reaches 110 years
# datu <- idl2016$numdays - u110
# ind <- datu > 0
# slow <- pmax(0, (idl2016$slow[ind] - u110)/365.25)
# supp <- (idl2016$supp[ind] - u110)/365.25
# datu <- datu[ind]/365.25
# gender <- factor(idl2016$sex[ind], labels = c("female","male"))
# idl_df <- data.frame(datu = datu, 
#                       slow = slow, 
#                       supp = supp,
#                       gender = gender)[idl2016$countrydeath[ind] != "FRA",]
opt_IDL_all <- with(idl_df,
                    optim(fn = gpd_dtrunc, 
                          par = 1.4, 
                          method = "Brent",
                          lower = 0.1,
                          upper = 5,
                          hessian = TRUE,
                          dat = datu, 
                          supp = supp, 
                          slow = slow, 
                          expo = TRUE))
opt_IDL_female <- with(idl_df, 
                       optim(fn = gpd_dtrunc, 
                             par = 1.4, 
                             method = "Brent",
                             lower = 0.1, 
                             upper = 5, 
                             hessian = TRUE,
                             dat = datu[gender == "female"], 
                             supp = supp[gender == "female"], 
                             slow = slow[gender == "female"],
                             expo = TRUE))
opt_IDL_male <- with(idl_df, 
                     optim(fn = gpd_dtrunc, 
                           par = 1.4, 
                           method = "Brent",
                           lower = 0.1, 
                           upper = 5, 
                           hessian = TRUE,
                           dat = datu[gender == "male"], 
                           supp = supp[gender == "male"], 
                           slow = slow[gender == "male"],
                           expo = TRUE))
mle_IDL_exp <- c(opt_IDL_all$par,
                    opt_IDL_female$par,
                    opt_IDL_male$par)
nsamp_IDL <- with(idl_df, 
                  c(length(gender), 
                    sum(gender == "female"), 
                    sum(gender == "male")))
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
  setwd(table_dir)
  print(xtable::xtable(tab),
        file = "Table2.tex",
        booktabs = TRUE, 
        sanitize.text.function = identity)
  setwd(code_dir)
}

# Pooled estimates combining three database
# but excluding French from IDL2016 to avoid
# overlaps
pooled_loglik_exp <- 
  function(par, france = NULL, italy = NULL, idl = NULL){
    l1 <- l2 <- l3 <- 0
    if(!is.null(france)){
   l1 <- gpd_dtrunc( 
        par = par,
        dat = france$datu, 
        supp = france$supp, 
        slow = france$slow,
        expo = TRUE) 
    }
    if(!is.null(idl)){
   l2 <- gpd_dtrunc( 
         par = par,
         dat = idl$datu, 
         supp = idl$supp, 
         slow = idl$slow,
         expo = TRUE) 
    }
    if(!is.null(italy)){
  l3 <- gpd_cens( 
       par = par,
       dat = italy$datu, 
       rightcens = italy$rcens, 
       slow = italy$slow,
       expo = TRUE)
    }
   l1 + l2 + l3
  }

# Compute maximum likelihood estimate
mle_pooled <- optim(par = 1.4, 
      fn = pooled_loglik_exp,
      lower = 0.5, 
      upper = 3, 
      method = "Brent",
      france = fran_df,
      italy = ital_df,
      idl = idl_df)

# Same, without French men
mle_pooled_2 <- optim(par = 1.4, 
                    fn = pooled_loglik_exp,
                    lower = 0.5, 
                    upper = 3, 
                    method = "Brent",
                    france = fran_df[fran_df$gender == "female",],
                    italy = ital_df,
                    idl = idl_df)
# This time without Jeanne Calment
mle_pooled_3 <- optim(par = 1.4, 
                      fn = pooled_loglik_exp,
                      lower = 0.5, 
                      upper = 3, 
                      method = "Brent",
                      france = fran_df[!which.max(fran_df$datu),],
                      italy = ital_df,
                      idl = idl_df)

# Likelihood ratio confidence intervals
# These are invariant to reparametrization
confint_exp_pooled <- 
  list(psi = seq(0.5,3, by = 0.01),
       pll = - sapply(seq(0.5,3, by = 0.01), function(scale){
         pooled_loglik_exp(scale, 
                           france = fran_df,
                           italy = ital_df,
                           idl = idl_df)}),
       psi.max = mle_pooled$par,
       maxpll = -mle_pooled$value
  )

conf <- confint_int(object = confint_exp_pooled, 
                    parm = "profile")
# 95% confidence intervals for the scale
round(conf, 2)
# For probability of survival one year
round(exp(-1/conf), 2)

# Difference between scales excluding some French
mle_pooled_2$par - mle_pooled$par
mle_pooled_3$par - mle_pooled$par


# Pooling only IDL and IStat
# 
mle_pooled_4 <- optim(par = 1.4, 
                      fn = pooled_loglik_exp,
                      lower = 0.5, 
                      upper = 3, 
                      method = "Brent",
                      france = NULL,
                      italy = ital_df,
                      idl = idl_df)
confint_exp_pooled_4 <- 
list(psi = seq(1, 3, by = 0.01),
     pll = - sapply(seq(1,3, by = 0.01), function(scale){
       pooled_loglik_exp(scale, 
                         france = NULL,
                         italy = ital_df,
                         idl = idl_df)}),
     psi.max = mle_pooled_4$par,
     maxpll = -mle_pooled_4$value
)


mle_pooled_5 <- optim(par = 1.4, 
                      fn = pooled_loglik_exp,
                      lower = 0.5, 
                      upper = 3, 
                      method = "Brent",
                      france = NULL,
                      italy = ital_df[ital_df$gender == "male",],
                      idl = idl_df[idl_df$gender == "male",])
mle_pooled_6 <- optim(par = 1.4, 
                      fn = pooled_loglik_exp,
                      lower = 0.5, 
                      upper = 3, 
                      method = "Brent",
                      france = NULL,
                      italy = ital_df[ital_df$gender == "female",],
                      idl = idl_df[idl_df$gender == "female",])

# Confidence interval for the Istat + IDL
conf <- confint_int(object = confint_exp_pooled_4, 
                    parm = "profile")

# Check that pooling doesn't lead to a significant decrease in fit
lrt_exp_dbase <- 2*(mle_pooled_4$value - (opt_IDL_all$value + pooled_loglik_exp(par = mle_istat_exp[1], italy = ital_df)))
pchisq(q = lrt_exp_dbase, 
       df = 1, 
       lower.tail = FALSE)

lrt_exp_sex <- 2*(mle_pooled_4$value - (mle_pooled_5$value + mle_pooled_6$value))
pchisq(q = lrt_exp_sex, 
       df = 1, 
       lower.tail = FALSE)
