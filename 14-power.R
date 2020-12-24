#########################################################################################
# Power plots (Figure 3)
#########################################################################################
# Power of GP for H0: xi>=0 versus Ha:xi < 0 (finite endpoint)
# Three database: IDL, France (2019) and IStat
# We condition for all three on entry points
# Right-censor for Istat data, else simulate doubly truncated data

# Compute both directed likelihood root and Wald statistics and compare power

# Power analysis for one-sided test with Ha: gamma>=0 versus H0: gamma < 0
# and for two-sided tests with Ha gamma=0 versus H0: gamma != 0 (not shown in text)

set.seed(20200618)
B <- 1e4
library(progress)

u108 <- 39447L #108 years
xis <- seq(-0.25, 0.25, length.out = 101L)
nxis <- length(xis)
save <- FALSE # Uncomment to run the power study

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
datu <- datu[ind]/365.25
slow <- italcent$slow[ind]/365.25

itex <- data.frame(datu = datu, slow = slow, rcens = rcens, uplim = uplim)
profile_italcent <- t(sapply(1:nxis, function(i){
  prof_gpd_cens_xi(xi = xis[i], 
                   dat = itex$datu, 
                   rightcens = itex$rcens, 
                   slow = itex$slow)}))[,3:2]

load("francent.rda")
# Calendar date at which individual reaches 105 years
xcal <- francent$birth + u108
# Calendar time for sampling frame
c1a <- lubridate::dmy("01-01-1978")
c1b <- lubridate::dmy("01-01-1987")
c2 <- lubridate::dmy("31-12-2017")

# Lower truncation level, zero if individual reached 110 between c1 and c2
francent$slow <- ifelse(francent$numdays > 110*365.24, 
                        as.numeric(pmax(0, c1b - xcal)),
                        as.numeric(pmax(0, c1a - xcal)))
francent$supp <- as.numeric(pmax(0, c2 - xcal))    
datu <- francent$numdays - u108
ind <- which(datu > 0)
slow <- francent$slow[ind]/365.25
supp <- francent$supp[ind]/365.25
datu <- (francent$numdays - u108)[ind]/365.25

frex <- data.frame(datu = datu, slow = slow, supp = supp)
profile_francent <- prof_gpd_dtrunc_xi(xi = xis,
                                       dat = frex$datu, 
                                       supp = frex$supp, 
                                       slow = frex$slow)$param



load("IDL2016.rda")
idl2016 <- idl2016[!idl2016$countrydeath == "FRA",]
u110 <- min(idl2016$numdays)-1 #some people are 2-3 days from 110...
# Calendar date at which individual reaches 105 years
datu <- idl2016$numdays - u110
slow <- pmax(0, (idl2016$slow - u110)/365.25)
supp <- (idl2016$supp - u110)/365.25
datu <- datu/365.25

idlex <- data.frame(datu = datu, slow = slow, supp = supp)
# Why is lower truncation level sometimes higher than the excess lifetime (in days)?

bootsamp <- matrix(0, ncol = length(datu), nrow = B)
prof_xi_idl <- prof_gpd_dtrunc_xi(xi = xis,
                                  dat = idlex$datu, 
                                  supp = idlex$supp, 
                                  slow = idlex$slow)$param
profile_idl <- prof_xi_idl





if(save){
  set.seed(20200824)
  power_italcent <- array(0, dim = c(B, nxis, 2))
  power_idl <- array(0, dim = c(B, nxis, 2))
  power_francent <- array(0, dim = c(B, nxis, 2))
  shapes <- array(0, dim = c(B, 3, nxis, 2))
  bootsampital <- matrix(0, ncol = length(itex$datu), nrow = B)
  bootsampfr <- matrix(0, ncol = length(frex$datu), nrow = B)
  bootsampidl <- matrix(0, ncol = length(idlex$datu), nrow = B)
  rightcensb <- matrix(FALSE, ncol = length(itex$datu), nrow = B)
  
  
  for(i in 1:nxis){
    #Generate new data
    for(j in 1:ncol(bootsampfr)){
      bootsampfr[,j] <- revddtrunc(n = B, param = profile_francent[i,], 
                                   lower = frex$slow[j], 
                                   upper = frex$supp[j])
    }
    for(j in 1:ncol(bootsampidl)){
      bootsampidl[,j] <- revddtrunc(n = B, param = profile_idl[i,], 
                                    lower = frex$slow[j], 
                                    upper = frex$supp[j])
    }
    for(j in 1:ncol(bootsampital)){
      bootsampital[,j] <- evd::qgpd(evd::pgpd(itex$slow[j], loc = 0, scale = profile_italcent[i,1], shape = profile_italcent[i,2]) + runif(B)*
                                      (1-evd::pgpd(itex$slow[j], loc = 0, scale = profile_italcent[i,1], shape = profile_italcent[i,2])), 
                                    loc = 0, scale = profile_italcent[i,1], shape = profile_italcent[i,2])
      rightcensb[,j] <- bootsampital[,j] > itex$uplim[j]
      bootsampital[rightcensb[,j],j] <- itex$uplim[j]
    }
    for(b in 1:B){
      #Find maximum likelihood estimates
      
      mleboot_it <- try(optim(par =  profile_italcent[i,],
                              fn = gpd_cens,
                              method = "N",
                              dat = bootsampital[b,],
                              rightcens = rightcensb[b,],
                              slow = itex$slow,
                              expo = FALSE,
                              control = list(fnscale = c(1.5,0.01), reltol = 1e-12, maxit = 1e5),
                              hessian = TRUE))
      mleboot_fr <- try(optim(par =  profile_francent[i,],
                              fn = gpd_dtrunc,
                              method = "N",
                              dat = bootsampfr[b,],
                              supp = frex$supp,
                              slow = frex$slow,
                              expo = FALSE,
                              control = list(fnscale = c(1.5,0.01), reltol = 1e-12, maxit= 1e5),
                              hessian = TRUE))
      mleboot_idl <- try(optim(par =  profile_idl[i,],
                               fn = gpd_dtrunc,
                               method = "N",
                               dat = bootsampidl[b,],
                               slow = idlex$slow,
                               supp = idlex$supp,
                               expo = FALSE,
                               control = list(fnscale = c(1.5,0.01), reltol = 1e-12, maxit= 1e5),
                               hessian = TRUE))
      
      shapes[b,1,i,] <- c(mleboot_it$par[2], sqrt(diag(solve(mleboot_it$hessian))[2]))
      shapes[b,2,i,] <- c(mleboot_fr$par[2], sqrt(diag(solve(mleboot_fr$hessian))[2]))
      shapes[b,3,i,] <- c(mleboot_idl$par[2], sqrt(diag(solve(mleboot_idl$hessian))[2]))
      
      mbootfr <- mean(bootsampfr[b,])
      mbootidl <- mean(bootsampidl[b,])
      pexpboot_idl <- try(optim(par = mbootidl, 
                                fn = gpd_dtrunc, 
                                method = "Brent", 
                                lower = 0.01*mbootidl, 
                                upper = 3*mbootidl, 
                                dat = bootsampidl[b,],
                                supp = idlex$supp, 
                                slow = idlex$slow, 
                                expo = TRUE))
      pexpboot_fr <- try(optim(par = mbootfr, 
                               fn = gpd_dtrunc, 
                               method = "Brent", 
                               lower = 0.01*mbootfr, 
                               upper = 3*mbootfr, 
                               dat = bootsampfr[b,],
                               supp = frex$supp, 
                               slow = frex$slow, 
                               expo = TRUE))
      rate_it <- exp_mle_lt_rc(dat =  bootsampital[b,], 
                               rightcens = rightcensb[b,], 
                               slow = itex$slow)
      pexpboot_it <- list(par = rate_it, 
                          value = gpd_cens(par = c(rate_it, 0), 
                                           dat =  bootsampital[b,], 
                                           rightcens = rightcensb[b,], 
                                           slow = itex$slow))
      
      if(is.character(mleboot_it)){
        power_italcent[b,i,] <- rep(NA, 2)
      } else{
        #Compute directed  profile likelihood ratio root test statistic
        power_italcent[b, i, 1] <- sign(mleboot_it$par[2])*
          sqrt(2*(pexpboot_it$value - mleboot_it$value))
        #Compute Wald-test
        if(mleboot_it$par[2] > -0.5){
          power_italcent[b,i,2] <- mleboot_it$par[2]/sqrt(diag(solve(mleboot_it$hessian))[2])
        } else{ 
          power_italcent[b,i,2] <- NA
        }
      }
      
      if(is.character(mleboot_fr) || is.character(pexpboot_fr)){
        power_francent[b,i,] <- rep(NA, 2)
      } else{
        #Compute directed profile likelihood ratio root test statistic
        power_francent[b, i, 1] <- sign(mleboot_fr$par[2])*
          sqrt(2*(pexpboot_fr$value - mleboot_fr$value))
        #Compute Wald-test
        if(mleboot_fr$par[2] > -0.5){
          power_francent[b,i,2] <- mleboot_fr$par[2]/sqrt(diag(solve(mleboot_fr$hessian))[2])
        } else{ 
          power_francent[b,i,2] <- NA
        }
      }
      
      if(is.character(mleboot_idl) || is.character(pexpboot_idl)){
        power_idl[b,i,] <- rep(NA, 2)
      } else{
        #Compute directed profile likelihood ratio root test statistic
        power_idl[b, i, 1] <- sign(mleboot_idl$par[2])*
          sqrt(-2*(mleboot_idl$value - pexpboot_idl$value))
        #Compute Wald-test
        if(mleboot_idl$par[2] > -0.5){
          power_idl[b,i,2] <- mleboot_idl$par[2]/sqrt(diag(solve(mleboot_idl$hessian))[2])
        } else{ 
          power_idl[b,i,2] <- NA
        }
      }
    }
  }
  save(xis, power_italcent, power_idl, power_francent, shapes, file = "power.RData")
} else{
  load("power.RData") 
}

# Power for the endpoint

set.seed(20200618)
B <- 1e4
#Set directory from Rstudio to current document
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
#setwd("/home/lbelzile/Documents/Dropbox/SuperCentenarian/code")
source("Semi-supercentenarian_fn.R")
figures <- TRUE
fig_dir <- paste0(getwd(), "/figure")
library(progress)

u108 <- 39447L #108 years
endpoints <- c(seq(116.2, 125, by = 0.2), seq(126, 150, by = 1))
nendpoints <- length(endpoints)
save <- FALSE

load("italcent.rda")
# Calendar date at which individual reaches 108 years
italcent$xcal <- italcent$birth + u108
# Calendar time for sampling frame
c1 <- lubridate::dmy("01-01-2009")
c2 <- lubridate::dmy("01-01-2016")

# Lower truncation level, zero if individual reached 108 between c1 and c2
italcent$slow <- as.numeric(pmax(0, c1 - italcent$xcal))
maxital <- max(italcent$numdays/365.25)
# Exceedances
datu <- italcent$numdays - u108
# Fit a generalized Pareto model
ind <- which(datu > 0)
uplim <- as.numeric(c2 - italcent$birth)
# Keep only exceedances above 108
uplim <- (uplim[ind] - u108)/365.25
rcens <- italcent$rightcens[ind]
datu <- datu[ind]/365.25
slow <- italcent$slow[ind]/365.25

itex <- data.frame(datu = datu, slow = slow, rcens = rcens, uplim = uplim)
profile_italcent <- t(sapply(1:nendpoints, function(i){
  prof_gpd_cens_endpoint(endpoint = (endpoints[i] - 108), 
                         dat = itex$datu, 
                         rightcens = itex$rcens, 
                         slow = itex$slow)}))[,2:3]

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
maxfran <- max(francent$numdays/365.25)
frex <- data.frame(datu = datu, slow = slow, supp = supp)
profile_francent <- t(sapply(1:nendpoints, function(i){
  prof_gpd_dtrunc_endpoint(endpoint = (endpoints[i] - 108), 
                           dat = frex$datu, slow = frex$slow, supp = frex$supp)}))[,2:3]


load("IDL2016.rda")
idl2016 <- idl2016[!idl2016$countrydeath == "FRA",]
u110 <- min(idl2016$numdays)-1 #some people are 2-3 days from 110...
# Calendar date at which individual reaches 105 years
datu <- idl2016$numdays - u110
slow <- pmax(0, (idl2016$slow - u110)/365.25)
supp <- (idl2016$supp - u110)/365.25
datu <- datu/365.25
maxidl <- max(idl2016$numdays/365.25)
idlex <- data.frame(datu = datu, slow = slow, supp = supp)
# Why is lower truncation level sometimes higher than the excess lifetime (in days)?
profile_idl <- t(sapply(1:nendpoints, function(i){
  prof_gpd_dtrunc_endpoint(endpoint = (endpoints[i] - 110), 
                           dat = idlex$datu, slow = idlex$slow, supp = idlex$supp)}))[,2:3]


if(save){
  set.seed(20200824)
  power_ep_combo <- array(NA, dim = c(B, nendpoints))
  power_ep_francent <- array(NA, dim = c(B, nendpoints))
  power_ep_idl <- array(NA, dim = c(B, nendpoints))
  power_ep_italcent <- array(NA, dim = c(B, nendpoints))
  
  
  bootsampital <- matrix(0, ncol = length(itex$datu), nrow = B)
  bootsampfr <- matrix(0, ncol = length(frex$datu), nrow = B)
  bootsampidl <- matrix(0, ncol = length(idlex$datu), nrow = B)
  rightcensb <- matrix(FALSE, ncol = length(itex$datu), nrow = B)
  
  combllep <- function(par,  maxfr, maxidl,  maxital, datidl, datfr, datital, rcensit){
    endpoint <- par[1]
    shape <- par[2:4]
    scale <- -shape*(endpoint - c(108,108,110))
    par_it <- c(scale[1], shape[1])
    par_fr <- c(scale[2], shape[2])
    par_idl <- c(scale[3], shape[3])
    ub <- endpoint - c(108,108,110)
    gpd_dtrunc(par_fr, dat = datfr, slow = frex$slow, supp = frex$supp, expo = FALSE) +
      gpd_dtrunc(par_idl, dat = datidl, slow = idlex$slow, supp = idlex$supp, expo = FALSE) +  
      gpd_cens(par_it, dat = datital, rightcens = rcensit, slow = itex$slow, expo = FALSE)
  }
  combllep_scale <- function(par,  maxfr, maxidl,  maxital, datidl, datfr, datital, rcensit){
    endpoint <- par[1]
    scale <- par[2:4]
    shape <- -scale/(endpoint - c(108,108,110))
    par_it <- c(scale[1], shape[1])
    par_fr <- c(scale[2], shape[2])
    par_idl <- c(scale[3], shape[3])
    ub <- endpoint - c(108,108,110)
    gpd_dtrunc(par_fr, dat = datfr, slow = frex$slow, supp = frex$supp, expo = FALSE) +
      gpd_dtrunc(par_idl, dat = datidl, slow = idlex$slow, supp = idlex$supp, expo = FALSE) +  
      gpd_cens(par_it, dat = datital, rightcens = rcensit, slow = itex$slow, expo = FALSE)
  }
  ineq <- function(par, ...){
    args <- list(...)
    c(par[1]-pmax(args$maxital, args$maxfr, args$maxidl), par[2:4] + 1, -par[2:4])}
  ineq_scale <- function(par, ...){
    args <- list(...)
    endpoint <- par[1]
    scale <- par[2:4]
    shape <- -scale/(endpoint - c(108,108,110))
    c(endpoint-pmax(args$maxital, args$maxfr, args$maxidl), shape + 1, -shape)}
  
  for(i in 1:nendpoints){
    #Generate new data
    if(endpoints[i] > maxfran){
      for(j in 1:ncol(bootsampfr)){
        bootsampfr[,j] <- revddtrunc(n = B, param = profile_francent[i,], 
                                     lower = frex$slow[j], 
                                     upper = frex$supp[j])
      }
    }
    if(endpoints[i] > maxidl){
      for(j in 1:ncol(bootsampidl)){
        bootsampidl[,j] <- revddtrunc(n = B, param = profile_idl[i,], 
                                      lower = frex$slow[j], 
                                      upper = frex$supp[j])
      }
    }
    if(endpoints[i] > maxital){
      for(j in 1:ncol(bootsampital)){
        bootsampital[,j] <- evd::qgpd(evd::pgpd(itex$slow[j], loc = 0, scale = profile_italcent[i,1], shape = profile_italcent[i,2]) + runif(B)*
                                        (1-evd::pgpd(itex$slow[j], loc = 0, scale = profile_italcent[i,1], shape = profile_italcent[i,2])), 
                                      loc = 0, scale = profile_italcent[i,1], shape = profile_italcent[i,2])
        rightcensb[,j] <- bootsampital[,j] > itex$uplim[j]
        bootsampital[rightcensb[,j],j] <- itex$uplim[j]
      }
    }
    for(b in 1:B){
      #Find maximum likelihood estimates
      if(endpoints[i] > maxfran){
        maxital = max(bootsampital[b,])
        maxfr = max(bootsampfr[b,])
        maxidl = max(bootsampidl[b,])
        mleboot_combo <- try(alabama::auglag(par = c(pmax(maxital, maxfr, maxidl, endpoints[i]) + 0.5, profile_italcent[i,1], profile_francent[i,1], profile_idl[i,1]),
                                             fn = combllep_scale,
                                             hin = ineq_scale,
                                             maxital = maxital,
                                             maxfr = maxfr,
                                             maxidl = maxidl,
                                             datital = bootsampital[b,],
                                             datfr = bootsampfr[b,],
                                             datidl = bootsampidl[b,],
                                             rcensit = rightcensb[b,],
                                             control.outer=list(method="BFGS", trace=0),
                                             control.optim = list(parscale = c(endpoints[i], profile_italcent[i,1], profile_francent[i,1], profile_idl[i,1]))))
      }
      if(endpoints[i] > maxital){
        mleboot_it <- try(optim(par =  profile_italcent[i,],
                                fn = gpd_cens,
                                method = "N",
                                dat = bootsampital[b,],
                                rightcens = rightcensb[b,],
                                slow = itex$slow,
                                expo = FALSE,
                                control = list(fnscale = c(1.5,0.01), reltol = 1e-12, maxit = 1e5),
                                hessian = TRUE))
        rate_it <- exp_mle_lt_rc(dat =  bootsampital[b,], 
                                 rightcens = rightcensb[b,], 
                                 slow = itex$slow)
        pexpboot_it <- list(par = rate_it, 
                            value = gpd_cens(par = c(rate_it, 0), 
                                             dat =  bootsampital[b,], 
                                             rightcens = rightcensb[b,], 
                                             slow = itex$slow))
        if(!is.character(mleboot_it)){
          #Compute directed profile likelihood ratio root test statistic
          power_ep_italcent[b, i] <- sign(mleboot_it$par[2])*sqrt(2*(pexpboot_it$value - mleboot_it$value))
        }
        
      }
      if(endpoints[i] > maxfran){
        mleboot_fr <- try(optim(par =  profile_francent[i,],
                                fn = gpd_dtrunc,
                                method = "N",
                                dat = bootsampfr[b,],
                                supp = frex$supp,
                                slow = frex$slow,
                                expo = FALSE,
                                control = list(fnscale = c(1.5,0.01), reltol = 1e-12, maxit= 1e5),
                                hessian = TRUE))
        mbootfr <- mean(bootsampfr[b,])
        pexpboot_fr <- try(optim(par = mbootfr, 
                                 fn = gpd_dtrunc, 
                                 method = "Brent", 
                                 lower = 0.01*mbootfr, 
                                 upper = 3*mbootfr, 
                                 dat = bootsampfr[b,],
                                 supp = frex$supp, 
                                 slow = frex$slow, 
                                 expo = TRUE))
        if(!(is.character(mleboot_fr) || is.character(pexpboot_fr))){
          #Compute directed profile likelihood ratio root test statistic
          power_ep_francent[b, i] <- sign(mleboot_fr$par[2])*sqrt(2*(pexpboot_fr$value - mleboot_fr$value))
        }
      }
      if(endpoints[i] > maxidl){
        mleboot_idl <- try(optim(par =  profile_idl[i,],
                                 fn = gpd_dtrunc,
                                 method = "N",
                                 dat = bootsampidl[b,],
                                 slow = idlex$slow,
                                 supp = idlex$supp,
                                 expo = FALSE,
                                 control = list(fnscale = c(1.5,0.01), reltol = 1e-12, maxit= 1e5),
                                 hessian = TRUE))
        mbootidl <- mean(bootsampidl[b,])
        pexpboot_idl <- try(optim(par = mbootidl, 
                                  fn = gpd_dtrunc, 
                                  method = "Brent", 
                                  lower = 0.01*mbootidl, 
                                  upper = 3*mbootidl, 
                                  dat = bootsampidl[b,],
                                  supp = idlex$supp, 
                                  slow = idlex$slow, 
                                  expo = TRUE))
        if(!(is.character(mleboot_idl) || is.character(pexpboot_idl))){
          power_ep_idl[b, i] <- sign(mleboot_idl$par[2])*sqrt(-2*(mleboot_idl$value - pexpboot_idl$value))
        }
      }
      if(endpoints[i] > maxfran){
        if(!any(is.character(mleboot_combo), 
                is.character(pexpboot_idl), 
                is.character(pexpboot_fr))){
          power_ep_combo[b, i] <- -sqrt(pmax(0,-2*(mleboot_combo$value - (pexpboot_fr$value + pexpboot_it$value + pexpboot_idl$value))))
        }
      } 
    }
  }
  save(endpoints, power_ep_italcent, power_ep_idl, power_ep_combo, power_ep_francent, file = "power_endpoint.RData")
} else{
  load("power_endpoint.RData") 
}
powerIDL_ep <- colMeans(power_ep_idl < qnorm(0.05), na.rm = TRUE)
powerIstat_ep <- colMeans(power_ep_italcent < qnorm(0.05), na.rm = TRUE)
powerFrance_ep <- colMeans(power_ep_francent < qnorm(0.05), na.rm = TRUE)
powerCombo_ep <- colMeans(power_ep_combo < qnorm(0.05), na.rm = TRUE)


# Smoothing of the power curve to get rid of the artifacts of Monte-Carlo variability
library(cobs)
library(ggplot2)
library(poorman)
library(viridis)
library(patchwork)
library(tikzDevice)
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
             constraint = "decrease",
             nknots = nknots),
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

load("power.RData") 
nxis <- length(xis)
power_italcent[,xis == 0,]

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
  xlab("upper limit to human lifespan") 


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
round(rbind(
  powersmoothIstat[which(endpts %in% c(125, 130, 135))],
  powersmoothFrance[which(endpts %in% c(125, 130, 135))],
  powersmoothIDL[which(endpts %in% c(125, 130, 135))],
  powersmoothCombo[which(endpts %in% c(125, 130, 135))]
),2)