#########################################################################################
# Power plots (Figure 3, left)
#########################################################################################
# Compute the power for the endpoint
set.seed(20200618)
B <- 1e4

u108 <- 39447L #108 years
endpoints <- c(seq(116.2, 125, by = 0.2), seq(126, 150, by = 1))
nendpoints <- length(endpoints)
save <- TRUE

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
c1 <- lubridate::dmy("01-01-1987")
c2 <- lubridate::dmy("31-12-2017")

# Lower truncation level, zero if individual reached 110 between c1 and c2
francent$slow <- as.numeric(pmax(0, c1 - xcal))
francent$supp <- as.numeric(c2 - xcal)
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


# load("IDL2016.rda")
# idl2016 <- idl2016[!idl2016$countrydeath == "FRA",]
# u110 <- min(idl2016$numdays)-1 #some people are 2-3 days from 110...
# # Calendar date at which individual reaches 105 years
# datu <- idl2016$numdays - u110
# slow <- pmax(0, (idl2016$slow - u110)/365.25)
# supp <- (idl2016$supp - u110)/365.25
# datu <- datu/365.25
# maxidl <- max(idl2016$numdays/365.25)
# idlex <- data.frame(datu = datu, slow = slow, supp = supp)

load(IDL2021.rda)
maxidl <- max(idlex$datu) + 110

profile_idl <- t(sapply(1:nendpoints, function(i){
  prof_gpd_dtrunc_endpoint(endpoint = (endpoints[i] - 110), 
                           dat = idlex$datu, slow = idlex$slow, supp = idlex$supp)}))[,2:3]


if(save){
  prbar <- progress_bar$new(total = nendpoints,
                            format = "(:spin) [:bar] :percent",
                            clear = FALSE)
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
    prbar$tick()
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
    save(endpoints, 
         power_ep_italcent, 
         power_ep_idl, 
         power_ep_combo, 
         power_ep_francent, 
         file = "power_endpoint.RData")
  }
  prbar$terminate()
  save(endpoints, 
       power_ep_italcent, 
       power_ep_idl, 
       power_ep_combo, 
       power_ep_francent, 
       file = "power_endpoint.RData")
} 
