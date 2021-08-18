#########################################################################################
# Power plots (Figure 3, right)
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
save <- TRUE #FALSE turn this switch to run the power study

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

frex <- data.frame(datu = datu, slow = slow, supp = supp)
profile_francent <- prof_gpd_dtrunc_xi(xi = xis,
                                       dat = frex$datu, 
                                       supp = frex$supp, 
                                       slow = frex$slow)$param




# load("IDL2016.rda")
# idl2016 <- idl2016[!idl2016$countrydeath == "FRA",]
# # Lower truncation level sometimes higher than the excess lifetime (in days)?
# idl2016 <- idl2016[-which(idl2016$slow > idl2016$numdays),]
# u110 <- min(idl2016$numdays)-1 #some people are 2-3 days from 110...
# # Calendar date at which individual reaches 110 years
# datu <- idl2016$numdays - u110
# slow <- pmax(0, (idl2016$slow - u110)/365.25)
# supp <- (idl2016$supp - u110)/365.25
# datu <- datu/365.25
# idlex <- data.frame(datu = datu, slow = slow, supp = supp)
load(idl2021.rda)

bootsamp <- matrix(0, ncol = length(datu), nrow = B)
prof_xi_idl <- prof_gpd_dtrunc_xi(xi = xis,
                                  dat = idlex$datu, 
                                  supp = idlex$supp, 
                                  slow = idlex$slow)$param
profile_idl <- prof_xi_idl
if(save){
  set.seed(20200824)
  prbar <- progress_bar$new(total = B*nxis,
                            format = "(:spin) [:bar] :percent",
                            clear = FALSE)
  power_combo <- array(0, dim = c(B, nxis, 2))
  power_italcent <- array(0, dim = c(B, nxis, 2))
  power_idl <- array(0, dim = c(B, nxis, 2))
  power_francent <- array(0, dim = c(B, nxis, 2))
  shapes <- array(0, dim = c(B, 4, nxis, 2))
  bootsampital <- matrix(0, ncol = length(itex$datu), nrow = B)
  bootsampfr <- matrix(0, ncol = length(frex$datu), nrow = B)
  bootsampidl <- matrix(0, ncol = length(idlex$datu), nrow = B)
  rightcensb <- matrix(FALSE, ncol = length(itex$datu), nrow = B)
  
  combll <- function(par, datidl, datfr, datital, rcensit, maxidl, maxfr, maxital){
    par_it <- par[2:1]
    par_fr <- par[c(3,1)]
    par_idl <- par[c(4,1)]
    gpd_dtrunc(par_fr, dat = datfr, slow = frex$slow, supp = frex$supp, expo = FALSE) +
      gpd_dtrunc(par_idl, dat = datidl, slow = idlex$slow, supp = idlex$supp, expo = FALSE) +  
      gpd_cens(par_it, dat = datital, rightcens = rcensit, slow = itex$slow, expo = FALSE)
  }
  ineq <- function(par, datidl, datfr, datital, rcensit, maxidl, maxfr, maxital){
    if(par[1] < 0){
      ubi <- c(maxital, maxfr, maxidl)-par[2:4]/par[1]
    } else{
      ubi <- rep(0.02, 3)
    }
    c(ubi, par - c(-1-1e8,1e-8,1e-8,1e-8), c(1.5,100,100,100) - par)
  }
  
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
      prbar$tick()
      #Find maximum likelihood estimates
      mleboot_combo <- try(alabama::auglag(par = c(profile_italcent[i,2:1], profile_francent[i,1], profile_idl[i,1]),
                                           fn = combll,
                                           hin = ineq,
                                           maxital = max(bootsampital[b,]),
                                           maxfr = max(bootsampfr[b,]),
                                           maxidl = max(bootsampidl[b,]),
                                           datital = bootsampital[b,],
                                           datfr = bootsampfr[b,],
                                           datidl = bootsampidl[b,],
                                           rcensit = rightcensb[b,],
                                           control.outer=list(method="nlminb", trace=0)))
      
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
      shapes[b,4,i,] <- c(mleboot_combo$par[1], sqrt(diag(solve(mleboot_combo$hessian))[1]))
      
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
      
      if(is.character(mleboot_combo)){
        power_combo[b,i,] <- rep(NA, 2)
      } else{
        #Compute directed profile likelihood ratio root test statistic
        power_combo[b, i, 1] <- sign(mleboot_combo$par[1])*
          sqrt(-2*(mleboot_combo$value - (pexpboot_fr$value + pexpboot_it$value + pexpboot_idl$value)))
        #Compute Wald-test
        if(mleboot_combo$par[2] > -0.5){
          power_combo[b,i,2] <- mleboot_combo$par[1]/sqrt(diag(solve(mleboot_combo$hessian))[1])
        } else{ 
          power_combo[b,i,2] <- NA
        }
      }
    }
  }
  prbar$terminate()
  save(xis, 
       power_italcent, 
       power_idl, 
       power_combo, 
       power_francent, 
       shapes, file = "power.RData")
} else{
  load("power.RData") 
}
