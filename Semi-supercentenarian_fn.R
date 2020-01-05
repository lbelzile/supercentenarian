#' Likelihood for left-truncated and right-censored generalized Pareto variates
#'
#' Computes the log-likelihood for generalized Pareto or exponential observations.
#'
#' @param par vector of scale and shape
#' @param dat vector of threshold exceedances
#' @param rightcens logical indicating right-censoring (\code{TRUE} for censored)
#' @param slow lower truncation limit
#' @param expo logical; should an exponential model be fitted? Default to \code{FALSE}
#' @return log-likelihood value
gpd_cens <- function(par, dat, rightcens, slow, expo = FALSE){
  if(expo){
    stopifnot(length(par) == 1L)
    shape <- 0
  } else{
    shape <- par[2]
  }
  if(par[1] < 0 || shape < (-1+1e-8)){
    return(1e10)
  }
  g1 <- intersect(which(!rightcens), which(slow > 0))
  g2 <- intersect(which(rightcens), which(slow > 0))
  ll <- 0
  #Contribution from observations in sampling frame (density divided by truncation interval)
  if(sum(!rightcens)>0){
    ll <- sum(eva::dgpd(loc = 0, scale = par[1], shape = shape, x = dat[!rightcens], log = TRUE))
    if(length(g1) > 0){
      ll <- ll - sum(log(1-eva::pgpd(slow[g1], loc = 0, scale = par[1], shape = shape)))  #right censored individuals
    }
  }
  if(sum(rightcens)>0){
    ll <- ll +  sum(log(1-eva::pgpd(dat[rightcens], loc = 0, scale = par[1], shape = shape)))
    if(length(g2) > 0){
      ll <- ll - sum(log(1-eva::pgpd(slow[g2], loc = 0, scale = par[1], shape = shape)))
    }
  }
  if (!is.finite(ll)) {
    return(1e10)
  }  else {
    return(-ll)
  }
}


#' Observed information matrix 
#' 
#' Observed information matrix for left-truncated right-censored generalized Pareto observations 
#' @inheritParams gpd_cens
#' @param theta vector of endpoing and shape parameters for the generalized Pareto
#' @return a two by two matrix
obs.infomat <- function(theta, dat, rightcens, slow){
  endpt = theta[1]; xi = theta[2];
  if(endpt < 0){
    return(matrix(rep(Inf, 4), nrow = 2))
  }
  j11 <- (endpt^2*xi - dat^2 + 2*dat*endpt)/((dat^2*endpt^2 - 2*dat*endpt^3 + endpt^4)*xi)
  j12 <- -dat/((dat*endpt - endpt^2)*xi^2)
  j22 <- (xi - 2*log(-(dat - endpt)/endpt))/xi^3
  
  jc11 <- -(dat^2 - 2*dat*endpt)/((dat^2*endpt^2 - 2*dat*endpt^3 + endpt^4)*xi)
  jc12 <- -dat/((dat*endpt - endpt^2)*xi^2)
  jc22 <- -2*log(-(dat - endpt)/endpt)/xi^3
  
  jt11 <- (slow^2 - 2*slow*endpt)/((slow^2*endpt^2 - 2*slow*endpt^3 + endpt^4)*xi)
  jt12 <- slow/((slow*endpt - endpt^2)*xi^2)
  jt22 <- 2*log(-(slow - endpt)/endpt)/xi^3
  
  info11 <- sum(j11[!rightcens])+sum(jc11[rightcens])+sum(jt11)
  info12 <- sum(j12[!rightcens])+sum(jc12[rightcens])+sum(jt12)
  info22 <- sum(j22[!rightcens])+sum(jc22[rightcens])+sum(jt22)
  -cbind(c(info11, info12), c(info12, info22))
}

#' Profile likelihood for endpoint
#' 
#' Profile likelihood for the endpoint of the generalized Pareto distribution
#' for left-truncted and right-censored observations.
#' 
#' @inheritsParam gpd_cens
#' @param psi value of the endpoint at which to compute the profile
#' @return a vector of length three containing twice the negative log-likelihood value, the endpoint value and the maximum of the nuisance lambda (i.e., the shape parameter).
prof_gpd_cens <- function(psi, dat, rightcens, slow){
  opt <- optim(fn = function(lambda){gpd_cens(par = c(-lambda*psi, lambda),
                                              dat = dat,  rightcens = rightcens, slow = slow)},
               method = "Brent", par = -0.01, lower = -0.98, upper = -1e-3, control = list(reltol=1e-12))  
  res <- -2*opt$value
  attributes(res) <- list("param" = c(psi, opt$par))
  #res
  return(c(res, psi, opt$par))
}


#' Likelihood root function
#' 
#' This function returns the likelihood root of the profile log-likelihood for the endpoint of the generalized Pareto distribution with left-truncated and right-censored data.
#' Specifically, \eqn{-r^2/2} is the profile likelihood and the two-sided p-value is\code{qchisq(p, 1)/2}.
#' 
#' @param psi parameter of the endpoint at which to compute the p-value
#' @param thetahat maximum likelihood estimates of the endpoint and the shape parameters
#' @inheritParams gpd_cens
#' @return a p-value
rfun <- function(psi, thetahat, dat, rightcens, slow){
  # Log-likelihood function in terms of parameters (endpt, xi)
  loglik <- function(theta, dat, rightcens, slow){
    gpd_cens(par = c(-theta[2]*theta[1], theta[2]), dat = dat, rightcens = rightcens, slow = slow)
  }
  llp <- prof_gpd_cens(psi = psi, dat = dat, rightcens = rightcens, slow = slow)
  sign(thetahat[1]-llp[2])*sqrt(-2*loglik(thetahat, dat = dat, rightcens = rightcens, slow = slow) - llp[1])
}


#' Confidence bands for Nelson-Aalen estimator
#' 
#' This function is adapted from code by David Diez
#' in OIsurv::confBands, licensed under GPL (>= 2)

cumhazbands <- function(time, n, cumhaz, type = c("ptwise", "band"), sigma.sq, transfo = c("none", "log", "arcsin"), confLevel = c(0.9, 0.95, 0.99)) {
  transfo <- match.arg(transfo)
  type <- match.arg(type)
  t.U <- last(time)
  t.L <- first(time)
  n <- length(time)
  if(type == "band"){
    a.L <- n * sigma.sq[1]/(1 + n * sigma.sq[1])
    a.U <- n * sigma.sq[n]/(1 + n * sigma.sq[n])
    aU <- format(c(round(50 * a.U)/50, 0.01))[1]
    aL <- format(c(round(50 * a.L)/50, 0.01))[1]
    if (confLevel[1] == 0.9) {
      test <- try(data(hw.k10))
      if(is.matrix(test)){
        input <- hw.k10[aU, aL]
      } else{
        input <- OIsurv:::hw.k10 [aU, aL]
      }
    }
    else if (confLevel[1] == 0.95) {
      test <- try(data(hw.k05))
      if(is.matrix(test)){
        input <- hw.k05[aU, aL]
      } else{
        input <- OIsurv:::hw.k05 [aU, aL]
      }
    }
    else if (confLevel[1] == 0.99) {
      test <- try(data(hw.k01))
      if(is.matrix(test)){
        input <- hw.k01[aU, aL]
      } else{
        input <- OIsurv:::hw.k01 [aU, aL]
      }
    }  else {
      stop("Only confidence levels of 0.90, 0.95 and 0.99\nare allowed.")
    }
    m.fact <- input * (1 + n * sigma.sq)/sqrt(n)
  } else{
    m.fact <- qnorm((1+confLevel[1])/2)*sqrt(sigma.sq)
  }
  CI <- matrix(NA, length(time), 2)
  if(transfo == "none"){
    CI[, 1] <- cumhaz - m.fact
    CI[, 2] <- cumhaz + m.fact
  } else if(transfo == "log"){
    CI[, 1] <- cumhaz*exp(-m.fact*sqrt(sigma.sq)/cumhaz)
    CI[, 2] <- cumhaz*exp(m.fact*sqrt(sigma.sq)/cumhaz)
  } else if(transfo == "arcsin"){
    CI[,1] <- -2*log(sin(pmin(pi/2, asin(exp(-cumhaz/2))+0.5*m.fact*sqrt(sigma.sq/(exp(cumhaz)-1)))))
    CI[,2] <- -2*log(sin(pmax(0, asin(exp(-cumhaz/2))-0.5*m.fact*sqrt(sigma.sq/(exp(cumhaz)-1)))))                         
  }
  CI[CI[, 1] < 0, 1] <- 0
  tR <- list(time = time, lower = CI[, 1], upper = CI[, 2])
  class(tR) <- "confBands"
  return(tR)
}

#' Maximum likelihood estimate of scale
#' 
#' Maximum likelihood estimate of the shape for the exponential distribution with left-truncated right-censored observations
#' @inheritParams gpd_cens
#' @return scale parameter estimtae
exp_mle_lt_rc <- function(dat, slow, rightcens){
  sum(dat-slow)/ sum(!rightcens) 
}

#' Confidence intervals based on (modified) profile likelihoods
#'
#' This code is adapted from the mev package (mev:::confint.eprof)
#' @param object a list containing informations about the profile likelihood in the same format as the \code{hoa} package
#' @param parm string or numerical vector giving the type of interval to consider
#' @param level probability level of the confidence interval
#' @param prob vector of length 2 containing the bounds, by default double-sided
#' @param print logical indicating whether the intervals are printed to the console
#' @param ... additional arguments passed to the function
#' @return a table with confidence intervals.
confint_int <- function (object, parm, level = 0.95, prob = c((1 - level)/2,  1 - (1 - level)/2), print = FALSE, ...) 
{
  if (!isTRUE(all.equal(diff(prob), level, check.attributes = FALSE))) {
    warning("Incompatible arguments: `level` does not match `prob`.")
  }
  args <- list(...)
  if ("warn" %in% names(args) && is.logical(args$warn)) {
    warn <- args$warn
  }  else {
    warn <- TRUE
  }
  if (length(prob) != 2) {
    stop("`prob` must be a vector of size 2")
    prob <- sort(prob)
  }
  if (is.numeric(parm)) {
    ind <- parm
    parm <- c("profile", "tem", "modif.tem", "modif.empcov")[ind]
  } else {
    parm <- match.arg(arg = parm, choices = c("profile", 
                                              "tem", "modif.tem", "modif.empcov", "r", "rstar"), 
                      several.ok = TRUE)
    parm[parm %in% "r"] <- "profile"
    parm[parm %in% "rstar"] <- "tem"
    ind <- which(c("profile", "tem", "modif.tem", "modif.empcov") %in% 
                   parm)
  }
  parm <- unique(parm)
  ind <- unique(ind[ind %in% 1:4])
  if (length(ind) == 0) {
    stop("Invalid `parm` argument.")
  }
  qulev <- qnorm(1 - prob)
  conf <- matrix(ncol = 4, nrow = 3)
  i = 1
  if (is.null(object$pll) && is.null(object$r)) {
    break
  }
  if (is.null(object$r)) {
    object$r <- sign(object$psi.max - object$psi) * 
      sqrt(2 * (object$maxpll - object$pll))
  } else {
    object$r[is.infinite(object$r)] <- NA
  }
  if (is.null(object$normal)) {
    object$normal <- c(object$psi.max, object$std.error)
  }
  if (requireNamespace("cobs", quietly = TRUE)) {
    fit.r <- cobs::cobs(x = object$r, y = object$psi, 
                        constraint = "decrease", lambda = 0, ic = "SIC", 
                        pointwise = cbind(0, 0, object$normal[1]), 
                        knots.add = TRUE, repeat.delete.add = TRUE, 
                        print.mesg = FALSE, print.warn = FALSE)
    pr <- predict(fit.r, c(0, qulev))[, 2]
  } else {
    fit.r <- stats::smooth.spline(x = na.omit(cbind(object$r, 
                                                    object$psi)), cv = FALSE)
    pr <- predict(fit.r, c(0, qulev))$y
    pr[1] <- object$normal[1]
  }
  conf[, i] <- pr
  if (warn) {
    if (!any(object$r > qnorm(prob[1]))) {
      warning("Extrapolating the lower confidence interval for the profile likelihood ratio test")
    }
    if (!any(object$r < qnorm(prob[2]))) {
      warning("Extrapolating the upper confidence interval for the profile likelihood ratio test")
    }
  }
  
  if (!is.null(conf)) {
    colnames(conf) <- c("Profile", "TEM", "Severini (TEM)", 
                        "Severini (emp. cov.)")
    rownames(conf) <- c("Estimate", "Lower CI", "Upper CI")
    wrong_below <- which(conf[2, ] > conf[1, ])
    if (length(wrong_below) > 0) {
      conf[2, ][wrong_below] <- NA
    }
    wrong_above <- which(conf[3, ] < conf[1, ])
    if (length(wrong_above) > 0) {
      conf[3, ][wrong_above] <- NA
    }
    if (print) {
      cat("Point estimate for the parameter of interest psi:\n")
      cat("Maximum likelihood          :", round(object$psi.max, 
                                                 3), "\n")
      cat("\n")
      cat("Confidence intervals, levels :", prob, "\n")
      cat("Wald intervals               :", round(object$psi.max + 
                                                    sort(qulev) * object$std.error, 3), "\n")
      cat("Profile likelihood           :", round(conf[2:3, 
                                                       1], 3), "\n")
    }
    return(invisible(conf[, ind]))
  }
}


#' Profile likelihood for shape parameter of the generalized Pareto distribution
#' 
#' #' This function implements the profile likelihood with left-truncated and right-censored observations
#' and must be applied repeatedly for each threshold
#' @inheritParams gpd_cens
#' @return point estimate and confidence limits for the shape parameter
prof_gpd_cens_xi_confint <- function(dat, rightcens, slow){
  xis <- seq(-0.5, 0.5, length = 200)
  mdat <- max(dat)
  dev <- sapply(xis, function(xi){
    opt <- optimize(f = function(lambda){gpd_cens(par = c(lambda, xi),
                                                  dat = dat,  rightcens = rightcens, slow = slow)},
                    interval = c(ifelse(xi < 0, mdat*abs(xi), 1e-8), 1e9), tol = 1e-10)
    c(-2*opt$objective, opt$minimum)
  })
  ms <- which.max(dev[1,])
  opt_mle <- optim(f = function(theta){gpd_cens(par = theta, dat = dat,  rightcens = rightcens, slow = slow)},
                   par = c(dev[2,ms],xis[ms]), method = "N", hessian = TRUE, 
                   control = list(parscale = c(500, 0.01), reltol = 1e-10, maxit = 1000))
  
  mle <- opt_mle$par       
  maxll <- -2*opt_mle$value
  prof <- list(psi = xis, pll = dev[1,]-maxll, maxpll = 0, mle = mle, 
               psi.max = mle[2], std.error = sqrt(solve(opt_mle$hessian)[2,2]))
  confint <- confint_int(prof, parm = "profile")
  confint
  # return(c(res, psi, opt$minimum))
}

#' Profile likelihood for the scale of the exponential distribution 
#' 
#' This function implements the profile likelihood with left-truncated and right-censored observations
#' and must be applied repeatedly for each threshold
#' @inheritParams gpd_cens
#' @return point estimate and confidence limits for the scale parameter
prof_exp_cens <- function(dat, rightcens, slow){
  opt_mle <- optim(par = 1.5, fn = gpd_cens, method = "Brent", lower = 0.1, upper = 3, dat = dat,
                   rightcens = rightcens, slow = slow, expo = TRUE, hessian = TRUE)
  grid_theta <- c(opt_mle$par, seq(0.4, 3, length = 250L))
  dev <- -2*sapply(grid_theta, function(theta){
    gpd_cens(par = theta, dat = dat, rightcens = rightcens, slow = slow, expo = TRUE)})
  confint_int(list(psi = sort(grid_theta), pll = dev[order(grid_theta)], 
                   maxpll = -2*opt_mle$value, psi.max = opt_mle$par, 
                   std.error = sqrt(solve(opt_mle$hessian)[1,1]),
                   mle = c(opt_mle$par,sqrt(solve(opt_mle$hessian)[1,1]))),parm = "profile")
}

#' Profile deviance of exponential sub-model
#' 
#' This function computes the profile when the shape is zero and returns
#' the value of twice the negative log-likelihoodat that value
#' @inheritParams gpd_cens
#' @return twice the negative profile log-likelihood
nllhxizero <- function(dat, rightcens, slow){
  opt <- optimize(f = function(lambda){gpd_cens(par = c(lambda, 0),
                                                dat = dat,  rightcens = rightcens, slow = slow)},
                  interval = c(1e-8, 1e9), tol = 1e-10)
  c(-2*opt$objective)
}



#' Local generalized Pareto distribution
#'
#' This function breaks down the likelihood contribution from the sample in intervals
#' so that the total contribution is only counted once
#' @inheritParams gpd_cens
#' @param low lower bound under which data are censored or left-truncated
#' @param upp upper bound beyond which observations are right-censored
#' @return negative of local likelihood
gpd_intcens <- function(par, dat, rightcens, slow, low, upp, expo = FALSE){
  ab <- as.logical(I(dat > low)*I(slow < upp)) #survive beyond age under consideration
  datu <- dat[ab]
  slowu <- pmax(low, slow[ab])
  rcens <- rightcens[ab]
  aboveupp <- datu > upp
  if(length(aboveupp) > 0){
    rcens[aboveupp] <- TRUE 
  }
  if(expo){
    stopifnot(length(par) == 1L)
    shape <- 0
  } else{
    shape <- par[2]
  }
  if(par[1] < 0 || shape < (-1+1e-8)){
    return(1e10)
  }
  ll <- 0
  #Contribution from observations in sampling frame (density divided by truncation interval)
  if(sum(!rcens) > 0){
    ll <- ll + sum(dgpd(loc = 0, scale = par[1], shape = shape, x = datu[!rcens], log = TRUE))
  }
  if(sum(rcens) > 0){
    ll <- ll + sum(log(pgpd(lower.tail = FALSE, loc = 0, scale = par[1], shape = shape, q = pmin(datu[rcens], rep(upp, sum(rcens))))))
  }
  ll <- ll -  sum(log(pgpd(loc = 0, scale = par[1], shape = shape, q = slowu, lower.tail = FALSE)))
  if (!is.finite(ll)) {
    return(1e10)
  }  else {
    return(-ll)
  }
}

#' Profile likelihood for hazard of generalized Pareto
#' 
#' This function returns the profile likelihood-based 95% confidence intervals 
#' of the hazard of the generalized Pareto distribution with left-truncated and right-censored observations.
#' @inheritParams gpd_cens
#' @return point estimate along with 95% confidence interval
prof_gpd_hazard_confint <- function(dat, rightcens, slow, thresh){
  confintmat <- matrix(0, nrow = length(thresh), ncol = 4)
  th <- thresh
  k <- 200L
  datd <- dat/365.25
  slowd <- slow/365.25
  b <- 0
  for(i in 1:length(thresh)){
    b <- b + 1L
    dth <- (th[i] - th[1])
    low <- (th[i] - th[1])
    upp <- ifelse(i > length(thresh)-1, 20, (th[i] - th[1])+1)
    par_est <- optim(par = c(0.7,-0.3), fn = function(x){gpd_intcens(par = c(1/x[1]-x[2]*dth, x[2]), 
                                                                     dat = datd, 
                                                                     rightcens = rightcens, low = low, upp = upp,
                                                                     slow = slowd, expo = FALSE)}, control = list(reltol = 1e-11))
    hazard_mle <- par_est$par[1]
    infoi <- numDeriv::hessian(x = par_est$par,
                               func =  function(x){gpd_intcens(par = c(1/x[1]-x[2]*dth, x[2]), 
                                                               dat = datd, 
                                                               rightcens = rightcens, low = low, upp = upp,
                                                               slow = slowd, expo = FALSE)}) 
    stderr.transfo <- sqrt(diag(solve(infoi)))[1]
    if(stderr.transfo < 1e-6){
      grid_psi <-  par_est$par[1] + seq(-0.4, 0.4, length = k)  
    } else{
      grid_psi <-  par_est$par[1] + seq(-3*stderr.transfo, 4*stderr.transfo, length = k)  
    }
    
    prof_vals <- xi_sigma_vals <- rep(0, k)
    for (j in 1:k) {
      opt_prof <- optimize(f = function(xi, haza, dat, slow, rightcens, dth){
        gpd_intcens(par = c(1/haza-xi*dth, xi), dat = dat, rightcens = rightcens, 
                    slow = slow, low = low, upp = upp)},
        upper = min(1.1,1/(grid_psi[j]*dth)), lower = -1,
        haza = grid_psi[j], dat = datd, 
        rightcens = rightcens, slow = slowd, dth = dth)
      xi_sigma_vals[j] <- opt_prof$minimum
      prof_vals[j] <- opt_prof$objective
    }
    infin <- which(prof_vals  > 1e10-1)
    if(length(infin) > 0){
      prof_vals[infin] <- NA
    }
    prof <- structure(list(psi = grid_psi, psi.max = hazard_mle[1], 
                           pll = -prof_vals, maxpll = -par_est$value, std.err = stderr.transfo), 
                      class = "eprof")
    confintmat[b,] <- c(confint_int(prof, parm = "profile",print = FALSE)[1:3], stderr.transfo)
  }
  return(confintmat)
}


#' Likelihood for left-truncated and right-censored generalized Pareto variates
#'
#' Computes the log-likelihood for generalized Pareto or exponential observations.
#'
#' @param par vector of scale and shape
#' @param dat vector of threshold exceedances
#' @param supp upper truncation limit
#' @param slow lower truncation limit
#' @param expo logical; should an exponential model be fitted? Default to \code{FALSE}
#' @return log-likelihood value
gpd_dtrunc <- function(par, dat, slow, supp, expo = FALSE){
  if(expo){
    stopifnot(length(par) == 1L)
    shape <- 0
  } else{
    shape <- par[2]
  }
  if(par[1] < 0 || shape < (-1+1e-8)){
    return(1e10)
  }
  ll <- 0
  #Contribution from observations in sampling frame (density divided by truncation interval)
  ll <- sum(eva::dgpd(loc = 0, scale = par[1], shape = shape, x = dat, log = TRUE))
  if(shape < 0){
    ll <- ll - sum(log(eva::pgpd(pmin(-par[1]/shape-1e-12, supp), loc = 0, scale = par[1], shape = shape) - eva::pgpd(slow, loc = 0, scale = par[1], shape = shape))) 
  } else {
    ll <- ll - sum(log(eva::pgpd(supp, loc = 0, scale = par[1], shape = shape) - eva::pgpd(slow, loc = 0, scale = par[1], shape = shape))) 
  }
  if (!is.finite(ll)) {
    return(1e10)
  }  else {
    return(-ll)
  }
}
#' Fit doubly truncated generalized Pareto distribution
#' 
#' @param thresh vector of thresholds over which to apply the function
#' @inheritsParam gpd_dtrunc
fitgpddtrunc <- function(dat, thresh, supp, slow, expo = FALSE){
  res <- sapply(seq_along(thresh), function(i){
    #For each threshold, compute new threshold exceedances
    datu <- dat - thresh[i]
    # Keep only exceedances
    ind <- which(datu > 0)
    if(expo){
      vals <- optim(par = 1, fn = gpd_dtrunc, 
                    method = "Brent", dat = (datu[ind])/365.25, lower = 0.01, upper = 100,
                    supp = (supp[ind]-thresh[i])/365.25, 
                    slow = (pmax(0, slow[ind]-thresh[i]))/365.25, expo = TRUE, 
                    control = list(reltol = 1e-10, maxit= 1e5),
                    hessian = TRUE)
    } else{
      vals <- optim(par = fit.gpd(datu/365.25)$est, fn = gpd_dtrunc, 
                    method = "N", dat = (datu[ind])/365.25,
                    supp = (supp[ind]-thresh[i])/365.25, 
                    slow = (pmax(0, slow[ind]-thresh[i]))/365.25, expo = FALSE, 
                    control = list(parscale = c(1,0.01), 
                                   reltol = 1e-10, maxit= 1e5),
                    hessian = TRUE)
    }
    c(loc = (u + thresh[i])/365.25, vals$par, 
      -2*vals$value, sqrt(diag(solve(vals$hessian))),
      length(ind))
    
  })
  res <- t(as.matrix(res))
  colnames(res) <- switch(expo,
                          c("loc","scale","deviance", "scale.stderr","nu"),
                          c("loc","scale","shape","deviance","scale.stderr","shape.stderr","nu"))
  rownames(res) <- round((u + thresh)/365.25, 1)
  return(res)
}


#' Profile likelihood of doubly truncated generalized Pareto distribution
#' 
#' This function returns 95% profile likelihood-based confidence intervals for 
#' the shape parameter of the generalized Pareto distribution
#' @inheritParams gpd_dtrunc
#' @return estimate and 95% confidence interval for the shape parameter
prof_gpd_dtrunc_xi_confint <- function(dat, supp, slow){
  xis <- seq(-0.5, 0.5, length = 200)
  mdat <- max(dat)
  dev <- sapply(xis, function(xi){
    opt <- optimize(f = function(lambda){gpd_dtrunc(par = c(lambda, xi),
                                                    dat = dat,  supp = supp, slow = slow)},
                    interval = c(ifelse(xi < 0, mdat*abs(xi), 1e-8), 1e9), tol = 1e-10)
    c(-2*opt$objective, opt$minimum)
  })
  ms <- which.max(dev[1,])
  opt_mle <- optim(f = function(theta){gpd_dtrunc(par = theta, dat = dat,  supp = supp, slow = slow)},
                   par = c(dev[2,ms],xis[ms]), method = "N", hessian = TRUE, 
                   control = list(parscale = c(500, 0.01), reltol = 1e-10, maxit = 1000))
  
  mle <- opt_mle$par       
  maxll <- -2*opt_mle$value
  prof <- list(psi = xis, pll = dev[1,]-maxll, maxpll = 0, mle = mle, 
               psi.max = mle[2], std.error = sqrt(solve(opt_mle$hessian)[2,2]))
  confint <- confint_int(prof, parm = "profile")
  confint
}

#' Log-likelihood of doubly truncated exponential distribution
#' 
#' This function returns 95% likelihood-based confidence intervals for 
#' the scale parameter of the exponential distribution
#' @inheritParams gpd_dtrunc
#' @return estimate and 95% confidence interval for the scale parameter
prof_exp_dtrunc <- function(dat, supp, slow){
  opt_mle <- optim(par = 1.5, fn = gpd_dtrunc, method = "Brent", lower = 0.1, upper = 3, dat = dat,
                   supp = supp, slow = slow, expo = TRUE, hessian = TRUE)
  grid_theta <- c(opt_mle$par, seq(0.4, 3, length = 250L))
  dev <- -2*sapply(grid_theta, function(theta){
    gpd_dtrunc(par = theta, dat = dat, supp = supp, slow = slow, expo = TRUE)})
  confint_int(list(psi = sort(grid_theta), pll = dev[order(grid_theta)], 
                   maxpll = -2*opt_mle$value, psi.max = opt_mle$par, 
                   std.error = sqrt(solve(opt_mle$hessian)[1,1]),
                   mle = c(opt_mle$par,sqrt(solve(opt_mle$hessian)[1,1]))),parm = "profile")
}

