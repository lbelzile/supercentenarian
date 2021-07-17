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
    ll <- sum(evd::dgpd(loc = 0, scale = par[1], shape = shape, x = dat[!rightcens], log = TRUE))
    if(length(g1) > 0){
      ll <- ll - sum(log(1-evd::pgpd(slow[g1], loc = 0, scale = par[1], shape = shape)))  #right censored individuals
    }
  }
  if(sum(rightcens)>0){
    ll <- ll +  sum(log(1-evd::pgpd(dat[rightcens], loc = 0, scale = par[1], shape = shape)))
    if(length(g2) > 0){
      ll <- ll - sum(log(1-evd::pgpd(slow[g2], loc = 0, scale = par[1], shape = shape)))
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
#' for left-truncated and right-censored observations.
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
  ll <- sum(evd::dgpd(loc = 0, scale = par[1], shape = shape, x = dat, log = TRUE))
  if(shape < 0){
    ll <- ll - sum(log(evd::pgpd(pmin(-par[1]/shape-1e-12, supp), loc = 0, scale = par[1], shape = shape) - evd::pgpd(slow, loc = 0, scale = par[1], shape = shape))) 
  } else {
    ll <- ll - sum(log(evd::pgpd(supp, loc = 0, scale = par[1], shape = shape) - evd::pgpd(slow, loc = 0, scale = par[1], shape = shape))) 
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
#' @param xis vector of shape parameters
#' @param confint logical; if \code{TRUE}, returns confidence interval rather than estimates
#' @return estimate and 95% confidence interval for the shape parameter
prof_gpd_dtrunc_xi <- function(dat, supp, slow, confint = FALSE, xis = NULL){
  if(is.null(xis)){
    xis <- seq(-0.5, 0.5, length = 200)
  }
  mdat <- max(dat)
  dev <- sapply(xis, function(xi){
    opt <- optimize(f = function(lambda){gpd_dtrunc(par = c(lambda, xi),
                                                    dat = dat,  supp = supp, slow = slow)},
                    interval = c(ifelse(xi < 0, mdat*abs(xi), 1e-8), 1e4), tol = 1e-10)
    c(-2*opt$objective, opt$minimum)
  })
  ms <- which.max(dev[1,])
  opt_mle <- optim(f = function(theta){
    gpd_dtrunc(par = theta, dat = dat,  supp = supp, slow = slow)},
                   par = c(dev[2,ms],xis[ms]), method = "N", hessian = TRUE, 
                   control = list(parscale = c(1, 0.01), reltol = 1e-10, maxit = 1000))
  
  mle <- opt_mle$par       
  maxll <- -2*opt_mle$value
  std_error <- try(sqrt(solve(opt_mle$hessian)[2,2]))
  if(is.character(std_error)){
    std_error <- NA
  }
  prof <- list(psi = xis, pll = dev[1,]-maxll, 
               param = cbind(dev[2,], xis),
               maxpll = 0, mle = mle, 
               psi.max = mle[2], std.error = std_error)
  if(confint){
    confint_int(prof, parm = "profile")
  } else{
    prof
  }
}

#' Sample observations from a doubly truncated generalized Pareto distribution
#' 
#' @param n sample size
#' @param param scale and shape parameters of the generalized Pareto distribution
#' @param lower scalar lower bound
#' @param upper scalar upper bound
#' @return a vector of \code{n} observations 
revddtrunc <- function(n, param, lower, upper){
  evd::qgpd(evd::pgpd(lower, scale = param[1], shape = param[2]) + 
              runif(n)*diff(evd::pgpd(c(lower, upper), scale = param[1], shape = param[2])),
            scale = param[1], shape = param[2])
}

#' Log-likelihood of doubly truncated exponential distribution
#' 
#' This function returns 95% likelihood-based confidence intervals for 
#' the scale parameter of the exponential distribution
#' @inheritParams gpd_dtrunc
#' @return estimate and 95% confidence interval for the scale parameter
prof_exp_dtrunc <- function(dat, supp, slow){
  opt_mle <- optim(par = 1.5, fn = gpd_dtrunc, method = "Brent", lower = 0.01*mean(dat), upper = 50*mean(dat), dat = dat,
                   supp = supp, slow = slow, expo = TRUE, hessian = TRUE)
  grid_theta <- c(opt_mle$par, seq(0.4*opt_mle$par, 3*opt_mle$par, length = 250L))
  dev <- -2*sapply(grid_theta, function(theta){
    gpd_dtrunc(par = theta, dat = dat, supp = supp, slow = slow, expo = TRUE)})
  confint_int(list(psi = sort(grid_theta), pll = dev[order(grid_theta)], 
                   maxpll = -2*opt_mle$value, psi.max = opt_mle$par, 
                   std.error = sqrt(solve(opt_mle$hessian)[1,1]),
                   mle = c(opt_mle$par,sqrt(solve(opt_mle$hessian)[1,1]))),
              parm = "profile")
}


#' Likelihood for left-truncated and right-censored generalized Gompertz variates
#'
#' Computes the log-likelihood for Gompertz distribution.
#'
#' @param par vector of scale and shape
#' @param dat vector of threshold exceedances
#' @param rightcens logical indicating right-censoring (\code{TRUE} for censored)
#' @param slow lower truncation limit
#' @param expo logical; should an exponential model be fitted? Default to \code{FALSE}
#' @return log-likelihood value
gomp_cens <- function(par, dat, rightcens, slow, expo = FALSE){
  scale <- par[1]
  shape <- par[2]
  if(scale <= 0 || shape <= 0){
    return(1e10)
  }
  # Density of Gompertz
  dgomp <- function(x, scale, shape, log = FALSE){
    ll <- -log(scale) + (1-exp(shape/scale*x))/shape + shape/scale*x
    if(log){
      return(ll)
    } else{
      return(exp(ll))
    }
  }
  # Survival function of Gompertz
  sgomp <- function(x, scale, shape){
    exp((1-exp(shape/scale*x))/shape)
  }
  
  g1 <- intersect(which(!rightcens), which(slow > 0))
  g2 <- intersect(which(rightcens), which(slow > 0))
  ll <- 0
  #Contribution from observations in sampling frame (density divided by truncation interval)
  if(sum(!rightcens)>0){
    ll <- sum(dgomp(scale = scale, shape = shape, x = dat[!rightcens], log = TRUE))
    if(length(g1) > 0){
      ll <- ll - sum(log(sgomp(x = slow[g1], scale = scale, shape = shape)))  #right censored individuals
    }
  }
  if(sum(rightcens)>0){
    ll <- ll +  sum(log(sgomp(x = dat[rightcens], scale = scale, shape = shape)))
    if(length(g2) > 0){
      ll <- ll - sum(log(sgomp(slow[g2], scale = scale, shape = shape)))
    }
  }
  if (!is.finite(ll)) {
    return(1e10)
  }  else {
    return(-ll)
  }
}

#' Likelihood for doubly-truncated generalized Gompertz variates
#'
#' Computes the log-likelihood for Gompertz distribution.
#'
#' @param par vector of scale and shape
#' @param dat vector of threshold exceedances
#' @param supp upper truncation limit
#' @param slow lower truncation limit
#' @param expo logical; should an exponential model be fitted? Default to \code{FALSE}
#' @return log-likelihood value
gomp_dtrunc <- function(par, dat, supp, slow, expo = FALSE){
  scale <- par[1]
  shape <- par[2]
  if(scale <= 0 || shape <= 0){
    return(1e10)
  }
  if(isTRUE(any(slow < 0, dat < slow, dat > supp))){
    return(1e10)
  }
  # Density of Gompertz
  dgomp <- function(x, scale, shape, log = FALSE){
    ll <- -log(scale) + (1-exp(shape/scale*x))/shape + shape/scale*x
    if(log){
      return(ll)
    } else{
      return(exp(ll))
    }
  }
  # Survival function of Gompertz
  sgomp <- function(x, scale, shape){
    exp((1-exp(shape/scale*x))/shape)
  }
  
  ll <- 0
  #Contribution from observations in sampling frame (density divided by truncation interval)
  ll <- sum(dgomp(scale = scale, shape = shape, x = dat, log = TRUE))
  ll <- ll - sum(log(sgomp(slow, scale = scale, shape = shape) - sgomp(supp, scale = scale, shape = shape)))
  if (!is.finite(ll)) {
    return(1e10)
  }  else {
    return(-ll)
  }
}



#' Profile likelihood for shape of generalized Pareto distribution
#'
#' Profile likelihood for the shape parameter of the generalized Pareto distribution
#' for left-truncated and right-censored observations.
#'
#' @inheritsParam gpd_cens
#' @param psi value of the shape at which to compute the profile
#' @return a vector of length three containing twice the negative log-likelihood value, the shape parameter and the maximum of the nuisance lambda (i.e., the scale parameter).
prof_gpd_cens_xi <- function(xi, dat, rightcens, slow){
  if(abs(xi) > 1e-6){
    opt <- optim(fn = function(lambda){gpd_cens(par = c(lambda, xi),
                                                dat = dat,  rightcens = rightcens, slow = slow)},
                 method = "Brent", par = 1, lower = ifelse(xi>=0, 1e-2, max(dat)*abs(xi)+1e-8) , upper = 15,
                 control = list(reltol=1e-12))
    res <- -2*opt$value
    return(c(res, xi, opt$par))
  } else{
    rate <- exp_mle_lt_rc(dat = dat, rightcens = rightcens, slow = slow)
    return(c(-2*gpd_cens(par = c(rate, 0), dat = dat,  rightcens = rightcens, slow = slow), 0, rate))
  }
  #attributes(res) <- list("param" = c(xi, opt$par))
  #res
}

#' Profile likelihood for the endpoint of the generalized Pareto distribution
#'
#' Profile likelihood for the endpoint of the generalized Pareto distribution
#' for left-truncated and right-censored observations.
#'
#' @inheritsParam gpd_cens
#' @param endpoint value of the endpoint at which to compute the profile
#' @return a vector of length three containing twice the negative log-likelihood value, the endpoint value and the maximum of the nuisance lambda (i.e., the shape parameter).
prof_gpd_cens_endpoint <- function(endpoint, dat, rightcens, slow){
  stopifnot(length(endpoint) == 1L)
    opt <- optim(fn = function(xi){
      sigma <- -endpoint*xi
      gpd_cens(par = c(sigma, xi), dat = dat, rightcens = rightcens, slow = slow)},
                 method = "Brent", par = -0.1, 
      lower = -1, 
      upper = -1e-8,
                 control = list(reltol=1e-12))
    res <- -2*opt$value
    return(c(res, -opt$par*endpoint, opt$par))
}

#' Profile likelihood for the endpoint of the generalized Pareto distribution
#'
#' Profile likelihood for the endpoint of the generalized Pareto distribution
#' for left-truncated and right-censored observations.
#'
#' @inheritsParam gpd_dtrunc
#' @param endpoint value of the endpoint at which to compute the profile
#' @return a vector of length three containing twice the negative log-likelihood value, the endpoint value and the maximum of the nuisance lambda (i.e., the shape parameter).
prof_gpd_dtrunc_endpoint <- function(endpoint, dat, slow, supp){
  stopifnot(length(endpoint) == 1L)
  opt <- optim(fn = function(xi){
    sigma <- -endpoint*xi
    gpd_dtrunc(par = c(sigma, xi), dat = dat, supp = supp, slow = slow)},
    method = "Brent", par = -0.1, 
    lower = -1, 
    upper = -1e-8,
    control = list(reltol=1e-12))
  res <- -2*opt$value
  return(c(res, -opt$par*endpoint, opt$par))
}

#' Sample lifetime of semi-supercentenarians
#' 
#' Given parameters of a generalized Pareto distribution, sampling window and 
#' birth dates with excess lifetimes, sample new observations; excess lifetime
#' at \code{c1} are sampled from an exponential distribution, whereas
#' the birth dates are samples from a jittered histogram-based distribution
#' The new excess lifetime above the threshold are right-censored if they exceed
#' \code{c2}.
#' 
#' @param n sample size
#' @param pars vector of length 2 containing the scale and shape parameters of the generalized Pareto distribution
#' @param xcal date at which individual reaches \code{u} years
#' @param c1 date, first day of the sampling frame
#' @param c2 date, last day of the sampling frame
#' @param slow excess lifetime at \code{c1} whenever \code{xcal} precedes the latter.
#' @return list with new birthdates (\code{xcal}), excess lifetime at \code{c1} (\code{slow}),
#' excess lifetime above \code{u} (\code{dat}) and right-censoring indicator (\code{rightcens}).
sample_lifetime <- function(n, pars, xcal, c1, c2, slow){
  
  sample_dates <- function(n, xcal, c1, c2, slow){
    sample_slow <- function(n, slow){
      sort(round(rexp(n, rate = 1/mean(365.25*slow[slow>0])))/365.25, decreasing = TRUE)
    }
    nday <- as.numeric(xcal[ind]-c1)
    nday <- nday[nday>0]
    nmax <- as.numeric(c2-c1)
    sslow <- sample_slow(round(sum(slow>0)*n/length(slow)), slow = slow)
    xhist <- hist(nday, plot = FALSE)
    bins <- with(xhist, sample(length(mids), n-length(sslow), p=density, replace=TRUE)) # choose a bin
    result <- round(runif(length(bins), xhist$breaks[bins], pmin(xhist$breaks[bins+1], nmax-1)))
    list(xcal = as.Date(round(c(-sslow*365.25, sort(result))), origin = "2009-01-01"),
         slow = c(sslow, rep(0, n-length(sslow))))
  }
  
  traject <- evd::rgpd(n = n, loc = 0, scale = pars[1], shape = pars[2])
  sdates <- sample_dates(n = n, xcal, c1, c2, slow)
  lifeah <- pmax(1,pmin(round(365.25*(traject + sdates$slow)), as.numeric(c2-sdates$xcal)))
  rcens_new <- lifeah == as.numeric(c2-sdates$xcal)
  list(xcal = sdates$xcal, slow = sdates$slow, dat = lifeah/365.25, rightcens = rcens_new)
}

#' Likelihood root function for shape parameter of the generalized Pareto distribution
#'
#' This function returns the likelihood root of the profile log-likelihood for the shape of the generalized Pareto distribution with left-truncated and right-censored data.
#' Specifically, \eqn{-r^2/2} is the profile likelihood and the two-sided p-value is\code{qchisq(p, 1)/2}.
#'
#' @param psi value of the shape at which to compute the p-value
#' @param thetahat maximum likelihood estimates of the scale and shape parameters
#' @inheritParams gpd_cens
#' @return a p-value
rfun_xi_gpd_cens <- function(psi, thetahat, dat, rightcens, slow){
  if(abs(psi) > 1e-6){
    llp <- prof_gpd_cens_xi(xi = psi, dat = dat, rightcens = rightcens, slow = slow)
  } else{
    rate <- exp_mle_lt_rc(dat = dat, rightcens = rightcens, slow = slow)
    llp <- c(-2*gpd_cens(par = c(rate, 0), dat = dat,  rightcens = rightcens, slow = slow), 0)
  }
  -2*gpd_cens(par = thetahat, dat = dat, rightcens = rightcens, slow = slow) - llp[1]
}



#' Likelihood for left-truncated and right-censored extended generalized Pareto variates
#'
#' Computes the log-likelihood for extended generalized Pareto distribution.
#' This distribution has survival function
#' \deqn{(1+\xi\frac{\exp(\beta x/\sigma)-1}{\beta}\right)_{+}^{-1/\xi}}
#' where \eqn{\sigma} is the scale parameter, \eqn{\beta} is a shape parameter of the Gompertz distribution, and 
#' \eqn{\xi} is the shape parameter of the generalized Pareto distribution.
#' @param par vector of scale and shape parameters
#' @param dat vector of threshold exceedances
#' @param rightcens logical indicating right-censoring (\code{TRUE} for censored)
#' @param slow lower truncation limit
#' @param expo logical; should an exponential model be fitted? Default to \code{FALSE}
#' @return log-likelihood value
exggp_cens <- function(par, dat, rightcens, slow){
  #scale par[1]
  # beta = par[2]; if zero, recover generalized Pareto
  #xi =  par[3]; if zero, recover Gompertz
  maxdat <- max(dat)
  # both shape parameters can be negative, but the data
  # must lie inside the support of the distribution
  if(abs(par[3]) <= 1e-8 && abs(par[2]) <= 1e-8){
    bounds <- par[1]
  } else if(abs(par[3]) > 1e-8 && abs(par[2]) < 1e-8){
    bounds <- c(par[1], ifelse(par[3] < 0, 1+par[3]*maxdat/par[1], Inf))
  } else if(abs(par[3]) < 1e-8 && abs(par[2]) > 1e-8){
    bounds <- par[1:2]
  } else {
    bounds <- c(par[1:2], 1+par[3]*(exp(par[2]*maxdat/par[1])-1)/par[2])
  }
  if(isTRUE(any(bounds <= 0))){
    return(1e10)
  }
  g1 <- intersect(which(!rightcens), which(slow > 0))
  g2 <- intersect(which(rightcens), which(slow > 0))
  ll <- 0
  sexggp <- function(dat, par){
    if(abs(par[3]) < 1e-8 && abs(par[2]) < 1e-8){ #exponential
      exp(-dat/par[1])
    } else if(abs(par[3]) < 1e-8 && abs(par[2]) > 1e-8){ #Gompertz
      exp(-exp(par[2]*dat/par[1])/par[2] + 1/par[2])
    } else if(abs(par[3]) >= 1e-8 && abs(par[2]) < 1e-8){ #generalized Pareto
      (1+dat*par[3]/par[1])^(-1/par[3])
    } else if(abs(par[3]) >= 1e-8 && abs(par[2]) >= 1e-8){ #extended
      (1+par[3]*(exp(par[2]*dat/par[1])-1)/par[2])^(-1/par[3]) 
    }
  }
  dexggp <- function(dat, par){
    if(abs(par[3]) < 1e-8 && abs(par[2]) < 1e-8){ 
      exp(-dat/par[1])/par[1]
    } else if(abs(par[3]) < 1e-8 && abs(par[2]) > 1e-8){ #Gompertz
      exp(par[2]*dat/par[1] - exp(par[2]*dat/par[1])/par[2] + 1/par[2])/par[1]
    } else if(abs(par[3]) >= 1e-8 && abs(par[2]) < 1e-8){ #generalized Pareto
      (1+dat*par[3]/par[1])^(-1/par[3] - 1)/par[1]  
    } else if(abs(par[3]) >= 1e-8 && abs(par[2]) >= 1e-8){ #extended
      (par[3]*(exp(par[2]*dat/par[1]) - 1)/par[2] + 1)^(-1/par[3] - 1)*exp(par[2]*dat/par[1])/par[1]
    }
  }
  #Contribution from observations in sampling frame (density divided by truncation interval)
  if(sum(!rightcens)>0){
    ll <- sum(log(dexggp(par = par, dat = dat[!rightcens])))
    if(length(g1) > 0){
      ll <- ll - sum(log(sexggp(dat = slow[g1], par = par)))  #right censored individuals
    }
  }
  if(sum(rightcens)>0){
    ll <- ll + sum(log(sexggp(dat = dat[rightcens], par = par)))
    if(length(g2) > 0){
      ll <- ll - sum(log(sexggp(dat = slow[g2], par = par)))
    }
  }
  if (!is.finite(ll)) {
    return(1e10)
  } else {
    return(-ll)
  }
}

#' Nelson-Aalen estimator of the hazard for left-truncated right-censored observations
#' 
#' This estimator is an alternative to \code{survfit} from the \code{survival} package. It does not include a correction for ties.
#' 
#' @param time vector of failure or censoring time
#' @param lefttrunc vector of left truncation time
#' @param rightcens logical vector, \code{TRUE} if the observation is right-censored, \code{FALSE} otherwise
#' @param level confidence level, with default value of 0.95
#' @param prob percentiles, with default giving equisided confidence intervals
#' @return a list with arguments
#' \itemize{
#' \item{time}{unique failure times}
#' \item{n.event}{number of failures at \code{time}}
#' \item{n.risk}{size of the risk set at \code{time}}
#' \item{haz}{hazard at \code{time}}
#' \item{cumhaz}{cumulative hazard at \code{time}}
#' \item{var}{variance of the cumulative hazard at \code{time}}
#' \item{confint}{log-transformed \code{level} confidence intervals for the cumulative hazard}
#' \item{confband}{approximate log-transformed equal precision bands for the cumulative hazard}
#'}
nelson.aalen <- 
  function(time,
           lefttrunc,
           rightcens,
           tl,
           tu,
           level = 0.95,
           prob = c((1 - level)/2, 1 - (1 - level)/2)
  ){
    stopifnot(
      length(lefttrunc) == length(rightcens),
      length(lefttrunc) == length(time),
      isTRUE(all(time > lefttrunc)), 
      is.logical(rightcens),
      length(prob) == 2L,
      length(level) == 1L,
      level > 0 & level < 1,
      all(prob > 0, prob < 1))

    prob <- sort(prob)
    # unique failure times
    utimes <- unique(sort(time))
    if(missing(tl)){
      tl <- utimes[1] 
      tl_index <- 1L
    } else{
      if(!tl %in% utimes){
        stop("Lower bound for time interval must be one of the observed failure times.")
      }
      tl_index <- which(utimes == tl)
      if(!tl > 0){
        stop("Failure time must be strictly positive")
      }
    }
    if(missing(tu)){
      tu_index <- length(utimes) - 1L
      tu <- utimes[tu_index]
    } else{
      if(!tu %in% utimes){
        stop("Upper bound for time interval must be one of the observed failure times.")
        tu_index <- which(utimes == tu)
      }
    }
    stopifnot(tl < tu)
    risk <- fail <- vector(mode = "numeric",
                           length = length(utimes))
    for(i in seq_along(utimes)){
      ti <- utimes[i]
      #risk set: not failed or censored
      risk[i] <- sum(I(lefttrunc < ti)*I(ti <= time))
      #failures: how many death at ti
      fail[i] <- sum(time[!rightcens] == ti)
    }
    # Nelson-Aalen estimator
    na_est <- cumsum(fail/risk)
    na_haz <- fail/risk
    # Variance estimator (Greenwood formula)
    na_var <- cumsum((risk - fail) * fail/
                       ((risk - 1) * risk^2))
    # Log-transform pointwise confidence intervals
    na_confint <- cbind(
      na_est*exp(qnorm(prob[1])*na_var/na_est),
      na_est*exp(qnorm(prob[2])*na_var/na_est))
    # Equal precision bands
    ep_critical <- function(d){
      4*dnorm(d)/d+2*dnorm(d)*(d-1/d)*
        0.5*(log(na_var[tu_index])-log(na_var[tl_index])) - (1-level)
    }
    dalpha <- uniroot(ep_critical, 
                      interval = c(qnorm(level), 100))$root
    na_confband <- cbind(
      na_est*exp(-dalpha*na_var/na_est),
      na_est*exp(dalpha*na_var/na_est))
    return(list(
      time = utimes,
      n.event = fail,
      n.risk = risk,
      haz = na_haz,
      cumhaz = na_est,
      var = na_var,
      confint = na_confint,
      confband = na_confband
    ))
  }

#' Confidence band, adapted from code by David Diez
# cumhazbands <- function(time, 
#                         cumhaz, 
#                         type = c("ptwise", "band"), 
#                         sigma.sq, 
#                         transfo = c("none", "log", "arcsin"),
#                         confLevel = c(0.9, 0.95, 0.99)) {
#   transfo <- match.arg(transfo)
#   type <- match.arg(type)
#   t.U <- last(time)
#   t.L <- first(time)
#   n <- length(time)
#   if(type == "band"){
#     a.L <- n * sigma.sq[1]/(1 + n * sigma.sq[1])
#     a.U <- n * sigma.sq[n]/(1 + n * sigma.sq[n])
#     aU <- format(c(round(50 * a.U)/50, 0.01))[1]
#     aL <- format(c(round(50 * a.L)/50, 0.01))[1]
#     if (confLevel[1] == 0.9) {
#       test <- try(data(hw.k10))
#       if(is.matrix(test)){
#         input <- hw.k10[aU, aL]
#       } else{
#         input <- OIsurv:::hw.k10 [aU, aL]
#       }
#     }
#     else if (confLevel[1] == 0.95) {
#       test <- try(data(hw.k05))
#       if(is.matrix(test)){
#         input <- hw.k05[aU, aL]
#       } else{
#         input <- OIsurv:::hw.k05[aU, aL]
#       }
#     }
#     else if (confLevel[1] == 0.99) {
#       test <- try(data(hw.k01))
#       if(is.matrix(test)){
#         input <- hw.k01[aU, aL]
#       } else{
#         input <- OIsurv:::hw.k01 [aU, aL]
#       }
#     }  else {
#       stop("Only confidence levels of 0.90, 0.95 and 0.99\nare allowed.")
#     }
#     m.fact <- input * (1 + n * sigma.sq)/sqrt(n)
#   } else{
#     m.fact <- qnorm((1+confLevel[1])/2)*sqrt(sigma.sq)
#   }
#   CI <- matrix(NA, length(time), 2)
#   if(transfo == "none"){
#     CI[, 1] <- cumhaz - m.fact
#     CI[, 2] <- cumhaz + m.fact
#   } else if(transfo == "log"){
#     CI[, 1] <- cumhaz*exp(-m.fact*sqrt(sigma.sq)/cumhaz)
#     CI[, 2] <- cumhaz*exp(m.fact*sqrt(sigma.sq)/cumhaz)
#   } else if(transfo == "arcsin"){
#     CI[,1] <- -2*log(sin(pmin(pi/2, asin(exp(-cumhaz/2))+0.5*m.fact*sqrt(sigma.sq/(exp(cumhaz)-1)))))
#     CI[,2] <- -2*log(sin(pmax(0, asin(exp(-cumhaz/2))-0.5*m.fact*sqrt(sigma.sq/(exp(cumhaz)-1)))))                         
#   }
#   CI[CI[, 1] < 0, 1] <- 0
#   tR <- list(time = time, lower = CI[, 1], upper = CI[, 2])
#   class(tR) <- "confBands"
#   return(tR)
# }

cumhaz_gp <- function (x, 
                       time, 
                       time2 = NULL, 
                       event = NULL, 
                       status = NULL,
                       thresh = 0, 
                       ltrunc = NULL, 
                       rtrunc = NULL,
                       type = c("right", "left", "interval", "interval2"), 
                      weights = rep(1,length(time)), 
                      level = 0.95, 
                      psi = NULL, 
                      plot = FALSE) 
{
  stopifnot(level >  0 && level < 1,
            length(level) == 1L,
            is.numeric(x))
  type <- match.arg(type)
  mle <- longevity::fit_elife(time = time, 
                   time2 = time2, 
                   event = event, 
                   status = status, 
                   thresh = thresh, 
                   ltrunc = ltrunc, 
                   rtrunc = rtrunc, 
                   type = type, family = "gp",
                   weights = weights)
  gp_cumhaz <- function(par, x){
    log(1+x*par[2]/par[1])/par[2]
  }
  mle_cumhaz <- gp_cumhaz(par = mle$par, x = x)
    if (mle$par[2] < 0 && x > -mle$par[1]/mle$par[2]) {
      stop("Value of x is outside of the range of the distribution evaluated at the maximum likelihood estimate.")
    }
  # Given cumulative hazard, retrieve scale parameter
    inv_cumhaz <- function(cumhaz, xi, x) {
      as.numeric((xi*x)/(exp(xi*cumhaz)-1))
    }
    if (is.null(psi)) {
      haz_stderror <- try(sqrt(diag(solve(numDeriv::hessian(func = function(par) {
        longevity::nll_elife(
          par = c(inv_cumhaz(par, xi = mle$par[2], x = x),
                  mle$par[2]), 
          time = time, 
          time2 = time2, 
          event = event, 
          thresh = thresh, 
          type = type, 
          ltrunc = ltrunc, 
          rtrunc = rtrunc, 
          family = "gp",
          weights = weights)
      }, x = mle_cumhaz)))))
      if (is.character(haz_stderror)) {
        stop("Could not find a grid of values for the hazard confidence interval: please provide `psi` argument.")
      }
      psi <- mle_cumhaz + 
        seq(from = -4 * haz_stderror, 
            to = 4 * haz_stderror, 
            length.out = 101L)
    }
    psi <- psi[psi > 0]
    mdat <- max(time, time2, na.rm = TRUE) - thresh
    ubound <- ifelse(x > mdat, x, mdat)
    npll <- sapply(psi, function(cumhaz_i) {
      opt <- optimize(f = function(xi) {
        longevity::nll_elife(par = c(inv_cumhaz(cumhaz_i, xi = xi, x = x), xi),
                             time = time,
                             time2 = time2,
                             event = event,
                             ltrunc = ltrunc,
                             rtrunc = rtrunc, 
                             weights = weights, 
                             family = "gp",
                             type = type, 
                             thresh = thresh)
      }, upper = 1, lower = -1)
      opt$objective
    })
  profile_confint <- function(psi, npll, psi_mle, nll_mle, 
                              level = 0.95) {
    stopifnot(is.numeric(psi_mle),
              length(psi_mle) == 1L,
              min(psi) < psi_mle & max(psi) > psi_mle,
              length(psi) == length(npll),
              nll_mle <= min(npll))
    pred <- c(predict(smooth.spline(x = npll[psi < psi_mle] - 
                                      nll_mle, y = psi[psi < psi_mle]), qchisq(level, 1)/2)$y, 
              predict(smooth.spline(x = npll[psi > psi_mle] - nll_mle, 
                                    y = psi[psi > psi_mle]), qchisq(level, 1)/2)$y)
    return(pred)
  }
  mnpll <- which.min(npll)
  if (npll[mnpll] < -mle$loglik) {
    mle_haz <- psi[mnpll]
    nll_mle <- npll[mnpll]
  }  else {
    nll_mle <- -mle$loglik
  }
  pconfint <- profile_confint(psi = psi, 
                              psi_mle = mle_cumhaz,
                              npll = npll, 
                              nll_mle = nll_mle, 
                              level = level)
  shifted_npll <- -(npll + mle$loglik)
  keep <- shifted_npll > -7
  retv <- structure(list(hazards = psi[keep], 
                         pll = shifted_npll[keep], 
                         confint = pconfint, 
                         par = as.numeric(mle_cumhaz), 
                         level = level), 
                    class = "elife_cumhazard")
  if (plot) {
    plot(retv)
  }
  return(invisible(retv))
}

ep_critical <- function(d, cumhazstd, tL, tU){
  4*dnorm(d)/d+2*dnorm(d)*(d-1/d)*
    (log(cumhazstd[tU])-log(cumhazstd[tL])) - 0.05
}