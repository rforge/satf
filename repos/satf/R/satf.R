.packageName <- "satf"

Deviance <- function(LL) -2*LL
AIC <- function(LL, k) Deviance(LL) + 2*k

p2logodds <- function(p) log(p/(1-p))
logodds2p <- function(lodds) exp(lodds)/(1+exp(lodds))

NumSmallest <- -.Machine$double.xmax
NumLargest <- .Machine$double.xmax

compute_satf_criterion <- function(data, dv, bias, cnames, optim.digits=optim.digits, control=optim.control) {
  bias.cnames <- bias$bias
  
  if('response' %in% names(dv)) 
  {
      eval.fn.bias <- function(p, time, response, trial.id=as.logical(c())) {
          par <- c(lower=p[['lower']], upper=p[['upper']], center=p[['center']], stretch=p[['stretch']])
          -rcpp_compute_criterion_logLikFn(par, time, response, trial.id ) 
      }
      
      if(length(bias.cnames) == 0) {
          is.noise <- (data[[ cnames['signal'] ]] == 0)
          data.noise <- subset(data, is.noise)
          time <- data.noise[[ cnames['time'] ]]
          response <- data.noise[[ dv['response'] ]]
          res <- optim(c(lower=-1, upper=1, center=0, stretch=1), eval.fn.bias, time=time, response=response)
          predicted.crit <- rcpp_compute_criterion(res$par, data[[ cnames['time'] ]])
          
      } else {
          # TODO: Test this part of the code
          data$index <- 1:nrow(data)
          data.crit <- ddply(data, bias.cnames, function(d) {
              is.noise <- (d[[ cnames['signal'] ]] == 0)
              time <- (d[[ cnames['time'] ]])[is.noise]
              response <- (d[[ dv['response'] ]])[is.noise]
              res <- optim(c(lower=-1, upper=1, center=0, stretch=1), eval.fn.bias, time=time, response=response)
              d$crit <- rcpp_compute_criterion(res$par, time)
              d
          })
          predicted.crit <- data.crit[ordered(data.crit$index), 'crit']
      }
      
  } else if( all(c('n.responses.yes','n.responses') %in% names(dv)) )
  {
    data$index <- 1:nrow(data)
    data.crit <- ddply(data, c(bias.cnames, cnames[['time']]), function(d) {
        is.noise <- (d[,cnames['signal']] == 0)
        d.noise <- d[is.noise,]
        
        # compute criterion with slightly adjusted numbers (to avoid infinite criteria)
        # TODO: Figure out if this is really appropriate, since the addition
        #       of a constant disregards sample size.
        n.yes = 0.25 + sum(d.noise[[ dv['n.responses.yes'] ]])
        n = 0.5 + sum(d.noise[[ dv['n.responses'] ]])
        d$crit <- qnorm(1-n.yes/n )
        d
    })
    predicted.crit <- data.crit[ordered(data.crit$index), 'crit']
  }
  reportifnot(!any(is.na(predicted.crit)) && all(is.finite(predicted.crit)), "Invalid criterion predicted." )
  predicted.crit
}

                  
satf  <- function(dv, signal, start, contrasts, bias, data, time, metric, trial.id=NULL, 
                  constraints=list(), optim.control=list(), optim.digits=1, plot=FALSE,  
                  method="Nelder-Mead", likelihood.byrow=FALSE, ...)
{
  metric.permissible <- c('RMSD','R2','adjR2','logLik','logLikRaw')
  reportifnot(metric %in% metric.permissible, sprintf("'metric' has to be one of: %s", paste(metric.permissible, collapse=", ")))

  # Check 'data' parameter
  reportifnot(nrow(data) > 0, "Parameter 'data' needs to consist of several rows.")

  
  coreparams.satf <- c('asymptote', 'invrate', 'intercept')
  coreparams.bias <- c('bias.min', 'bias.max', 'bias.invrate', 'bias.intercept')
  
  # init default parameters for optimization
  optim.control$fnscale <- -1
  optim.control <- default(optim.control, 'maxit', 10^6)
  
  # set defaults for start values and constraints
  start = set_start_defaults(start, set.corr.mrsat=!is.null(trial.id))
  constraints = set_constraints_defaults(constraints)
  
  # check parameters
  params <- translate.parameters(data=data, dv=dv, contrasts=contrasts, bias=bias,
                                 signal=signal, time=time, trial.id=trial.id)
  
  # initialize start parameters and create constraint matrix  
  dm <- init_designmatrix(data=data, contrasts=params$contrasts, bias=params$bias, cnames=params$cnames,
                          coreparams.satf=coreparams.satf, coreparams.bias=coreparams.bias)
  
  # initialize start parameters and create constraint matrix  
  coefs <- init_coefs_and_constraints(coefnames=colnames(dm$dm), start=start, constraints=constraints,
                                      coreparams=c(coreparams.satf, coreparams.bias))
  
  # initialize the C++ optimization routine
  rcpp_initialize_logLikFn(params$dv, dm$dm, dm$dm.coef.cnt, coefs$constraints, 
                           data, params$cnames)

  start.tmp <- rcpp_unconstrain_coefs( coefs$start )
  
  if(length(start) == 0) {
    res <- rcpp_compute_logLikFn(start.tmp, FALSE)
    rcpp_deinitialize_logLikFn()
    return( res )
  }
  
  if(likelihood.byrow) {
    stopifnot(length(start) == 0)
    res <- rcpp_compute_logLikFn(start.tmp, TRUE)
    rcpp_deinitialize_logLikFn()
    return( res )
  }
  
  # make sure the start values yield a valid likelihood
  res <- rcpp_compute_logLikFn(start.tmp)  
  if(res == -Inf || res == Inf) {
    rcpp_deinitialize_logLikFn()
    stop("Got Inf or -Inf on first iteration.")
  }
  
  # optimization consists of several steps
  signal <- data [[ params$cnames['signal'] ]]
  start <- coefs$start
  
  optimize.subset <- function(fixed.coefs, constrained.start) 
  {
      rcpp_set_coef_values( fixed.coefs )
      
      # get unconstrained start parameters
      start <- rcpp_unconstrain_coefs( start )
      res <- optim.to.precision(start=start, optim.digits=optim.digits, control=optim.control, 
                                method=method, fn=rcpp_compute_logLikFn)
      
      res$par <- rcpp_constrain_coefs( res$par )
      rcpp_reset_coef_values( names(fixed.coefs) )
      res
  }

  # STEP 1: optimize over criterion parameters (with correlation, if specified)
  #         use only noise trials
  
  # fix satf parameters
  fixed.coef.names <- dm$contrasts.coefs
  fixed.coefs <- rep(0, length(fixed.coef.names))
  names(fixed.coefs) <- fixed.coef.names
  obtained.coefs.names <- (names(start))[!(names(start) %in% fixed.coef.names)]

  # optimize bias paramters (and correlation, if specified)
  rcpp_select_subset( !signal )
  res <- optimize.subset(fixed.coefs=fixed.coefs, constrained.start=start) 
  rcpp_reset_selection()
  
  start[ obtained.coefs.names ] <- res$par[ obtained.coefs.names ] 

  
  # step 2: optimize over d' parameters (with correlation, if specified)
  fixed.coef.names <- dm$bias.coefs
  if('corr.mrsat' %in% names(start)) {
    fixed.coef.names <- c(fixed.coef.names, 'corr.mrsat')
  }
  fixed.coefs <- start[fixed.coef.names]
  names(fixed.coefs) <- fixed.coef.names
  
  rcpp_select_subset( signal )
  res.step2 <- optimize.subset(fixed.coefs=fixed.coefs, constrained.start=start) 
  rcpp_reset_selection()

  res <- res.step2

  
##  # step 3: fit the parameters by maximizing the likelihood
##  res <- optim.to.precision(start=start, optim.digits=optim.digits, control=optim.control, 
##                            method=method, fn=rcpp_compute_logLikFn)
##res$par <- rcpp_constrain_coefs( res$par )

  # transform the parameters into their proper (constrained) form
  estimates <- res$par
  
  # recompute the likelihood for the entire dataset
  logLik <- rcpp_compute_logLikFn( rcpp_unconstrain_coefs(estimates), FALSE, FALSE)

  res <- list(estimates=estimates, LL=logLik)
  
  # free all memory reserved by C++
  rcpp_deinitialize_logLikFn()
  
  return(res)
}
