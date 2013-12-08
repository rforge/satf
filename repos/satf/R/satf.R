.packageName <- "satf"

Deviance <- function(LL) -2*LL
AIC <- function(LL, k) Deviance(LL) + 2*k

p2logodds <- function(p) log(p/(1-p))
logodds2p <- function(lodds) exp(lodds)/(1+exp(lodds))

NumSmallest <- -.Machine$double.xmax
NumLargest <- .Machine$double.xmax


                  
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

  start <- rcpp_unconstrain_coefs( coefs$start )
  names(start) <- names( coefs$start )
  
  if(all(coefs$fixed.coefs)) {
    if(likelihood.byrow) {
      res <- rcpp_compute_logLikFn(start, TRUE, TRUE)
      rcpp_deinitialize_logLikFn()
      return( res )
      
    } else {
      res <- rcpp_compute_logLikFn(start, FALSE, TRUE)
      rcpp_deinitialize_logLikFn()
      return( res )
    }
  }
  
  # make sure the start values yield a valid likelihood
  res <- rcpp_compute_logLikFn(start)
  if( is.na(res) ) {
    rcpp_deinitialize_logLikFn()
    stop("Got NA on first iteration. Adjust start parameters.")
    
  } else if( is.infinite(res) ) {
    rcpp_deinitialize_logLikFn()
    stop("Got Inf or -Inf on first iteration. Adjust start parameters.")
  }
  
  # optimization consists of several steps
  fixed.coefs <- coefs$fixed.coefs
  dprime.coefs <- names(start) %in% dm$contrasts.coefs
  bias.coefs <- names(start) %in% dm$bias.coefs
  
  # STEP 1: optimize over criterion parameters (with correlation, if specified) use only noise trials
  # fix satf parameters
  signal <- data [[ params$cnames['signal'] ]]
  rcpp_select_subset( !signal )
  res <- maxLik(logLik=rcpp_compute_logLikFn, start=start, fixed=(fixed.coefs | dprime.coefs), method="Nelder-Mead", iterlim=10^6)
  rcpp_reset_selection()
  start <- coef(res)
  
  # STEP 2: optimize over d' parameters (with correlation, if specified)
  rcpp_select_subset( signal )
  res <- maxLik(logLik=rcpp_compute_logLikFn, start=start, fixed=(fixed.coefs | !dprime.coefs), method="Nelder-Mead", iterlim=10^6)
  rcpp_reset_selection()
  start <- coef(res)
  
  # STEP 3: optimize over all parameters
  res <- maxLik(logLik=rcpp_compute_logLikFn, start=start, fixed=(fixed.coefs), method="Nelder-Mead", iterlim=10^6)

  if(likelihood.byrow) {
    logLik <- rcpp_compute_logLikFn(coef(res), TRUE, TRUE)
    rcpp_deinitialize_logLikFn()
    return( logLik )
  } 

  # recompute the likelihood for the entire dataset (this time, do not tolerate approximation errors for extreme values in the C++ code)
  logLik <- rcpp_compute_logLikFn( coef(res), FALSE, TRUE)

  # transform the parameters into their proper (constrained) form
  estimates <- rcpp_constrain_coefs( coef(res) )
  
  res <- list(estimates=estimates, LL=logLik)
  
  # free all memory reserved by C++
  rcpp_deinitialize_logLikFn()
  
  return(res)
}
