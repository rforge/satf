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
  
  satf.coefnames.core <- c('asymptote', 'invrate', 'intercept')
  bias.coefnames.core <- c('bias.min', 'bias.max', 'bias.invrate', 'bias.intercept')
  
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
                          satf.coefnames.core=satf.coefnames.core, bias.coefnames.core=bias.coefnames.core)
  
  # initialize start parameters and create constraint matrix  
  coefs <- init_coefs_and_constraints(coefnames=colnames(dm$dm), start=start, constraints=constraints,
                                      coreparams=c(satf.coefnames.core, bias.coefnames.core))
  
  # initialize the C++ optimization routine
  rcpp_initialize_logLikFn(params$dv, dm$dm, dm$dm.coef.cnt, coefs$constraints, 
                           data, params$cnames)

  start <- rcpp_unconstrain_coefs( coefs$start )
  names(start) <- names( coefs$start )
  
  if(all(names(start) %in% coefs$fixed.coefs)) {
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
  
  create_selection <- function(value, coefnames){
    selection <- rep(F, length(coefnames))
    names(selection) <- coefnames
    selection
  }

  optimize_subset <- function(selection, variable, start) {
    if(length(selection) > 0) {
      rcpp_select_subset_by_zero_dm_columns_all( selection )
    }
    variable.coefnames = setdiff(variable, coefs$fixed.coefs)
    fixed.coefnames = setdiff(names(start), variable.coefnames)
    fixed = names(start) %in% fixed.coefnames
    res = maxLik(logLik=rcpp_compute_logLikFn, start=start, fixed=fixed, method="Nelder-Mead", iterlim=10^6)
    rcpp_reset_selection()
    res
  }
  
  satf.coefnames.all = dm$contrasts.coefs
  bias.coefnames.all = dm$bias.coefs
  satf.coefnames.noncore = setdiff(satf.coefnames.all, satf.coefnames.core)
  bias.coefnames.noncore = setdiff(bias.coefnames.all, bias.coefnames.core)
  
  # optimization consists of several steps
print(sprintf("step A1: %s", date()))  
  # STEP A1: optimize over core criterion parameters, use only noise trials
  selection = create_selection(TRUE, union(satf.coefnames.all, bias.coefnames.noncore))
  res = optimize_subset(selection=selection, variable=bias.coefnames.core, start=start) 
print(coef(res))
  
print(sprintf("step A2: %s", date()))  
  # STEP A2: optimize over non-core criterion parameters, still using only noise trials
  selection = c(create_selection(TRUE, satf.coefnames.all), create_selection(FALSE, bias.coefnames.noncore))
  res = optimize_subset(selection=selection, variable=bias.coefnames.noncore, start=coef(res)) 
print(coef(res))
  
print(sprintf("step A3: %s", date()))  
  # STEP A3: optimize over all criterion parameters, still using only noise trials
  selection = create_selection(TRUE, satf.coefnames.all)
  res = optimize_subset(selection=selection, variable=bias.coefnames.all, start=coef(res)) 
print(coef(res))
stop()
      
print(sprintf("step A3: %s", date()))  
  # STEP A3: optimize over non-core criterion parameters, still using only noise trials
  res = optimize_subset(selection=selection, variable=dm$bias.coefs, start=coef(res)) 
print(summary(res))
  
  
  satf.coefnames.core <- c('asymptote', 'invrate', 'intercept')
  bias.coefnames.core <- c('bias.min', 'bias.max', 'bias.invrate', 'bias.intercept')
  
print(sprintf("step 2: %s", date()))  
  # STEP 2: optimize over d' parameters (with correlation, if specified)
  selection = create_selection(FALSE, satf.coefnames.core)
  res = optimize_subset(selection=selection, variable=c(dm$bias.coefs, 'corr.mrsat'), start=coef(res))
print(summary(res))
  
print(sprintf("step 3: %s", date()))  
  # STEP 3: optimize over all parameters
  res <- maxLik(logLik=rcpp_compute_logLikFn, start=coef(res), fixed=(fixed.coefs), method="Nelder-Mead", iterlim=10^6)
print(sprintf("/step 3: %s", date()))  
  
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
