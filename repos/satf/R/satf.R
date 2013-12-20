.packageName <- "satf"

Deviance <- function(LL) -2*LL
AIC <- function(LL, k) Deviance(LL) + 2*k

p2logodds <- function(p) log(p/(1-p))
logodds2p <- function(lodds) exp(lodds)/(1+exp(lodds))

NumSmallest <- -.Machine$double.xmax
NumLargest <- .Machine$double.xmax


compute_logLikFn <- function(coefs, by_row=FALSE, tolerate_imprecision=TRUE) {
	res = rcpp_compute_logLikFn(coefs, by_row, tolerate_imprecision)
#print(res)
  res
}

compute_logLikFn_gradient <- function(coefs, by_row=FALSE, tolerate_imprecision=TRUE) {
  res = rcpp_compute_logLikFn_gradient(coefs, by_row, tolerate_imprecision)
  names(res) = names(coefs)
#  print(res)
  res
}

satf  <- function(dv, signal, start, contrasts, bias, data, time, metric, trial.id=NULL, 
                  constraints=list(), optim.control=list(), optim.digits=1, plot=FALSE,  
                  likelihood.byrow=FALSE, method="Nelder-Mead", stepwise=TRUE, debug=F, ...)
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
  rcpp_initialize_logLikFn(params$dv, dm$dm, dm$dm.coef.cnt, coefs$constraints, data, params$cnames)
  start <- rcpp_unconstrain_coefs( coefs$start )
  names(start) <- names( coefs$start )

  if(all(names(start) %in% coefs$fixed.coefs)) {
    if(likelihood.byrow) {
      res <- compute_logLikFn(start, TRUE, TRUE)
      rcpp_deinitialize_logLikFn()
      return( res )
      
    } else {
      res <- compute_logLikFn(start, FALSE, TRUE)
      rcpp_deinitialize_logLikFn()
      return( res )
    }
  }
  
  # make sure the start values yield a valid likelihood
  res <- compute_logLikFn(start)
  if( is.na(res) ) {
    if(likelihood.byrow) {
   	res <- compute_logLikFn(start, TRUE, FALSE)
    }
    rcpp_deinitialize_logLikFn()
    warning("Got NA on first iteration. Adjust start parameters.")
    return(res)
    
  } else if( is.infinite(res) ) {
    if(likelihood.byrow) {
   	res <- compute_logLikFn(start, TRUE, FALSE)
    }
    rcpp_deinitialize_logLikFn()
    warning("Got Inf or -Inf on first iteration. Adjust start parameters.")
    return(res)
  }

  select_zero <- function(value, coefnames){
    selection <- rep(value, length(coefnames))
    names(selection) <- coefnames
    selection
  }

  #print( compareDerivatives(f=rcpp_compute_logLikFn, grad=compute_logLikFn_gradient, t0=start) )
  #stop()

  optimize_subset <- function(selection, variable, method, start)
  {
    # select the required data subset
    if(length(selection) > 0) {
      rcpp_select_subset_by_zero_dm_columns_all( selection )
    }
    
    # extract start values if necessary
    original.start = start
    if(!is.vector(start))
      start = coef(start)
    
    # determine which variables should vary, and which are fixed
    variable.coefnames = setdiff(variable, coefs$fixed.coefs)
    fixed.coefnames = setdiff(names(start), variable.coefnames)
    fixed = names(start) %in% fixed.coefnames

    # if all are fixed, return the original coefs
    if( all(fixed) )
      return(original.start)
    
    # set gradient, etc.
    # transform the parameters into their proper (constrained) form
    fnLogLik = compute_logLikFn
    fnLogLikGradient = NULL
    
    start.constrained <- rcpp_constrain_coefs( start )
    corr.mrsat.nonzero = ( 'corr.mrsat' %in% names(start) ) && !( 'corr.mrsat' %in% fixed.coefnames && start.constrained['corr.mrsat'] == 0 )
    if(corr.mrsat.nonzero) {
      # print('corr.mrsat.nonzero is TRUE')
      method = "Nelder-Mead"
      
    } else {
      # print('corr.mrsat.nonzero is FALSE')
      fnLogLikGradient = compute_logLikFn_gradient
      if( is.null(method) )
        method = "BFGS"
    }
    
    res = maxLik(logLik=fnLogLik, grad=compute_logLikFn_gradient, start=start, fixed=fixed, iterlim=10^6, method=method)
    if(debug) {
      print(summary(res))
      old.LL = compute_logLikFn( coefs=start, by_row=FALSE, tolerate_imprecision=TRUE)
      new.LL = compute_logLikFn( coefs=coef(res), by_row=FALSE, tolerate_imprecision=TRUE)
      print(sprintf("LL improved to %.2f (old LL = %.2f)", new.LL, old.LL))
      print("---")
    }
    
    rcpp_reset_selection()
    res
  }
  
  satf.coefnames.all = dm$contrasts.coefs
  bias.coefnames.all = dm$bias.coefs
  satf.coefnames.noncore = setdiff(satf.coefnames.all, satf.coefnames.core)
  bias.coefnames.noncore = setdiff(bias.coefnames.all, bias.coefnames.core)

  log_step_n = function(n) {
    if(debug) {
      print("-------------------")
      print(sprintf("optimization step %d", n))
      print("-------------------")
    }
  }
  
  # optimizate parameters by subsets
  if(stepwise)
  {
    log_step_n(1)
    # STEP A1: optimize over core criterion parameters, use noise trials for which all specified contrasts are zero
    selection = select_zero(TRUE, union(satf.coefnames.all, bias.coefnames.noncore))
    free.variables = c(bias.coefnames.core, 'corr.mrsat')
    res = optimize_subset(selection=selection, variable=free.variables, method=method, start=start) 

    log_step_n(2)
    # STEP A2: optimize over non-core criterion parameters, use noise trials for which all specified contrasts are non-zero
    selection = c(select_zero(TRUE, satf.coefnames.all), select_zero(FALSE, bias.coefnames.noncore))
    free.variables = bias.coefnames.noncore
    res = optimize_subset(selection=selection, variable=free.variables, method=method, start=res)
  
    log_step_n(3)
    # STEP B1: optimize over core d' parameters, use signal trials for which all specified contrasts are zero
    selection = c( select_zero(TRUE, satf.coefnames.noncore) ) # c(select_zero(FALSE, satf.coefnames.core), select_zero(TRUE, satf.coefnames.noncore))
    free.variables = satf.coefnames.core
    res = optimize_subset(selection=selection, variable=free.variables, method=method, start=res)

    log_step_n(4)
    # STEP B2: optimize over non-core d' parameters, use signal trials for which all specified contrasts are non-zero
    selection = c() #c(select_zero(FALSE, satf.coefnames.core), select_zero(FALSE, satf.coefnames.noncore))
    free.variables = satf.coefnames.noncore
    res = optimize_subset(selection=selection, variable=free.variables, method=method, start=res) 
      
    log_step_n(5)
    # STEP C1: optimize over all parameters, use all data
    free.variables = c(satf.coefnames.all, bias.coefnames.all, 'corr.mrsat')
    res = optimize_subset(selection=c(), variable=free.variables, method=method, start=res)
    
  } else {
    # optimize over all parameters, use all data
    free.variables = c(satf.coefnames.all, bias.coefnames.all, 'corr.mrsat')
    res = optimize_subset(selection=c(), variable=free.variables, start=start, method=method) 
    
  }
  
  if(likelihood.byrow) {
    logLik <- compute_logLikFn(coefs=coef(res), by_row=TRUE, tolerate_imprecision=TRUE)
    rcpp_deinitialize_logLikFn()
    return( logLik )
  } 

  # recompute the likelihood for the entire dataset (this time, do not tolerate approximation errors for extreme values in the C++ code)
  logLik <- compute_logLikFn( coefs=coef(res), by_row=FALSE, tolerate_imprecision=TRUE)

  # transform the parameters into their proper (constrained) form
  estimates <- rcpp_constrain_coefs( coef(res) )

  res <- list(estimates=estimates, LL=logLik)
  
  # free all memory reserved by C++
  rcpp_deinitialize_logLikFn()
  
  return(res)
}
