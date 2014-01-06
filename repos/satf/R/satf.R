.packageName <- "satf"

Deviance <- function(LL) -2*LL
AIC <- function(LL, k) Deviance(LL) + 2*k

p2logodds <- function(p) log(p/(1-p))
logodds2p <- function(lodds) exp(lodds)/(1+exp(lodds))

NumSmallest <- -.Machine$double.xmax
NumLargest <- .Machine$double.xmax


compute_logLikFn <- function(coefs, by_row=FALSE, tolerate_imprecision=TRUE, force_update=FALSE) {
	res = rcpp_compute_logLikFn(coefs, by_row, tolerate_imprecision, force_update)
  res
}

compute_logLikFn_gradient <- function(coefs, by_row=FALSE, tolerate_imprecision=TRUE) {
  res = rcpp_compute_logLikFn_gradient(coefs, by_row, tolerate_imprecision)
  names(res) = names(coefs)
  res
}

parameter_combinations <- function(namedlist) {
  membernames <- names(namedlist)
  cross_prod = cset(namedlist[[1]])
  for(i in 2:length(namedlist))
    cross_prod = cross_prod*cset(namedlist[[i]])
  llply(as.list(cross_prod), function(x) { x= unlist(x); names(x) = membernames; x })
}

satf_gridsearch <- function(start, constraints, ...) {
  rcpp_deinitialize_logLikFn()
  start.vec <- parameter_combinations(start)
  logLik.vec <- laply(start.vec, function(start) {
    constraints[names(start)] = start
    res = satf(start=start, constraints=c(constraints,start), ..., stepwise=F, .internal.init.optional=TRUE, .internal.cleanup=FALSE)
    reportifnot(!is.list(res), "All coefficient values have to be specified in a grid search.")
    res
  })
  rcpp_deinitialize_logLikFn()
  start.vec[[which.max(logLik.vec)]]
}

satf  <- function(dv, signal, start, contrasts, bias, data, time, metric, trial.id=NULL, 
                  constraints=list(),
                  reoptimize.criterion=FALSE, reoptimize.corr=FALSE, 
                  likelihood.byrow=FALSE, debug=F, method="Nelder-Mead",
                  .internal.init.optional=FALSE, .internal.cleanup=TRUE)
{
  metric.permissible <- c('RMSD','R2','adjR2','logLik','logLikRaw')
  reportifnot(metric %in% metric.permissible, sprintf("'metric' has to be one of: %s", paste(metric.permissible, collapse=", ")))

  deinitialize_logLikFn <- function() T
  if(.internal.cleanup)
    deinitialize_logLikFn <- rcpp_deinitialize_logLikFn 
  
  # Check 'data' parameter
  reportifnot(nrow(data) > 0, "Parameter 'data' needs to consist of several rows.")
  
  satf.coefnames.core <- c('asymptote', 'invrate', 'intercept')
  bias.coefnames.core <- c('bias.max', 'bias.invrate', 'bias.intercept','bias.min')
  
  # set defaults for start values and constraints
  start = set_start_defaults(start, set.corr.mrsat=!is.null(trial.id))
  constraints = set_constraints_defaults(constraints)
  
  # check parameters
  params <- translate.parameters(data=data, dv=dv, contrasts=contrasts, bias=bias,
                                 signal=signal, time=time, trial.id=trial.id)

  skip.initialization <- .internal.init.optional && rcpp_is_initialized_logLikFn()
  
  # initialize start parameters and create constraint matrix if necessary
  if(!skip.initialization) {
    dm <- init_designmatrix(data=data, contrasts=params$contrasts, bias=params$bias, cnames=params$cnames,
                            satf.coefnames.core=satf.coefnames.core, bias.coefnames.core=bias.coefnames.core)
    coefnames <- colnames(dm$dm)
  } else {
    coefnames <- rcpp_get_coef_names()
  }

  # initialize start parameters and create constraint matrix  
  coefs <- init_coefs_and_constraints(coefnames=coefnames, start=start, constraints=constraints,
                                      coreparams=c(satf.coefnames.core, bias.coefnames.core))
  # initialize the C++ optimization routine
  if( !skip.initialization ) {
    rcpp_initialize_logLikFn(params$dv, dm$dm, dm$dm.coef.cnt, coefs$constraints, data, params$cnames)
  } else {
    rcpp_update_constraints_logLikFn(coefs$constraints)
  }

  start <- rcpp_unconstrain_coefs( coefs$start )
  names(start) <- names( coefs$start )

  if(all(names(start) %in% coefs$fixed.coefs)) {
    if(likelihood.byrow) {
      logLik <- compute_logLikFn(start, TRUE, TRUE)
      deinitialize_logLikFn()
      return( logLik )
      
    } else {
      logLik <- compute_logLikFn(start, FALSE, TRUE)
      deinitialize_logLikFn()
      return( logLik )
    }
  }

  # make sure the start values yield a valid likelihood
  logLik = compute_logLikFn(start)
  if( is.na(logLik) ) {
    if(likelihood.byrow)
     	res = compute_logLikFn(start, TRUE, FALSE)
    deinitialize_logLikFn()
    warning("Got NA on first iteration. Adjust start parameters.")
    return(logLik)
    
  } else if( is.infinite(logLik) ) {
    if(likelihood.byrow)
      logLik = compute_logLikFn(start, TRUE, FALSE)
    deinitialize_logLikFn()
    warning("Got Inf or -Inf on first iteration. Adjust start parameters.")
    return(logLik)
  }

  optimize_subset <- function(variable, start, all=TRUE, selection=c()) {
    .optimize_subset(variable, coefs$fixed.coefs, method, start, debug, data, selection)
  }
  
  # TODO: Implement the following.
  reportifnot(!.internal.init.optional, "The package does not support stepwise=T with .internal.init.optional=T. It's on my to-do list.")

  log_step_n = function(n) {
    if(debug) {
      print("-------------------")
      print(sprintf("optimization step %d", n))
      print("-------------------")
    }
  }

  # ignore correlation parameter for now, if specified
  rcpp_set_coef_values( c(corr.mrsat=0) )

  coeforder = append(dm$params.criterion, dm$params.dprime)

  # optimize in steps on subsets of the data
  for(i in 1:length(coeforder)) {
    log_step_n(i)
    selection.variables = c(unlist(coeforder[1:i]),  'corr.mrsat' )
    start = optimize_subset(variable=coeforder[[i]], start=start, selection=selection.variables) 
  }

  # optimize over the correlation coefficient
  log_step_n(i+1)
  rcpp_reset_coef_ranges( 'corr.mrsat' )
  start = optimize_subset(variable='corr.mrsat', start=start) 

  # reoptimize dprime or more parameters
  log_step_n(i+2)
  free.variables = unlist(dm$params.dprime)
  if(reoptimize.criterion) free.variables = c(free.variables, unlist(dm$params.criterion))
  if(reoptimize.corr)      free.variables = c(free.variables, 'corr.mrsat') 
  estimates = optimize_subset(variable=free.variables, start=start) 
  if(likelihood.byrow) {
    logLik <- compute_logLikFn(coefs=estimates, by_row=TRUE, tolerate_imprecision=TRUE)
    deinitialize_logLikFn()
    return( logLik )
  } 

  # recompute the likelihood for the entire dataset (this time, do not tolerate approximation errors for extreme values in the C++ code)
  logLik <- compute_logLikFn( coefs= estimates, by_row=FALSE, tolerate_imprecision=TRUE)
  res = list(estimates=rcpp_constrain_coefs(estimates), LL=logLik)
  deinitialize_logLikFn()
  return(res)
}

.optimize_subset <- function(variable, fixed, method, start, debug, data, selection)
{
  # NOTE: 'fixed' overrides 'variable'
  # extract start values if necessary
  original.start = start
  if(!is.vector(start))
    start = coef(start)
  
  # determine which variables should vary, and which are fixed
  variable.coefnames = setdiff(variable, fixed)
  fixed.coefnames = setdiff(names(start), variable.coefnames)
  fixed = names(start) %in% fixed.coefnames
  
  # if all are fixed, return the original coefs
  if( all(fixed) )
    return(original.start)
  
  # set gradient, etc.
  # transform the parameters into their proper (constrained) form
  fnLogLik = compute_logLikFn
  fnLogLikGradient = NULL
  
##  start.constrained <- rcpp_constrain_coefs( start )
##  corr.mrsat.nonzero = ( 'corr.mrsat' %in% names(start) ) && !( 'corr.mrsat' %in% fixed.coefnames && start.constrained['corr.mrsat'] == 0 )
##  if(corr.mrsat.nonzero) {
##    # print('corr.mrsat.nonzero is TRUE')
##    method = "Nelder-Mead"
##  } else {
##    # print('corr.mrsat.nonzero is FALSE')
##    fnLogLikGradient = compute_logLikFn_gradient
##    if( is.null(method) )
##      method = "BFGS"
##  }
  
  if(debug) print.level=0
  else      print.level=0
  
  # select the required data subset
  if(length(selection) > 0) {
    rcpp_select_coef_subset( selection )
#    print("--")
#    print(selection)
#    print(rcpp_return_selection( ))
#    print(data[rcpp_return_selection( ),])
#    print("--")
  }

  n.free = sum(!fixed)
  if(n.free == 1) {
    res <- start
    fnLogLikTmp = function(p) {res[!fixed] = p; ll = fnLogLik(res); ifelse(ll == -Inf, -.Machine$double.xmax, ll);  }
    res.optim = optimize(fnLogLikTmp, lower=-6, upper=6, maximum=TRUE)
    res[!fixed] = res.optim$maximum
    
  } else if(method == "Nelder-Mead") {
    if(any(start==0.0)) parscale = rep(1, length(start))
    else                parscale = abs(start)
    res = maxLik(logLik=fnLogLik, grad=compute_logLikFn_gradient, start=start, fixed=fixed,
                 parscale=parscale, iterlim=10^6, method="Nelder-Mead", print.level=print.level)
    if(debug) {
      variable.coefnames <- variable.coefnames[variable.coefnames%in%names(start)]
      print(sprintf("optimizing: %s", paste(variable.coefnames, collapse=', ') ))
      print(sprintf("method: %s", method))
      print(sprintf("code: %d", res$code))
      print(sprintf("iterations: %d", res$iterations))    
      print("coefs")
      print( rcpp_constrain_coefs( coef(res) ) )
      old.LL = compute_logLikFn( coefs=start, by_row=FALSE, tolerate_imprecision=TRUE)
      new.LL = compute_logLikFn( coefs=coef(res), by_row=FALSE, tolerate_imprecision=TRUE)
      print(sprintf("LL improved by %.2f (old LL = %.2f, new LL = %.2f)", new.LL-old.LL, new.LL, old.LL))
      print("---")
    }
    res = coef(res)
    
  } else {
    stop("This optimization method is not implemented.")
    ## res = maxLik(logLik=fnLogLik, grad=compute_logLikFn_gradient, start=start, fixed=fixed, iterlim=10^6, 
    ##             method="Nelder-Mead", print.level=print.level)
  }
  rcpp_reset_selection()
  res
}
