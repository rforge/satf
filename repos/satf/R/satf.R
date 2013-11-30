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

                  
satf  <- function(dv, signal, start, contrasts, data, time, metric, bias, trial.id=NULL, 
                  constraints=list(), optim.control=list(), optim.digits=1, plot=FALSE,  
                  method="Nelder-Mead", likelihood.byrow=FALSE, ...)
{
  metric.permissible <- c('RMSD','R2','adjR2','logLik','logLikRaw')
  reportifnot(metric %in% metric.permissible, sprintf("'metric' has to be one of: %s", paste(metric.permissible, collapse=", ")))

  # Check 'data' parameter
  reportifnot(nrow(data) > 0, "Parameter 'data' needs to consist of several rows.")

  # check parameters
  params <- translate.parameters(data=data, dv=dv, start=start, contrasts=contrasts, 
                                 constraints=constraints, bias=bias, signal=signal, 
                                 time=time, trial.id=trial.id, summarize=summarize)
  
  # init default parameters
  optim.control$fnscale <- -1
  optim.control <- default(optim.control, 'maxit', 10^6)

##  control <- default(control, 'correction.fraction', .01)
##  control <- default(control, 'plot.dprime', plot)
##  control <- default(control, 'plot.criterion', FALSE)
##  control <- default(control, 'plot.action', 'plot')
##  control <- default(control, 'condition', ~condition)
##  control <- default(control, 'interval', ~interval)
##  control <- default(control, 'time', ~time)

##  # init default control parameters
##  control$condition <- formula.terms(control$condition)
##  control$time <- formula.terms(control$time)
  
  predicted.criterion <- compute_satf_criterion(data, params$dv, params$bias, params$cnames,
                                                optim.digits=optim.digits, control=optim.control)
  
  # initialize the C++ optimization routine
  rcpp_initialize_logLikFn(params$dv, params$contrasts, params$constraints, 
                           data, predicted.criterion, params$cnames)
  
  start <- rcpp_unconstrain_coefs( params$start )
  
  if(length(start) == 0) {
    res <- rcpp_compute_logLikFn(start, FALSE)
    rcpp_deinitialize_logLikFn()
    return( res )
  }
  
  if(likelihood.byrow) {
    res <- rcpp_compute_logLikFn(start, TRUE)
    rcpp_deinitialize_logLikFn()
    return( res )
  }
  
  res <- rcpp_compute_logLikFn(start)
  
  if(res == -Inf || res == Inf) {
    rcpp_deinitialize_logLikFn()
    stop("Got Inf or -Inf on first iteration.")
  }

  # fit the parameters by maximizing the likelihood
  res <- optim.to.precision(start=start, optim.digits=optim.digits, control=optim.control, 
                            method=method, fn=rcpp_compute_logLikFn)
  
  # transform the parameters into their proper (constrained) form
  res$par <- rcpp_constrain_coefs( res$par )
  
  # free all memory reserved by C++
  rcpp_deinitialize_logLikFn()
  
  return(res)
}
