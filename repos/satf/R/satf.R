.packageName <- "satf"

Deviance <- function(LL) -2*LL
AIC <- function(LL, k) Deviance(LL) + 2*k

p2logodds <- function(p) log(p/(1-p))
logodds2p <- function(lodds) exp(lodds)/(1+exp(lodds))

NumSmallest <- -.Machine$double.xmax
NumLargest <- .Machine$double.xmax

compute_satf_criterion <- function(data, dv, bias, cnames) {
  library(plyr)
  formula.iv <- sprintf("poly(%s, %d)", cnames['time'], bias$poly.degree)
  bias.cnames <- bias$bias
  if(length(bias.cnames) > 0) {
    data$bias.id <- plyr::id(data[,bias.cnames])
    formula.iv <- sprintf("bias.id:%s", formula.iv)
  }
  
  if('response' %in% names(dv)) 
  {
      is.noise <- (data[,cnames['signal']] == 0)
      formula <- sprintf("%s==0 ~ %s", dv['response'], formula.iv)
      formula <- eval(parse(text=formula))
      m <- glm(formula=formula, family=binomial(link="probit"), data=data, 
               subset=is.noise)
      predicted.crit <- predict(m, data)

  } else if( all(c('n.responses.yes','n.responses') %in% names(dv)) )
  {    
    data.crit <- ddply(data, c(bias[['bias']], cnames[['time']]), function(d) {
        is.noise <- (d[,cnames['signal']] == 0)
        d.noise <- d[is.noise,]
        n.yes = sum(d.noise[[ dv['n.responses.yes'] ]])
        n = sum(d.noise[[ dv['n.responses'] ]])
        d$crit <- qnorm(1-n.yes/n)
        d
    })
    predicted.crit <- data.crit[rownames(data), 'crit']
  }
  reportifnot(!any(is.na(predicted.crit)) && all(is.finite(predicted.crit)), "Invalid criterion predicted." )
  predicted.crit
}

                  
satf  <- function(dv, signal, start, contrasts, data, time, metric, bias, trial.id=NULL, 
                  constraints=list(), optim.control=list(), optim.digits=1, plot=FALSE,  
                  method="Nelder-Mead", ...)
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
  
  predicted.criterion <- compute_satf_criterion(data, params$dv, params$bias, params$cnames)
  
  # initialize the C++ optimization routine
  rcpp_initialize_logLikFn(params$dv, params$contrasts, params$constraints, 
                           data, predicted.criterion, params$cnames)
  
  start <- rcpp_unconstrain_coefs( params$start )
  
  if(length(start) == 0) {
    return( rcpp_compute_logLikFn(start, FALSE) )
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
