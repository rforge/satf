.packageName <- "satf"


library(plyr)

Deviance <- function(LL) -2*LL
AIC <- function(LL, k) Deviance(LL) + 2*k

p2logodds <- function(p) log(p/(1-p))
logodds2p <- function(lodds) exp(lodds)/(1+exp(lodds))

NumSmallest <- -.Machine$double.xmax
NumLargest <- .Machine$double.xmax




library(reshape)

satf.criterion <- function(data, dv, bias, cnames) {
  library(plyr)
  formula.iv <- sprintf("poly(%s, %d)", cnames['time'], bias$poly.degree)
  bias.cnames <- bias$bias
  if(length(bias.cnames) > 0) {
    data$bias.id <- plyr::id(data[,bias.cnames])
    formula.iv <- sprintf("bias.id:%s", formula.iv)
  }
  
  is.noise <- (data[,cnames['signal']] == 0)
  if('response' %in% names(dv)) 
  {
    formula <- sprintf("%s==0 ~ %s", dv['response'], formula.iv)
    formula <- eval(parse(text=formula))
    m <- glm(formula=formula, family=binomial(link="probit"), data=data, 
             subset=is.noise)
    predicted.crit <- predict(m, data)

  } else if( all(c('n.responses.yes','n.responses') %in% names(dv)) )
  {
    data.crit <- ddply(data, c(bias[['bias']], cnames[['time']]), function(d) {
      d.noise <- subset(d, signal==0)
      n.yes = sum(d.noise[[ dv['n.responses.yes'] ]])
      n = sum(d.noise[[ dv['n.responses'] ]])
      d$crit <- qnorm(1-n.yes/n)
      d
    })
    predicted.crit <- data.crit[rownames(data), 'crit']
  }
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
  
  
  predicted.criterion <- satf.criterion(data, params$dv, params$bias, params$cnames)
  
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
  res <- optim(par=start, rcpp_compute_logLikFn, control=optim.control, method=method)  
  
  ##  # optimize a second time
  ##  res <- optim(par=start, rcpp_compute_logLikFn, control=optim.control, method=method)  
  
  # transform the parameters into their proper (constrained) form
  res$par <- rcpp_constrain_coefs( res$par )
  
  # free all memory reserved by C++
  rcpp_deinitialize_logLikFn()
  
  return(res)
}
