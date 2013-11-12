#library(pracma)
library(ggplot2)
library(reshape)

# TODO: Fix all hardcoded references to column names. They are easy to find, because
# the contain '$'.


# NOTE: This function is not used currently as it has been replaced by a C implementation

# What I call criterion c here, is what Wickens calls lambda in chapter 2
# (http://cognitrn.psych.indiana.edu/busey/q270/pdfs/wickensSDT.pdf)
Eval.SATF.RawBinomial <- function(params, fixedparams)
{
  response.yes <- fixedparams$response.yes
  last.response.yes <- fixedparams$last.response.yes
  last.response.na <- fixedparams$last.response.na
  data <- fixedparams$data
# point 0

  # process parameters
  params <- process.parameters(params=params, fixed=fixedparams$fixed, 
                               params.idx=fixedparams$params.idx, 
                               transform.fn=fixedparams$transform.fn,
                               fn.params.LL=fixedparams$fn.params.LL)

  # compute criterion predictions
  if(is.null(fixedparams$predicted.criterion)) {
     if(!is.null(fixedparams$fn.criterion)) {
       fixedparams$predicted.criterion <- fixedparams$fn.criterion(data$time, params$crit.params)
     } else {
       print( names(fixedparams) )
       stop("No estimate of the criterion provided.")
     }
  }

  # compute d' predictions from the SAT function 
  predicted.dprime <- predict.dprime(params$satf.params, fixedparams$design.matrix, data)
  
  if('trial.noise.s' %in% names(params$noise.params)) {
    
    # determine differences in d' and criterion
    diff.crit.psi <- predicted.dprime-predicted.criterion
    last.diff.crit.psi <- c(0, diff.crit.psi[-nrow(data)])
    last.diff.crit.psi[last.response.na==1] <- 0

    trial.noise.s <- params$noise.params[['trial.noise.s']]
    LLs <- sapply(fixedparams$prediction.idx, function(idx) {
      stopifnot(length(unique(predicted.dprime[idx]))==1)
      stopifnot(length(unique(predicted.criterion[idx]))==1)
      stopifnot(length(unique(last.diff.crit.psi[idx]))==1)
      stopifnot(length(unique(last.response.yes[idx]))==1)
      stopifnot(length(unique(last.response.na[idx]))==1)
      i <- idx[1]
      if( last.response.na[i] ) {
        probYes <- 1-pnorm(predicted.criterion[i], mean=predicted.dprime[i])
      } else {
        probYes <- prob.dependent.yes(diff.crit.psi[i], trial.noise.s, 
                                      last.diff.crit.psi[i], last.response.yes[i])
        stopifnot(probYes > 0 && probYes < 1)
      }
      probResponse <- response.yes[idx]*probYes + (1-response.yes[idx])*(1-probYes)
      #stopifnot(length(probResponse) == length(idx))
      c(LL=sum(log(probResponse)))
    })
    logLik <- sum(LLs)
    #print(logLik)
  } 
  else 
  {
    # # compute probabilities of a 'yes' response
    # probYes <- 1-pnorm(predicted.criterion, mean=predicted.dprime)
    # # compute probabilities of the given response
    # probResponse <- response.yes*probYes + (1-response.yes)*(1-probYes)
    probResponse <- response.yes -(2*response.yes-1)*pnorm(fixedparams$predicted.criterion, mean=predicted.dprime)
    logLik <- sum(log(probResponse))    
  }
  logLik <- logLik + params$LL.params
  logLik
}


Fit.SATF.RawBinomial.FixedCriterion <- function(dv.column, signal, 
        lambda, beta, delta, start, fixed, criterion, fn.criterion, data,
        transform.fn, transform.vec, params.minimum, params.maximum, corr.mrsat, control,
        optim.digits=1, fn.params.LL=NULL, optim.control=list(), method="Nelder-Mead", ...)
{
  stopifnot(criterion %in% c('fit','stepfun'))

  optim.control <- default(optim.control, 'maxit', 10^6)
  
  # make sure that all response-related parameters are as expected
  check.columns(dv.column, data)
  response.yes <- as.integer(data[,dv.column])
  stopifnot.binary( response.yes )
  
  # index different parameters
  params.idx <- index.params( c(start, fixed) )
  
  optim.control$fnscale <- -1
  dm <- design.matrix(lambda=lambda, beta=beta, delta=delta, data=data,
                      stimulus.signal=signal)

  # set all the parameters to be passed to the evaluating function
  fixedparams <- list(design.matrix=dm, fixed=fixed, data=data, 
                 fn.criterion=fn.criterion, response.yes=response.yes, 
                 transform.fn=transform.fn, params.idx=params.idx,
                 fn.params.LL=fn.params.LL)
  fixedparams$predicted.criterion <- fn.criterion(data$time)

  # set the relevant parameters necessary for optimization
  fp <- initialize.compute.logLik(design.matrix=dm, transform.vec=transform.vec,
                    params.minimum=params.minimum, params.maximum=params.maximum, 
                    time=data$time, predicted.criterion=fn.criterion(data$time),
                    trial.id=data$trial.id, response.yes=response.yes)
  res <- compute.logLik(start, fixedparams=fp)
  if(res == -Inf || res == Inf) {
    deinitialize.compute.logLik(fp)
    return(NULL)
  }
  # maximize the log-likelihood
  res <- optim(par=start, compute.logLik, fixedparams=fp, coefs.len=length(start),
               control=optim.control, method=method)  
 
  if(criterion == 'fit' && !any(params.idx$idx.crit.params)) {
    fn.params.start <- c(c1=1, c2=0, c3=0, c4=0)
    CriterionFnExp <- function(x, p) with(data.frame(t(p)), (c3*exp(c1*x+c2)/(1+exp(c1*x+c2)))+c4 )
    fit.tmp.ts <- seq(min(data$time), max(data$time), length.out=30)
    fit.tmp.fn <- function(p) {
      sum((CriterionFnExp(fit.tmp.ts, p)-fn.criterion(fit.tmp.ts))^2)
    }
    res.crit <- optim(fn.params.start, fit.tmp.fn, control=list(maxit=10^6))
    fn.criterion.new <- function(t) CriterionFnExp(t, res.crit$par)
   
    start <- c(res$par, res.crit$par)
    fixedparams$params.idx <- index.params( c(start, fixedparams$fixed) )
    fixedparams$fn.criterion <- CriterionFnExp
    fixedparams$predicted.criterion <- NULL
    # TODO: Fix the reference to Eval.SATF.RawBinomial
    res <- optim(par=start, Eval.SATF.RawBinomial, fixedparams=fixedparams,
                 coefs.len=length(start), control=optim.control, method=method)
    fn.criterion <- function(t) CriterionFnExp(t, res$par)
  } 
  
  # TODO: See if it's possible to find a closed-form solution for 'trial.noise.s'
  # given d' and criterion functions.
  if(corr.mrsat == TRUE)
  {
    start <- c(start, trial.corr=p2logodds(.5))
    res <- compute.logLik(start, fixedparams=fp)
    if(res == -Inf || res == Inf) {
      return(NULL)
    }
    
    res <- optim(start, fn=compute.logLik, fixedparams=fp,
                 coefs.len=length(start), control=optim.control,
                 method=method)
  }

  res <- optim.to.precision(start, optim.digits=optim.digits,
               fn=compute.logLik, fixedparams=fp, coefs.len=length(start),
               control=optim.control)
  
  # free all memory reserved by C for optimization
  deinitialize.compute.logLik(fp)
  
  res$par <- process.parameters(res$par, fixed, fixedparams$params.idx,
                                transform.fn, params.minimum, params.maximum)
  res$fn.criterion <- fn.criterion
  
  if(control$plot.dprime)
  { # TODO: move plotting to the main function in order to compare fits with other metrics
    min.time <- min(data$time)
    interval.size <- diff(range(data$time))/plot.fit$intervals
    time <- floor((data$time-min.time)/interval.size)*interval.size+min.time
    corr <- plot.fit$correction.fraction
      
    # compute d' values from the data for each time interval
    cond.column <- formula.terms(control$condition)
    psi <- function(x) qnorm(1-(sum(x)+corr)/(length(x)+2*corr))
    tmp <- data.frame(response=response.yes, signal=signal, time=time,
                      condition=data[,cond.column])
    tmp <- ddply(tmp, .(time), function(d) { if(!mean(d$signal)%in%c(0,1)) return(d)  })
    tmp$condition <- asc(tmp$condition)
    tmp.noise <- subset(tmp, !as.logical(signal))
    tmp.signal <- subset(tmp, as.logical(signal))
    (psi.noise <- with(tmp.noise, tapply(response, list(condition, time), psi)))
    (psi.signal <- with(tmp.signal, tapply(response, list(condition, time), psi)))
    emprical.dprime <-  as.vector(t(psi.noise)) - t(psi.signal)
    emprical.dprime <- melt(emprical.dprime)
    colnames(emprical.dprime) <- c('time','condition','dprime')

    # obtain predicted dprime values
    predicted.dprime <- predict.dprime(res$par$satf.params, dm, data)
    predicted.dprime <- data.frame(time=data$time, condition=data[,terms], dprime=predicted.dprime)
    predicted.dprime <- predicted.dprime[as.logical(signal),]

    dprime <- rbind(cbind(emprical.dprime, type="empirical"),
                    cbind(predicted.dprime, type="predicted"))
    p <- ggplot(predicted.dprime, aes(x=time, y=dprime, linetype=condition))+
         geom_line()+ geom_point(data=emprical.dprime, aes(x=time, y=dprime, linetype=condition))
    plot(p)
  }
  res  
}


# TODO: Determine the stepfun for every ungrammatical condition

# TODO: Merge this function with the function it calls.
Fit.SATF.RawBinomial <- function(dv, signal, lambda, beta, delta, start, 
                                 fixed, data, corr.mrsat=FALSE, criterion='stepfun',
                                 plot=FALSE, control, ...)
{
  # Check that all required columns and parameters are specified.
  stopifnot( is.formula(dv) )
  dv.column = formula.terms(dv); 
  stopifnot( length(dv.column) == 1 )

  stopifnot(nrow(data) > 0)
  check.columns(formula.terms(signal), data)
  signal <- as.integer(data[,formula.terms(signal)])
  stopifnot.binary( signal )
  fixed <- c()

  if(plot) plot.fit <- list(intervals=13, condition=~condition,
                            correction.fraction=control$correction.fraction)
  else plot.fit <- list()
  
  fn.crit <- ComputeEmpiricalBiasFn(data, signal, formula.terms(dv), control)
  est <- Fit.SATF.RawBinomial.FixedCriterion( dv.column=dv.column, signal=signal,
                              lambda=lambda, beta=beta, delta=delta, 
                              start=start, fixed=fixed, control=control,
                              criterion=criterion, fn.criterion=fn.crit, 
                              corr.mrsat=corr.mrsat, plot.fit=plot.fit, 
                              data=data, ...)
  return(est)
}
