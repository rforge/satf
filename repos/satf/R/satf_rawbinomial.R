#library(pracma)
library(ggplot2)
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

satf.rawbinomial <- function(dv, signal, start, contrasts, constraints, data,
                              control, bias, cnames,
                             corr.mrsat=FALSE, plot=FALSE, optim.digits=1, 
                             optim.control=list(), method="Nelder-Mead", ...)
{
  # Check 'data' and 'dv' parameters.
  reportifnot(nrow(data) > 0, "Parameter 'data' needs to consist of several rows.")

  predicted.criterion <- satf.criterion(data, dv, bias, cnames)
    
  # Create design matrix.
  optim.control$fnscale <- -1

  initialize.compute.logLik(dv=dv, contrasts=contrasts, constraints=constraints, 
                            data=data, predicted.criterion=predicted.criterion,
                            cnames=cnames)
  start = untransform.coefs(start)
  res <- compute.logLik(params=start)
  if(res == -Inf || res == Inf) {
    deinitialize.compute.logLik()
    return(NULL)
  }

  # maximize the log-likelihood
  res <- optim(par=start, compute.logLik, control=optim.control, method=method)  
  
  if(corr.mrsat == TRUE)
  {
    start <- c(start, trial.corr=p2logodds(.5))
    res <- compute.logLik(start)
    if(res == -Inf || res == Inf) {
      return(NULL)
    }
    res <- optim(start, fn=compute.logLik, fixedparams=fp, control=optim.control, method=method)
  }
  
  res <- optim.to.precision(start, optim.digits=optim.digits, fn=compute.logLik, control=optim.control)
  
  # free all memory reserved by C for optimization
  deinitialize.compute.logLik()
  
  res$par <- transform.coefs(res$par)
#  res$fn.criterion <- fn.criterion
  res
}
