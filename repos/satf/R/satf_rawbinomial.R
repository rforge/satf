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

satf.rawbinomial <- function(dv, signal, start, contrasts, constraints, data, control, bias, cnames,
                             plot=FALSE, optim.digits=1, optim.control=list(), method="Nelder-Mead", ...)
{
  # Check 'data' and 'dv' parameters.
  reportifnot(nrow(data) > 0, "Parameter 'data' needs to consist of several rows.")
  optim.control$fnscale <- -1

  predicted.criterion <- satf.criterion(data, dv, bias, cnames)
  
  # initialize the C++ optimization routine
  rcpp_initialize_logLikFn(dv, contrasts, constraints, 
                            data, predicted.criterion,
                            cnames)
  
  start <- rcpp_unconstrain_coefs( start )

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
  
  res
}
