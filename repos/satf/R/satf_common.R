
# TODO: Compute criterion for each ungrammatical condition separately
ComputeEmpiricalBiasFn <- function(data, signal, response, control) 
{
  check.columns('time', data)
  corr <- control$correction.fraction
  if(all(c('fas','crs') %in% colnames(data))) {
    crs <- (data[,'crs']+corr)/(data[,'fas']+data[,'crs']+2*corr)
    time <- data$time[!is.na(crs)]
    crs <- crs[!is.na(crs)]
    criterion <- qnorm(crs)
    fn <- function(t) (splinefun(time, criterion))(t)
  } else {
    stopifnot.binary( signal )
    check.columns('interval', data)
    check.columns('condition', data)
    check.columns(response, data)
    data <- subset(data, !as.logical(signal))

    Poly15 <- function(p, t) p[1]+p[2]*t+p[3]*t^2+p[4]*t^3+p[5]*t^4+p[6]*t^5+p[7]*t^6+p[8]*t^7+p[9]*t^8+p[10]*t^9+p[11]*t^10+p[12]*t^11+p[13]*t^12+p[14]*t^13+p[15]*t^14
    pad.poly <- function(p, n=15) c(p, rep(0, n-length(p)))
    Poly <- function(p, t) Poly15(pad.poly(p), t)

    EvalPoly <- function(p, time, response.yes) {
      predicted.p <- logodds2p( Poly(p, time) )
      LL <- log( predicted.p*(response.yes==1) + (1-predicted.p)*(response.yes==0) )
      BIC <- function(LL, k, n) Deviance(LL) + k*log(n)
      BIC(sum(LL), length(p), length(response.yes))
    }

    last.par <- c(0); last.BIC <- Inf
    for(i in 1:14) {
      res <- optim(c(last.par,0), fn=EvalPoly, time=data$time, response.yes=data[,response], control=list(maxit=10^6))
      if(res$value > last.BIC) {
        break;
      }
      last.par <- res$par
      last.BIC <- res$value
    }
    bias.pars <- pad.poly(last.par)
    fn <- function(t) logodds2p( Poly15(bias.pars, t) )
  }
  return(fn)
}

design.matrix <- function(lambda, beta, delta, stimulus.signal, data)
{
  library(plyr)
  stopifnot(nrow(data) == length(stimulus.signal))
  design.matrix <- matrix(c(as.double(stimulus.signal), rep(1, 2*length(stimulus.signal))),
                          ncol=3, nrow=length(stimulus.signal))
  lambda.names <- formula.terms(lambda)
  beta.names <- formula.terms(beta)
  delta.names <- formula.terms(delta)
  check.columns(c(lambda.names,beta.names,delta.names), data)
  
  lambda <- as.matrix(data[,lambda.names], colnames=NULL)
  beta <- as.matrix(data[,beta.names], dimnames=NULL)
  delta <- as.matrix(data[,delta.names], dimnames=NULL)
  
  lambda.idx <- seq(1, len=length(lambda.names), by=1)+ncol(design.matrix)
  beta.idx <- seq(1, len=length(beta.names), by=1)+max(lambda.idx,ncol(design.matrix))
  delta.idx <- seq(1, len=length(delta.names), by=1)+max(lambda.idx,beta.idx,ncol(design.matrix))
  
  lambda.idx <- c(1, lambda.idx)
  beta.idx   <- c(2, beta.idx)
  delta.idx  <- c(3, delta.idx)

  design.matrix <- cbind(design.matrix, lambda, beta, delta)
  contrast.ids <- id(data.frame(design.matrix))
  contrast.ids <- as.integer(as.factor(contrast.ids))

  design.matrix <- t(design.matrix)

  res <- list(design.matrix=design.matrix, ncol.dm=ncol(design.matrix),
          nrow.dm=nrow(design.matrix),
          contrast.ids=contrast.ids, n.contrast.ids=max(contrast.ids),
          lambda=as.integer(lambda.idx), beta=as.integer(beta.idx),
          delta=as.integer(delta.idx),
          n.lambda=as.integer(length(lambda.idx)),
          n.beta=as.integer(length(beta.idx)),
          n.delta=as.integer(length(delta.idx))
         )
  res
}

predict.dprime <- function(satf.params, design.matrix, data) {
  coef <- satf.params*design.matrix$design.matrix
  lambda <- col.sums(coef[design.matrix$lambda,], length(design.matrix$lambda))
  beta <- col.sums(coef[design.matrix$beta,], length(design.matrix$beta))
  delta <- col.sums(coef[design.matrix$delta,], length(design.matrix$delta))
  SATF(data$time, lambda, beta, delta)
}




initialize.compute.logLik <- function(design.matrix, transform.vec,
                                   params.minimum, params.maximum, time, predicted.criterion,
                                   trial.id=NULL, response.yes=NULL,
                                   response.hits=NULL, n.grammatical=NULL,
                                   response.fas=NULL, n.ungrammatical=NULL,
                                   fn.params.LL=NULL)
{
  # set the core parameters
  params <- list(fixedparams = NULL, dm = design.matrix, 
            transform.vec = transform.vec, params.min = params.minimum, params.max = params.maximum,
            predicted.criterion = predicted.criterion, time = time )
  
  # set hyperparams
  hyperparams1 = NULL; hyperparams2 = NULL
  if(!is.null(fn.params.LL)) {
    hyperparams1 <- fn.params.LL(action="return.hyperparams1")
    hyperparams2 <- fn.params.LL(action="return.hyperparams2")
    params$hyperparams <- cbind(hyperparams1, hyperparams2)
  }
    
  # set response variable
  if(all(!is.null(response.yes)))
  {
    stopifnot( length(response.yes) > 1  )
    stopifnot( length(response.yes) != length(trial.id)  )
    params$rv <- list(rv.type=1, response.yes=response.yes, )
    
  } else if( all(!is.null(response.hits)) && all(!is.null(n.grammatical)) &&
             all(!is.null(response.fas)) && all(!is.null(n.ungrammatical))  ) 
  {
    params$rv <- list(rv.type=2, 
                      n.response.hits=response.hits, n.grammatical=n.grammatical, 
                      n.response.fas=response.fas, n.ungrammatical=n.ungrammatical)
  } 
  else {
    stop("No response variable provided.")
  }

  rcpp_initialize_logLikFn( params )
}

deinitialize.compute.logLik <- function() {
  rcpp_deinitialize_logLikFn( params )
}

compute.logLik <- function(satf.params, fixed.params=NULL)
{
  rcpp_initialize_logLikFn( satf.params, fixed.params=NULL)
}



process.parameters <- function(params, fixed, params.idx,
                               transform.fn, params.minimum, 
                               params.maximum, fn.params.LL=NULL) {
  params <- c(params, fixed)
  LL <- 0
  
  # extract response repetition probability parameters
  noise.params <- params[params.idx$noise.params]
  if('trial.noise.s' %in% names(noise.params))
    noise.params[['trial.noise.s']] <- decode.trial.noise(noise.params[['trial.noise.s']])
  
  # extract response repetition probability parameters
  crit.params <- params[params.idx$crit.params]
  # extract satf parameters
  satf.params <- params[params.idx$satf.params]
  satf.params <- transform.fn(satf.params)
  satf.params <- pmax(pmin(satf.params, params.maximum), params.minimum)
  
  if(!is.null(fn.params.LL)) {
    LL <- fn.params.LL(list(satf.params=satf.params, 
                            noise.params=noise.params,
                            crit.params=crit.params))
  }
  list(satf.params=satf.params, 
       noise.params=noise.params,
       crit.params=crit.params, 
       LL.params=LL)
}

index.params <- function(params) {
  idx.noise.params  <- names(params) %in% c('trial.noise.s')
  idx.crit.params <- names(params) %in% paste('c', 1:20, sep='')
  idx.satf.params <- !idx.noise.params & !idx.crit.params
  data.frame(noise.params=idx.noise.params,
             crit.params=idx.crit.params,
             satf.params=idx.satf.params)
}
