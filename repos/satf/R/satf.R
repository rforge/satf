.packageName <- "satf"


library(plyr)

Deviance <- function(LL) -2*LL
AIC <- function(LL, k) Deviance(LL) + 2*k

p2logodds <- function(p) log(p/(1-p))
logodds2p <- function(lodds) exp(lodds)/(1+exp(lodds))

trial.noise.range <- c(.1, .85)
decode.trial.noise <- function(x) (exp(x)/(1+exp(x)) )*diff(trial.noise.range)+trial.noise.range[1]
encode.trial.noise <- function(p) log(( (p-trial.noise.range[1])/diff(trial.noise.range) )/(1- (p-trial.noise.range[1])/diff(trial.noise.range)) ) 

NumSmallest <- -.Machine$double.xmax
NumLargest <- .Machine$double.xmax

init.satf.params <- function(start, lambda, beta, delta,
                       pslope.lambda=FALSE, pslope.beta=FALSE, pslope.delta=FALSE, 
                       nslope.lambda=FALSE, nslope.beta=FALSE, nslope.delta=FALSE,
                       min.slope.lambda=NumSmallest, min.slope.beta=NumSmallest, min.slope.delta=NumSmallest,
                       max.slope.lambda=NumLargest, max.slope.beta=NumLargest, max.slope.delta=NumLargest,
                       ...)
{
  stopifnot(!is.null(start))
  stopifnot(!(pslope.lambda==T && nslope.lambda==T))
  stopifnot(!(pslope.beta  ==T && nslope.beta  ==T))
  stopifnot(!(pslope.delta ==T && nslope.delta ==T))
  
  cslope.lambda <- pslope.lambda || nslope.lambda
  cslope.beta <- pslope.beta  || nslope.beta
  cslope.delta <- pslope.delta || nslope.delta
  
  # convert into vector format if necessary
  if(is.list(start)) {
    start <- c(start$satf.params, start$crit.params, start$repeat.params)
  }

  start <- start[!is.na(start)]
  # init all necessary satf parameters
  newstart <- c(start[c('lambda','beta','delta')],
             add.params(lambda,'lambda', start, 0),
             add.params(beta,'beta', start, 0),
             add.params(delta,'delta', start, 0))
  remainingstart <- start[ -which(names(newstart) %in% names(start)) ]
  start <- newstart

  # set transformation flags
  constrained <- c(T, T, T, rep(cslope.lambda, n.terms(lambda)), rep(cslope.beta, n.terms(beta)), rep(cslope.delta, n.terms(delta)))
  sgn <- c(F, F, F, rep(nslope.lambda, n.terms(lambda)), rep(nslope.beta, n.terms(beta)), rep(nslope.delta, n.terms(delta)))
  sgn <- -(as.integer(sgn)-.5)*2

  minima <- c(rep(0, 3), rep(min.slope.lambda, n.terms(lambda)),
               rep(min.slope.beta, n.terms(beta)),
               rep(min.slope.delta, n.terms(delta)))
  maxima <- c(rep(NumLargest, 3), rep(max.slope.lambda, n.terms(lambda)),
               rep(max.slope.beta, n.terms(beta)),
               rep(max.slope.delta, n.terms(delta)))

  start.constrained.sign <- start!=0 & constrained
  if(!(all(sign(start[start.constrained.sign]) == sgn[start.constrained.sign]))) {
    print("One of the starting values contradicts the constraints.")
    stop()
    warning("One of the starting values contradicts the constraints.")
    return(NULL)
  }

  # transform the appropriate starting values
  start[constrained] <- start[constrained]*sgn[constrained]
  start[constrained] <- sqrt(start[constrained])
  
  # determine transformation
  transform.fn <- function(params) {
    n.extend <- length(params)-length(constrained)
    constrained <- c(constrained, rep(F, n.extend))
    sgn <- c(sgn, rep(1, n.extend))
    params[constrained] <- (params[constrained])^2
    params*sgn
  }
  
  # add the non-satf starting values
  if(length(remainingstart)) {
    overlap <- names(remainingstart) %in% names(start)
    if(any(overlap)) {
      overlap.txt <- paste("'", names(remainingstart)[overlap], "'", collapse=",")
      warning.msg <- paste("Ignoring start values for:", overlap.txt)
      warning(warning.msg)
      remainingstart <- remainingstart[!overlap]
    }
    start <- c(start, remainingstart)
  }

  list(start=start, transform.fn=transform.fn, transform.vec=constrained*sgn,
       params.minimum=minima, params.maximum=maxima)
}


satf  <- function(dv, signal, lambda, beta, delta, start, data, metric,
                  control=list(), plot=FALSE, ...)
{
  stopifnot(metric %in% c('RMSD','R2','adjR2','logLik','logLikRaw'))
  init <- init.satf.params(start=start, lambda, beta, delta, ...)
  if(is.null(init))
    return(NULL)
  
  # init default control parameters
  control <- default(control, 'correction.fraction', .01)
  control <- default(control, 'plot.dprime', plot)
  control <- default(control, 'plot.criterion', FALSE)
  control <- default(control, 'plot.action', 'plot')
  #control <- default(control, 'crit.spline.step', 0.1)
  control <- default(control, 'condition', ~condition)
  control <- default(control, 'interval', ~interval)
  control <- default(control, 'time', ~time)

  # init default control parameters
  control$condition <- formula.terms(control$condition)
  control$time <- formula.terms(control$time)
  
  if(metric=='logLikRaw') {
    res <- Fit.SATF.RawBinomial( dv=dv, signal=signal, start=init$start,
                                lambda=lambda, beta=beta, delta=delta,
                                transform.fn=init$transform.fn,
                                transform.vec=init$transform.vec,
                                params.minimum=init$params.minimum,
                                params.maximum=init$params.maximum,
                                data=data, control=control, ...)
  } else {
    res <- Fit.SATF.Aggregate(start=init$start, transform.fn=init$transform.fn,
                              transform.vec=init$transform.vec,
                              params.minimum=init$params.minimum,
                              params.maximum=init$params.maximum,
                              lambda=lambda, beta=beta, delta=delta,
                              data=data, metric=metric, control=control, ...)
  }
  return(res)
}

satf.group <- function(start, data, metric="logLik", multilevel=T, ...)
{
  stopifnot(metric %in% c("logLik", "logLikRaw"))
  library(plyr)
  init.params <- ddply(data, .(subject), function(d.cur) {
      print(paste("subject", d.cur$subject[1]))
      res <- satf(start=start, data=d.cur, metric=metric, lambda=~1, beta=~1, delta=~1)
      start <- res$par$satf.params
      res <- satf(start=start, data=d.cur, metric=metric, ...)
      res <- c(LL=res$value, unlist(res$par))
      res
  })
  if(!multilevel) {
    return(list(params=init.params, hyperparams=NULL))
  }
  init.params$LL.params <- NULL
  init <- init.params
  init$subject <- NULL
  init$LL <- NULL
  hyperparams <- data.frame()
  init.names <- colnames(init)
print(init)
  
  for(i in 1:ncol(init)) {
    name <- colnames(init)[i] 
    short.name <- strsplit(name, split=".params.")[[1]][2]
    # exclude outlying values
    x <- remove.outliers(init[,name], factor=3)
    m <- mean(x); v <- var(x)/2
    sgn <- 0
    if(all(init[,name] <= 0)) sgn = -1
    else if(all(init[,name] >= 0)) sgn = 1
    if(m == 0) m = 0.1
    scale <- v/m; shape <- abs(m/scale)
    hyperparams <- rbind(hyperparams, data.frame(param=name, param.short=short.name,
                                      hyperparam="scale", val=scale, constrained=sgn,
                                           offset=0.001))
    hyperparams <- rbind(hyperparams, data.frame(param=name, param.short=short.name,
                                      hyperparam="shape", val=shape+1, constrained=1,
                                           offset=1))
    #hyperparams <- rbind(hyperparams, data.frame(param=name, param.short=short.name,
    #                                  hyperparam="m", val=median(x), constrained=sgn,
    #                                  offset=0))
    #hyperparams <- rbind(hyperparams, data.frame(param=name, param.short=short.name,
    #                                  hyperparam="sd", val=sd(x)/4, constrained=1,
    #                                  offset=0))
  }
print(hyperparams)

  colnames(init.params) <- sapply(strsplit(colnames(init.params), split=".params."), function(x) x[length(x)])
  start.names <- (colnames(init.params))[-c(1:2)]
  start.params.bySubj <- dlply(init.params, .(subject),  function(d) unlist(d[1,-c(1:2)]) ) 
  
  transform.hyperparams <- function(hyperparams) within(hyperparams, {
    val[constrained!=0] <- val[constrained!=0] - offset[constrained!=0]
    val[constrained==0] <- val[constrained==0] - sign(val[constrained==0])*offset[constrained==0]
    val[constrained!=0] <- log(val[constrained!=0]*constrained[constrained!=0])
  })
  untransform.hyperparams <- function(hyperparams) within(hyperparams, {
    val[constrained!=0] <- constrained[constrained!=0]*exp(val[constrained!=0])
    val[constrained!=0] <- val[constrained!=0] + offset[constrained!=0]
    val[constrained==0] <- val[constrained==0] + sign(val[constrained==0])*offset[constrained==0]
  })
  
  Eval.SATF.MLV <- function(hyperparams.vec, params.indices=nrow(hyperparams.vec),
                            subj.params=F, plot.params=F) 
  {
    hyperparams$val[params.indices] <- hyperparams.vec
    hyperparams <- untransform.hyperparams(hyperparams)
 ###   
 #   print.params <- hyperparams$val
 #   names(print.params) <- hyperparams$param.short
 #   print("current params")
 #   print(print.params)
 ###    
    fn.hyperparams <- dlply(hyperparams, .(param), function(d) {
      if(all(d$hyperparam == c("scale","shape"))) {
        scale <- d$val[1]; shape <- d$val[2];
        fn <- function(x) dgamma(x, scale=scale, shape=shape, log=T)
      } else if(all(d$hyperparam == c("m","sd"))) {
        m <- d$val[1]; sd <- d$val[2];
        fn <- function(x) dnorm(x, mean=m, sd=sd, log=T)
      } else {
        stop("This distribution is not implemented.")
      }
    })
    hyperparams1 <- hyperparams$val[c(T,F)]
    hyperparams2 <- hyperparams$val[c(F,T)]
    parameter.sign <- sign(hyperparams1)
    fn.params.LL <- function(params=NULL, action="compute.LL") {
      if(action=="compute.LL") {
        params <- params$satf.params
        params <- params*sign(hyperparams1)
        hyperparams1 <- abs(hyperparams1)
        LLs <- dgamma(params, scale=hyperparams1, shape=hyperparams2, log=T)
        return( sum(LLs) )
      }
      else if(action=="compute.LLs") {
        params <- params$satf.params
        params <- params*sign(hyperparams1)
        hyperparams1 <- abs(hyperparams1)
        LLs <- dgamma(params, scale=hyperparams1, shape=hyperparams2, log=T)
        return( LLs )
      } else if(action=="return.hyperparams1") {
        return(hyperparams1)
      } else if(action=="return.hyperparams2") {
        return(hyperparams2)
      }
    }

    global.start <- start.params.bySubj[[ asc(d.cur$subject[1]) ]]
    fits <- ddply(data, .(subject), function(d.cur) {
      start <- start.params.bySubj[[ asc(d.cur$subject[1]) ]]
      #print("start")
      #print(start)
      #print("/start")
      
      wrong.sign <- which(parameter.sign != sign(start))
      start[wrong.sign] <- 0.01*parameter.sign[wrong.sign]
      ###startLL <- fn.params.LL(list(satf.params=cur.start),  action="compute.LLs")
      ###start.impossible <- which(abs(startLL)==Inf)
      ###start[start.impossible] <- start.impossible+rnorm(length(start.impossible), sd=5)
      ###print(start)
      res <- satf(start=start, data=d.cur, metric=metric,
                    fn.params.LL=fn.params.LL, optim.digits=NA, ...)
      if(is.null(res)) {
        stop("Failed inital evaluation even after resetting params.")
      }
   ###   print("param LL")
   ###   print( res$par$satf.params[2] )
   ###   print( fn.params.LL(res$par,  action="compute.LLs")[2] )
      return( c(LL=res$value, unlist(res$par)) )
      #print(res$par$satf.params-start)
    }, .parallel=T)
    
    #print("fits")
    #print(fits)
    if(nrow(fits) != nrow(init.params)) {
      return(-Inf)
    }
    if(plot.params) {
      print("params")
      print(sum(fits$LL.params))
      #print(round(fits$LL.params,2))
      p1 <- ggplot(data=fits, aes(x=satf.params.lambda, y=0))+geom_point()+ 
        stat_function(fun=function(x) exp((fn.hyperparams[['satf.params.lambda']])(x)))+
          scale_x_continuous(limits=c(0,10))
      p2 <- ggplot(data=fits, aes(x=satf.params.beta, y=0))+geom_point()+ 
        stat_function(fun=function(x) exp((fn.hyperparams[['satf.params.beta']])(x)) )+
          scale_x_continuous(limits=c(0,5))
      p3 <- ggplot(data=fits, aes(x=satf.params.delta, y=0))+geom_point()+ 
        stat_function(fun=function(x) exp((fn.hyperparams[['satf.params.delta']])(x)))+
          scale_x_continuous(limits=c(0,1))
      print(multiplot(p1,p2,p3))
    }
    
    #print( fits$LL )
    print( sum(fits$LL) )
    if(subj.params)
      res <- fits
    else
      res <- sum(fits$LL)
    return(res)
  }

  do.plot <- T  
  start.params <- function(hyperparams) {
    params <- Eval.SATF.MLV(hyperparams$val, subj.params=T, plot.params=do.plot)
    params <- dlply(params, .(subject),  function(d) {
      params <- unlist(d[1,-c(1:2,ncol(params))])
      names(params) <- start.names
      params
    })
    params
  }

  # determine new start parameters for participants
  hyperparams <- transform.hyperparams(hyperparams)
  indices.all <- seq(1, nrow(hyperparams), 1)
  start.params.bySubj <- start.params(hyperparams)
  
  # Fit parameters separately, first.
  for(i in seq(1, nrow(hyperparams), 2)) {
    indices <- c(i, i+1)
    print(paste("Fitting parameter ",(i+1)/2," of ",nrow(hyperparams)/2,".", sep=""))
    res <- optim(hyperparams$val[indices], Eval.SATF.MLV, params.indices=indices,
                 plot.params=F, control=list(fnscale=-1, maxit=100, reltol=1e-4))
    hyperparams$val[indices] <- res$par
    start.params.bySubj <- start.params(hyperparams)
    print(res)
  }
  
  print("Fitting all parameters, reltol=10e-4.")
  res <- optim(hyperparams$val[indices.all], Eval.SATF.MLV,
                params.indices=indices.all, 
		control=list(fnscale=-1, maxit=500, reltol=1e-4))
  print(res)
  hyperparams$val <- res$par
  start.params.bySubj <- start.params(hyperparams)
  convergence4 <- res$convergence 

  print("Fitting all parameters, reltol=10e-5.")
  res <- optim(hyperparams$val[indices.all], Eval.SATF.MLV,
                params.indices=indices.all, 
		control=list(fnscale=-1, maxit=500, reltol=1e-5))
  print(res)
  hyperparams$val <- res$par
  start.params.bySubj <- start.params(hyperparams)
  convergence5 <- res$convergence 

  print("Fitting all parameters, reltol=10e-6.")
  res <- optim(hyperparams$val[indices.all], Eval.SATF.MLV,
                params.indices=indices.all, 
		control=list(fnscale=-1, maxit=500, reltol=1e-6))
  print(res)
  hyperparams$val <- res$par
  start.params.bySubj <- start.params(hyperparams)
  convergence6 <- res$convergence 

  print("Fitting all parameters, reltol=10e-7.")
  res <- optim(hyperparams$val[indices.all], Eval.SATF.MLV,
                params.indices=indices.all, 
		control=list(fnscale=-1, maxit=500, reltol=1e-7))
  print(res)
  hyperparams$val <- res$par
  start.params.bySubj <- start.params(hyperparams)
  convergence7 <- res$convergence 

  print("Fitting all parameters, reltol=1e-8.")
  res <- optim(hyperparams$val[indices.all], Eval.SATF.MLV,
                params.indices=indices.all, 
		control=list(fnscale=-1, maxit=10^6, reltol=1e-8))
  print(res)
  hyperparams$val <- res$par
  convergence8 <- res$convergence 

  subj.params <- Eval.SATF.MLV(hyperparams$val, subj.params=T, plot.params=do.plot)
  hyperparams <- untransform.hyperparams(hyperparams)
  hyperparams$LL <- res$value  
  hyperparams$'convergence.1e-4' <- convergence4
  hyperparams$'convergence.1e-5' <- convergence5
  hyperparams$'convergence.1e-6' <- convergence6
  hyperparams$'convergence.1e-7' <- convergence7
  hyperparams$'convergence.1e-8' <- convergence8
  list(params=subj.params, hyperparams=hyperparams)
}
