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




init.satf.params <- function(start, contrasts, constraints) {

  reportifnot(is.vector(start) && !is.list(start), "Parameter 'start' needs to be a vector.")
  reportifnot(!is.null(start) || !any(is.na(start)), "Parameter 'start' was not provided or containts NAs.")

  start.coreparams <- start[names(contrasts)]
  reportifnot(!any(is.na(start.coreparams)), sprintf("Parameter 'start' needs to contain start values for: %s", 
                                        paste(names(contrasts), collapse=", ") ))

  # init all necessary satf parameters
  newstart <- start.coreparams
  for(name in names(contrasts))
    newstart <- c(newstart, add.params(contrasts[[name]], name, start, 0))

  # not sure what the following line is supposed to do
  # remainingstart <- start[ -which(names(newstart) %in% names(start)) ]
  start <- newstart

  # set default constraints for intercepts
  for(name in names(contrasts))
    constraints <- default(constraints, name, c(0,Inf))

  # initialize constraint matrix
  constraint.matrix <- matrix(c(-Inf, Inf), nrow=length(start), ncol=2, byrow=TRUE, 
                              dimnames=list(names(start), c('lower', 'upper')))

  is.match <- function(expr, name) regexpr(paste0('^',expr,'$'), text=name) > 0
  find.match <- function(name) {
    for(expr in names(constraints)) {
      if(is.match(expr, name)) {
        return(constraints[[expr]])
      }
    }
    return(NULL)
  }

  # fill the constraints matrix
  for(name in names(start)) {
    val <- find.match(name)
    if(!is.null(val)) {
      constraint.matrix[name,] <- val
    }
  }

  tmp <- cbind(start, constraint.matrix)
  allwithinlimits <- all( tmp[,'start'] >= tmp[,'lower'] & tmp[,'start'] <= tmp[,'upper'] )
  if( !allwithinlimits ) {
    print(tmp)
    stop("Start parameters are not within limits.")
  }
  return(list(constraints=constraint.matrix, start=start));
}


# translate parameters from formula notation to string notation,
# and check that all columns actually exist
# TODO: This function must be rather slow, optimize.

translate.parameters  <- function(data, dv, start, contrasts, constraints, 
                                  bias, signal, time, trial.id, summarize) {
  
  cnames <- list()
  satf.params <- init.satf.params(start=start, contrasts=contrasts, constraints=constraints)
  
  # TRANSLATE PARAMETER: dv
  if('response' %in% names(dv)) {
    dv[['response']] <- check.formula.for.colname(data, 'dv$response', dv$response)
    stopifnot.binary( data[, dv[['response']] ] )
    
  } else if( all(c('n.responses.yes','n.responses') %in% names(dv)) ) {
    for(name in names(dv))
      dv[[name]] = check.formula.for.colname(data, paste('dv', name, sep='$'), dv[[name]])
    
  } else {
    stop("No such dv implemented.")
    
  }
  dv <- unlist(dv)
  
  # TRANSLATE PARAMETER: contrasts
  for(name in names(contrasts))
    contrasts[[name]] = check.formula.for.colname(data, paste('contrasts', name, sep='$'), contrasts[[name]], n.cols=NA)

  # TRANSLATE PARAMETER: bias
  reportifnot(all(c('bias', 'poly.degree') %in% names(bias)), "Parameter 'bias' needs to containt 'bias' and 'poly.degree'.")
  bias[['bias']] = check.formula.for.colname(data, 'bias$bias', bias[['bias']], n.cols=NA)
  
  # CHECK COLUMN NAMES: signal, time, trial.id
  cnames$signal <- check.formula.for.colname(data, 'signal', signal)
  stopifnot.binary( data[,cnames$signal] )
    
  cnames$time <- check.formula.for.colname(data, 'time', time)

  if(!is.null(trial.id))
    cnames$trial.id <- check.formula.for.colname(data, 'trial.id', trial.id)

  # RETURN
  list(dv=dv, start=satf.params$start, contrasts=contrasts, bias=bias,
       constraints=satf.params$constraints, cnames=unlist(cnames) )
}
                  
satf  <- function(dv, signal, start, contrasts, data, time, metric, bias, 
                  trial.id=NULL, summarize=list(),
                  constraints=list(), control=list(), plot=FALSE, optim.control=list(), ...)
{
  metric.permissible <- c('RMSD','R2','adjR2','logLik','logLikRaw')
  reportifnot(metric %in% metric.permissible, sprintf("'metric' has to be one of: %s", paste(metric.permissible, collapse=", ")))

  # init default parameters
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
  
  # check parameters
  params <- translate.parameters(data=data, dv=dv, start=start, contrasts=contrasts, 
                                 constraints=constraints, bias=bias, signal=signal, 
                                 time=time, trial.id=trial.id, summarize=summarize)
  
  
  res <- satf.rawbinomial( dv=params$dv, start=params$start, contrasts=params$contrasts, 
                           constraints=params$constraints, bias=params$bias, data=data, control=control, 
                           optim.control=optim.control, cnames=params$cnames, ...)
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
