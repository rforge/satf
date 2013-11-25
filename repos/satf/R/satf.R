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




init.satf.params <- function(start, contrasts, constraints, trial.id) {

  reportifnot(is.vector(start) && !is.list(start), "Parameter 'start' needs to be a vector.")
  reportifnot(!is.null(start) || !any(is.na(start)), "Parameter 'start' was not provided or containts NAs.")

  start.original <- start
  start.coreparams <- start[names(contrasts)]
  reportifnot(!any(is.na(start.coreparams)), sprintf("Parameter 'start' needs to contain start values for: %s", 
                                             paste(names(contrasts), collapse=", ") ))

  # init all necessary satf parameters
  newstart <- start.coreparams
  for(name in names(contrasts))
    newstart <- c(newstart, add.params(contrasts[[name]], name, start, NA))

  # add a corr.mrsat parameter if necessary
  if(!is.null(trial.id))
    newstart <- c(newstart, corr.mrsat=NA)
  
  start <- newstart

  # set default constraints for intercepts
  for(name in names(contrasts))
    constraints <- default(constraints, name, c(0,Inf))

  # initialize constraints matrix
  constraint.matrix <- matrix(c(-Inf, Inf), nrow=length(start), ncol=2, byrow=TRUE, 
                              dimnames=list(names(start), c('lower', 'upper')))

  # define functions for finding coefficient constraints in the constraint list
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
  warn.about.start.params <- FALSE
  for(name in names(start)) {
    val <- find.match(name)
    if( is.null(val) )
      next

    # set value
    constraint.matrix[name,] <- val
    
    # if the constaint specifies a fixed value, change corresponding start parameter    
    if( length(val) != 1 )
      next
    if( is.na(start[[name]]) ) {
      start[[name]] <- val
      
    } else if(start[[name]] != val) 
    {
      start[[name]] <- val
      warn.about.start.params <- TRUE
    }
  }

  # make sure the bounds on 'corr.mrsat' are reasonable, i.e. within [0; 1]
  if('corr.mrsat' %in% rownames(constraint.matrix)) {
    if( constraint.matrix['corr.mrsat', 'lower'] < 0 )
      constraint.matrix['corr.mrsat', 'lower'] <- 0
    if( constraint.matrix['corr.mrsat', 'upper'] > 1 )
      constraint.matrix['corr.mrsat', 'upper'] <- 1
    
    lower <- constraint.matrix['corr.mrsat', 'lower']
    upper <- constraint.matrix['corr.mrsat', 'upper']
    if('corr.mrsat' %in% names(start.original))
      start[['corr.mrsat']] <- start.original[['corr.mrsat']]
    else
      start[['corr.mrsat']] <- (lower+(upper-lower)/2)
  }
  
  start[is.na(start)] <- 0
  if(warn.about.start.params) {
    start.str <- paste0(names(start), '=', start,  collapse=", ")
    warning(sprintf("Start values changed to <%s>.", start.str))
  }
  
  # check whether start parameters are within the limits specified by the constraints
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

  satf.params <- init.satf.params(start=start, contrasts=contrasts, constraints=constraints, trial.id=trial.id)
  
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
