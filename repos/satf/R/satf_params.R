
add.params <- function(formula, letter, startpars, start=0) {
  n <- n.terms(formula)
  if(n > 0) {
    names <- paste(letter, formula.terms(formula), sep='.')
    vec <- ifelse(!is.na(startpars[names]), startpars[names], start)
    names(vec) <- names
    return(vec)
  }
  return(c())
}

init.satf.params <- function(start, contrasts, constraints, trial.id) {

  reportifnot(is.vector(start) && !is.list(start), "Parameter 'start' needs to be a vector.")
  reportifnot(!is.null(start) || !any(is.na(start)), "Parameter 'start' was not provided or containts NAs.")

  # make sure start values were provided for the intercepts
  start.intercepts <- start[names(contrasts)]
  missing.intercepts <- paste(names(start.intercepts[is.na(start.intercepts)]), collapse=", ")
  reportifnot( missing.intercepts == "", sprintf("Parameter 'start' does not contain start values for: %s", 
                                                 missing.intercepts))
  
  
  # order arguments as expected in the C++ code (1=asymptote, 2=rate, 3=intercept)
  # NOTE/TODO: The reason this is hardcoded is that the start parameters are the design matrix are created in the
  #            R code, and the coefficient order needs to be aligned properly with what the C++ code expects.
  # TODO: Fix this problem by rewriting CCoefConstraints::Unconstrain() to use names for input (instead of indices), 
  #       and to output the vector in the right order (i.e., in agreement with the names in the design matrix and in 
  #       the constraint matrix).
  contrasts <- contrasts[c('asymptote', 'invrate', 'intercept')]
    
  # rearange start parameters, make sure the coefficients for every parameter are contiguous
  start.names <- c()
  for(cur.name in names(contrasts)) {
    start.names <- c(start.names, cur.name)
    if(length(contrasts[[cur.name]]))
      start.names <- c(start.names, paste(cur.name, contrasts[[cur.name]], sep='.'))  
  }
  start <- start[ start.names ]
  names(start) <- start.names
  
  # add a corr.mrsat parameter if necessary
  if(!is.null(trial.id))
    start <- c(start, corr.mrsat=NA)

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

  # make sure the bounds on 'corr.mrsat' are reasonable, given the approximation to the multivariate normal CDF we use
  # Albers and Kallenberg recommend [1/sqrt(2), 1] for their approximation
  corr.mrsat.lower <- 1/sqrt(2)
  corr.mrsat.upper <- 1
  if('corr.mrsat' %in% rownames(constraint.matrix)) {
    if( constraint.matrix['corr.mrsat', 'lower'] < corr.mrsat.lower )
      constraint.matrix['corr.mrsat', 'lower'] <- corr.mrsat.lower
    if( constraint.matrix['corr.mrsat', 'upper'] > corr.mrsat.upper )
      constraint.matrix['corr.mrsat', 'upper'] <- corr.mrsat.upper
    
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

  return(list(constraints=constraint.matrix, start=start, contrasts=contrasts));
}


# translate parameters from formula notation to string notation,
# and check that all columns actually exist
# TODO: This function must be rather slow, optimize.

translate.parameters  <- function(data, dv, start, contrasts, constraints, 
                                  bias, signal, time, trial.id, summarize) {
  
  cnames <- list()
  
  # translate parameter: contrasts
  for(name in names(contrasts)) 
    contrasts[[name]] = check.formula.for.colname(data, paste('contrasts', name, sep='$'), contrasts[[name]], n.cols=NA)

  # initialize start parameters and create constraint matrix  
  satf.params <- init.satf.params(start=start, contrasts=contrasts, constraints=constraints, trial.id=trial.id)
  
  # translate parameter: dv
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

  # translate parameter: bias
  reportifnot(all(c('bias', 'poly.degree') %in% names(bias)), "Parameter 'bias' needs to containt 'bias' and 'poly.degree'.")
  bias[['bias']] = check.formula.for.colname(data, 'bias$bias', bias[['bias']], n.cols=NA)
  
  # check column names: signal, time, trial.id
  cnames$signal <- check.formula.for.colname(data, 'signal', signal)
  stopifnot.binary( data[,cnames$signal] )
    
  cnames$time <- check.formula.for.colname(data, 'time', time)

  if(!is.null(trial.id))
    cnames$trial.id <- check.formula.for.colname(data, 'trial.id', trial.id)
  
  # return
  list(dv=dv, start=satf.params$start, contrasts=satf.params$contrasts, bias=bias,
       constraints=satf.params$constraints, cnames=unlist(cnames) )
}
