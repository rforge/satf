
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

init_designmatrix <- function(data, contrasts, bias, cnames, coreparams.satf, coreparams.bias) {

  reportifnot(all(coreparams.satf %in% names(contrasts)), 
                  "Parameter 'contrasts' needs to contain 'asymptote', 'invrate', and 'intercept'.")
  if( !is.list(bias) ) {
    bias.contrast <- bias
    bias <- list()
    bias[coreparams.bias] <- bias.contrast
    
  } else if( !all(coreparams.bias %in% names(bias)) ) {
    stop("Parameter 'bias' needs to contain formulae for 'bias.min', 'bias.max', 'bias.invrate', and 'bias.intercept', or to be a formula.")
    
  }

  # rearange start parameters, make sure the coefficients for every parameter are contiguous, 
  # and create the design matrix
  dm.colnames <- c();
  contrasts.coef.names <- c(); 
  # deal with satf params
  for(cur.name in names(contrasts)) {
    contrasts.coef.names <- c(contrasts.coef.names, cur.name)
    dm.colnames <- c(dm.colnames, cnames['signal'])
    if(length(contrasts[[cur.name]])) {
      contrasts.coef.names <- c(contrasts.coef.names, paste(cur.name, contrasts[[cur.name]], sep='.'))  
      dm.colnames <- c(dm.colnames, contrasts[[cur.name]])
    }
  }
  dm.coef.cnt <- unlist(lapply(contrasts, length))+1

  # deal with bias params
  bias.coef.names <- c(); 
  for(cur.name in names(bias)) {
    bias.coef.names <- c(bias.coef.names, cur.name)
    dm.colnames <- c(dm.colnames, '1')
    if(length(bias[[cur.name]])){
      bias.coef.names <- c(bias.coef.names, paste(cur.name, bias[[cur.name]], sep='.'))  
      dm.colnames <- c(dm.colnames, bias[[cur.name]])
    }
  }
  dm.coef.cnt <- c(dm.coef.cnt, unlist(lapply(bias, length))+1)
  coef.names <- c(contrasts.coef.names, bias.coef.names)

  # create a column consisting of only 1's, for use with the bias
  data$'1' <- 1
  
  # create the design matrix
  dm <- data[,dm.colnames]
  colnames(dm) <- coef.names
  
  # set all satf coefficients to 0 for noise trials
  satf.ncol <- sum(dm.coef.cnt[coreparams.satf])
  dm[,1:satf.ncol] <- dm[,1:satf.ncol]*data[, cnames['signal'] ]
  
  # add a corr.mrsat parameter if necessary
  if( 'trial.id' %in% names(cnames) ) {
    coef.names <- c(coef.names, "corr.mrsat")
    dm.coef.cnt <- c(dm.coef.cnt, corr.mrsat=1)
    trial.id <- data[[ cnames[['trial.id']] ]]
    dm$corr.mrsat <-  as.integer(c(F, trial.id[-1] == trial.id[-length(trial.id)]))
  }
  
  list(dm=as.matrix(dm), dm.coef.cnt=dm.coef.cnt, 
    contrasts.coefs=contrasts.coef.names, bias.coefs=bias.coef.names)
}



init_coefs_and_constraints <- function(coefnames, start, constraints, coreparams)
{
  reportifnot(is.vector(start) && !is.list(start), "Parameter 'start' needs to be a vector.")
  reportifnot(!is.null(start) || !any(is.na(start)), "Parameter 'start' was not provided or containts NAs.")
  
  # make sure start values were provided for all the necessary parameters
  missing.coreparams = coreparams[!(coreparams %in% names(start))]
  missing.coreparams.str = paste(missing.coreparams, collapse=", ")
  reportifnot( length(missing.coreparams) == 0, sprintf("Parameter 'start' does not contain start values for: %s", 
                                                        missing.coreparams.str))

  # warn about unused start values
  unused.startvalues = names(start[!(names(start) %in% coefnames)])
  unused.startvalues.str = paste(unused.startvalues, collapse=", ")
  if( length(unused.startvalues) > 0) {
      warning(sprintf("Parameter 'start' contains unused start values for: %s", unused.startvalues.str))
  }

  start = start[coefnames]
  names(start) = coefnames
  
  # initialize constraints matrix
  constraint.matrix = matrix(c(-Inf, Inf), nrow=length(start), ncol=2, byrow=TRUE, 
                             dimnames=list(coefnames, c('lower', 'upper')))
  
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
  
  # reverse the constraint list to make user-specified constraints take precedence over defaults
  # (the reason there can be several specifications for one parameter is that constraints names are
  #  specified by regex expressions)
  constraints = rev(constraints)

  # fill the constraints matrix
  warn.about.start.params <- FALSE
  for(name in names(start))
  {
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

translate.parameters  <- function(data, dv, contrasts, bias, signal, time, trial.id) {
  
  cnames <- list()

  # translate parameter: contrasts
  for(name in names(contrasts)) 
    contrasts[[name]] = check.formula.for.colname(data, paste('contrasts', name, sep='$'), contrasts[[name]], n.cols=NA)

  # translate parameter: contrasts
  if(is.list(bias)) {
    for(name in names(bias)) 
      bias[[name]] = check.formula.for.colname(data, paste('bias', name, sep='$'), bias[[name]], n.cols=NA)
  }
  
  # TODO: Make sure that trials are not discontinuous.
  # check column names: signal, time, trial.id
  cnames$signal <- check.formula.for.colname(data, 'signal', signal)
  stopifnot.binary( data[,cnames$signal] )
  
  cnames$time <- check.formula.for.colname(data, 'time', time)
  
  if(!is.null(trial.id))
    cnames$trial.id <- check.formula.for.colname(data, 'trial.id', trial.id)
  
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

  # return
  list(dv=dv, contrasts=contrasts, bias=bias, cnames=unlist(cnames) )
}

set_start_defaults <- function(start, set.corr.mrsat=FALSE) {
  start = default(start, 'bias.min', -1)
  start = default(start, 'bias.max', 1)
  start = default(start, 'bias.invrate', 1)
  start = default(start, 'bias.intercept', 0)
  if(set.corr.mrsat)
    start = default(start, 'corr.mrsat', .85)
  start
}

set_constraints_defaults <- function(constraints) {
  constraints = default(constraints, 'asymptote', c(0,Inf))
  constraints = default(constraints, 'invrate', c(0,Inf))
  constraints = default(constraints, 'intercept', c(0,Inf))
  # Albers and Kallenberg recommend [1/sqrt(2), 1] for their approximation
  constraints = default(constraints, 'corr.mrsat', c(1/sqrt(2), 1))
  constraints
}
