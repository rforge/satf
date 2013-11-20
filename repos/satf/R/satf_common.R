
# TODO: Use as.quoted() for formulas.

check.formula.for.colname <- function(data, name, formula, n.cols=1) {
  reportifnot( is.formula(formula), sprintf("Parameter '%s' needs to be a formula.", name) )
  colname <- formula.terms(formula)
  if(!is.na(n.cols))
    reportifnot( length(colname) == n.cols, sprintf("Parameter '%s' needs specify exactly %d column.", name, n.cols) )
  if(length(colname))
    check.columns( colname, data )
  colname
}

get.column <- function(data, name, formula, as.type) {
  colname = check.formula.for.colname(data, name, formula)
  as.type( data[, colname ] )
}

satf.summarize.n.yes <- function(data, ids, signal='signal',
                                 time='time', dv='response') {
  stopifnot.binary( data[[ dv[1] ]] )
  stopifnot(!(signal %in% ids))
  data <- ddply(data, c(ids, signal), function(d) {
    d$n.responses.yes <- sum(d[[ dv[1] ]])
    d$n.responses <- length(d[[ dv[1] ]])
    for(colname in colnames(d)) {
      if(length(unique(d[[colname]])) > 1) {
        d[[colname]] <- mean(d[[colname]])
      }
    }
    d[1,]
  })
  data[[ dv[1] ]] <- NULL
  ddply(data, c(ids), function(d) {
    d[[ time ]] <- mean(d[[ time ]])
    d
  })
}

satf.summarize.dprime <- function(data, ids, signal, dv=NULL) {
  if(is.null(dv))
    dv <- c('n.responses.yes', 'n.responses')
  
  res <- ddply(data, ids, function(d) {
    d.noise  <- d[d[[signal]] == 0, ]
    d.signal <- d[d[[signal]] != 0, ]
    reportifnot(nrow(d.signal) == 1, sprintf("ncol(d.signal)=%d", nrow(d.signal)))
    reportifnot(nrow(d.noise) == 1, sprintf("ncol(d.noise)=%d", nrow(d.noise)))
    
    hits = d.signal[[ dv[1] ]]
    misses = d.signal[[ dv[2] ]] - hits
    fas = d.noise[[ dv[1] ]]
    crs = d.noise[[ dv[2] ]] - fas
    d.signal[[ dv[1] ]] <- NULL
    d.signal[[ dv[2] ]] <- NULL
    
    d.signal$time <- (d.signal$time+d.noise$time)/2
    
    summary <- ldply(1:nrow(d.signal), function(i) {
      dprime.var(hits=hits[i], misses=misses[i], fas=fas, crs=crs)
    })
    cbind(d.signal, summary)
  })
  res
}

initialize.compute.logLik  <- function(dv, contrasts, constraints, data, predicted.criterion, cnames)
{
  rcpp_initialize_logLikFn(dv, contrasts, constraints, data, predicted.criterion, cnames)
}

deinitialize.compute.logLik <- function() {
  rcpp_deinitialize_logLikFn( )
}

compute.logLik <- function(params)
{
  rcpp_compute_logLikFn(params)
}

transform.coefs <- function( coefs ) {
  res <- rcpp_transform_coefs( coefs )
  names(res) <- names( coefs )
  res
}

untransform.coefs <- function( coefs ) {
  res <- rcpp_untransform_coefs( coefs )
  names(res) <- names( coefs )
  res
}