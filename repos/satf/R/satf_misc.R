
stopifnot.binary <- function(x) stopifnot( all(x %in% c(0,1)) )
stopifnot.binary.na <- function(x) stopifnot( all(x %in% c(0,1,NA)) )


reportifnot <- function(cond, text) {
  if(!cond)
    stop(text)
}

remove.outliers <- function(x, factor) {
  distance <- abs(mean(x)-x);
  x <- x[order(distance, decreasing=T)];
  distance <- sort(distance, decreasing=T);
  exclude <- c()
  for(i in 2:length(x)) {
    if(distance[i-1]/distance[i] > factor) {
      exclude <- c(exclude, i-1)
    } else {
      break;
    }
  }
  if(length(exclude) > length(x)) {
    stop("More than half the data were excluded as outliers.")
  }
  if(length(exclude)) {
    print(paste("Excluding:", paste(x[exclude], collapse="")))
    x <- x[(-1*exclude)]
  }
  x
}

check.columns <- function(columns, data) {
  present <-  columns %in% colnames(data)
  if(!all(present)) {
    missing <- columns[!present]
    missing <- paste("'",missing,"'",sep="", collapse=", ")
    msg <- paste("Missing columns in data:", missing)
    stop(msg)
  }
}

row.sums <- function(mat) {
  if(is.vector(mat))
    return(mat)
  else
    return(rowSums(mat))
}

col.sums <- function(x, len) {
  if(is.vector(x))
    return(x)
  dims = 1L
  dn <- dim(x)
  if (dims < 1L || dims > length(dn) - 1L) 
    stop("invalid 'dims'")
  n <- prod(dn[1L:dims])
  dn <- dn[-(1L:dims)]
  z <- .Internal(colSums(x, n, prod(dn), na.rm=FALSE))
  if (length(dn) > 1L) {
    dim(z) <- dn
    dimnames(z) <- dimnames(x)[-(1L:dims)]
  }
  else names(z) <- dimnames(x)[[dims + 1]]
  z
}

#col.sums <- function(mat, is.vec=is.vector(mat)) {
#  if(is.vec)
#    return(mat)
#  dn <- dim(mat)
#  n <- prod(dn[1L])
#  dn <- dn[-1L]
#  .Internal(colSums(mat, n, prod(dn), na.rm=F))
#}

default <- function(lst, name, val) {
  if(!name %in% names(lst)) lst[[name]] <- val; 
  lst
}


optim.to.precision <- function(start, optim.digits, control, ...)
{
  method = "Nelder-Mead"
  if("method" %in% names(control))
    method = control['method']
  run.optim <- function(start) optim(par=start, method=method, control=control, ...)
  res <- run.optim(start)
  if(method == "SANN")
    return(res)
  old.cnt <- NULL
  while(TRUE) {
    old.value <- res$value
    old.cnt <- paste(old.cnt, res$counts[['function']], sep=' ')
    res <- run.optim(res$par)
    if(is.na(optim.digits))
      return(res)
    if(round(old.value, optim.digits) == round(res$value, optim.digits)) {
      res$counts[['function']]  <- paste(old.cnt, res$counts[['function']], sep=' ')
      return(res)
    }
  }
}


Poly10 <- function(t, d) with(as.data.frame(t(d)), c1+c2*t+c3*t^2+c4*t^3+c5*t^4+c6*t^5+c7*t^6+c8*t^7+c9*t^8+c10*t^9)
Poly9 <- function(t, d) with(as.data.frame(t(d)), c1+c2*t+c3*t^2+c4*t^3+c5*t^4+c6*t^5+c7*t^6+c8*t^7+c9*t^8)
Poly8 <- function(t, d) with(as.data.frame(t(d)), c1+c2*t+c3*t^2+c4*t^3+c5*t^4+c6*t^5+c7*t^6+c8*t^7)
Poly7 <- function(t, d) with(as.data.frame(t(d)), c1+c2*t+c3*t^2+c4*t^3+c5*t^4+c6*t^5+c7*t^6)
Poly6 <- function(t, d) with(as.data.frame(t(d)), c1+c2*t+c3*t^2+c4*t^3+c5*t^4+c6*t^5)
Poly5 <- function(t, d) with(as.data.frame(t(d)), c1+c2*t+c3*t^2+c4*t^3+c5*t^4)
Poly4 <- function(t, d) with(as.data.frame(t(d)), c1+c2*t+c3*t^2+c4*t^3)
Poly3 <- function(t, d) with(as.data.frame(t(d)), c1+c2*t+c3*t^2)
Poly2 <- function(t, d) with(as.data.frame(t(d)), c1+c2*t)
Poly1 <- function(t, d) with(as.data.frame(t(d)), c1+c2*t)
BiasFunctions <- list(Poly1, Poly2, Poly3, Poly4, Poly5,
                      Poly6, Poly7, Poly8, Poly9, Poly10)


RepeatProbPoly <- function(t, d) with(as.data.frame(t(d)), r1 + r2*t )

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


formula.terms <- function(formula) { 
  terms.formula <- terms(formula)
  stopifnot(attr(terms.formula, "intercept")==1)
  attr(terms.formula, "term.labels")
}
n.terms <- function(formula) length(formula.terms(formula))

formula.terms <- function(formula) { 
  terms.formula <- terms(formula)
  stopifnot(attr(terms.formula, "intercept")==1)
  attr(terms.formula, "term.labels")
}
n.terms <- function(formula) length(formula.terms(formula))


dprime.var <- function(hits, fas, misses, crs, flat.min=.5)
{
  # 'correct' possible probabilites of 1  
  n.signal <- hits+misses; n.noise <- fas+crs
  if(n.signal < n.noise) {
    hits <- hits+flat.min; misses <- misses+flat.min;
    fas <- fas+flat.min*(n.noise/n.signal); crs <- crs+flat.min*(n.noise/n.signal);
  } else {
    hits <- hits+flat.min*(n.signal/n.noise); misses <- misses+flat.min*(n.noise/n.signal);
    fas <- fas+flat.min; crs <- crs+flat.min;
  }
  n.signal <- hits+misses; n.noise <- fas+crs

  p.hit <- hits/n.signal
  p.CR <- crs/n.noise
  dprime <- qnorm(p.hit)-qnorm(1-p.CR)
  c    <- -0.5*(qnorm(p.hit)+qnorm(1-p.CR))
  dprime.var <- ( p.hit*(1-p.hit) ) / ( n.signal*dnorm(qnorm(p.hit))^2 ) +
                ( p.CR*(1-p.CR) ) / ( n.noise*dnorm(qnorm(p.CR))^2 )
  
  c(dprime=dprime, dprime.var=dprime.var, c=c, c.var=dprime.var*0.25)
}


# This version is too inefficient: It requires 23 sec for 10^6 calls, where
# the new SATF (below) requires only 3 sec.
#SATF <- function(t, lambda, beta, delta)
#  ifelse(t >= delta, lambda*(1-exp(-beta*(t-delta))), 0)

SATF <- function(t, lambda, beta, delta)
  (t >= delta)*(lambda*(1-exp(-1/beta*(t-delta)))) 

