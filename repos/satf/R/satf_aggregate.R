
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


