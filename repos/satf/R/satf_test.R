

aggregate.dprime <- function(data, flat.min=.1, include.last=TRUE) {
  cnds <- unique(data$condition)
  lvls <- lapply(cnds, function(x) sapply(0:1, function(y) paste(x,y) ))
  lvls <- unlist(lvls)
  data$responseType <- factor(paste(data$condition, data$response), 
                              levels=lvls)
  res <- ddply(data, .(time), function(d) {
    x <- xtabs(~d$responseType)
    a <- data.frame(hits=x['signal1 1'], fas=x['noise 1'],
                    misses=x['signal1 0'], crs=x['noise 0'], condition="signal1")
    b <- with(a, dprime.var(hits=hits, fas=fas, misses=misses, crs=crs, flat.min=flat.min))
    res <- cbind(a, t(b))
    if("signal2" %in% d$condition) {
      a <- data.frame(hits=x['signal2 1'], fas=x['noise 1'],
                      misses=x['signal2 0'], crs=x['noise 0'], condition="signal2")
      b <- with(a, dprime.var(hits=hits, fas=fas, misses=misses, crs=crs, flat.min=flat.min))
      res2 <- cbind(a, t(b))
      res <- rbind(res, res2)
    }
    res
  })
  res <- res[order(res$condition,res$time),]
  res <- ddply(res, .(condition), function(d) {
    d$last.hits <- c(NA, d$hits[-nrow(d)])
    d$last.fas <- c(NA, d$fas[-nrow(d)])
    d$last.misses <- c(NA, d$misses[-nrow(d)])
    d$last.crs <- c(NA, d$crs[-nrow(d)])
    d
  })
  res
}

generate.sat <- function(fn.satf1, fn.satf2, fn.bias, time, n.signal, n.noise, p.repeat=0)
{
  # fn.bias specifies the bias between the noise and signal1 distributions.
  #  The location of the criterion is computed from this bias. 
  bias <- fn.bias(time)
  dprime1 <- fn.satf1(time)
  FAlarms = pnorm(-dprime1/2 - bias)
  predicted.criterion <- qnorm(1-FAlarms)
  intervals <- 1:length(time)
  
  all.intervals <- rep(intervals, n.noise)
  all.time <- rep(time, n.noise)
  all.trial.id <- rep(1:n.noise, each=length(time))
  last.trial.id <- tail(all.trial.id,1)
  all.responses <- rbinom(length(time)*n.noise, size=1, prob=rep(FAlarms, n.noise))
  noise <- data.frame(condition="noise", signal=0, interval=all.intervals,
                      time=all.time, trial.id=all.trial.id, 
                      response=all.responses)
    
  all.intervals <- rep(intervals, n.signal)
  all.time <- rep(time, n.signal)
  all.trial.id <- rep(last.trial.id+1:n.signal, each=length(time))
  last.trial.id <- tail(all.trial.id,1)
  Hits1 <- 1-pnorm(predicted.criterion, mean=dprime1)
  all.responses <- rbinom(length(time)*n.signal, size=1, prob=rep(Hits1, n.signal))
  signal1 <- data.frame(condition="signal1", signal=1, interval=all.intervals,
                      time=all.time, trial.id=all.trial.id, response=all.responses)
  data <- rbind(noise, signal1)
  
  if(!is.null(fn.satf2)) {
    all.trial.id <- rep(last.trial.id+1:n.signal, each=length(time))
    dprime2 <- fn.satf2(time)
    Hits2 <- 1-pnorm(predicted.criterion, mean=dprime2)
    all.responses <- rbinom(length(time)*n.signal, size=1, prob=rep(Hits2, n.signal))
    signal2 <- data.frame(condition="signal2", signal=1, interval=all.intervals,
                          time=all.time, trial.id=all.trial.id, response=all.responses)
    data <- rbind(data, signal2)
  }
  
  data$last.response <- c(NA, data$response[-nrow(data)])
  data$last.response[data$interval==1] <- NA
  
  repeated <- as.logical( rbinom(nrow(data), size=1, prob=p.repeat) )
  repeated[is.na(data$last.response)] <- FALSE
  data$response[repeated] <- data$last.response[repeated]
  
  data
}



generate.mrsat2 <- function(fn.satf1, fn.satf2, fn.bias, time, 
                           n.signal, n.noise, processing.noise.sd=0)
{
  library(plyr)
  # fn.bias specifies the bias between the noise and signal1 distributions.
  #  The location of the criterion is computed from this bias. 
  bias <- fn.bias(time)
  dprime <- fn.satf1(time)
  FAlarms = pnorm(-dprime/2 - bias)
  criterion <- qnorm(1-FAlarms)
  intervals <- 1:length(time)
  n.dpoints <- length(time)

  generate.noise <- function(condition, signal, n.trials, id.start=1) {
   # generate noise
    all.trial.id <- rep(id.start+1:n.trials, each=n.dpoints)
    sds <- c(1, rep(processing.noise.sd, n.dpoints-1))
    tmp <- rnorm(n.trials*n.dpoints, sd=rep(sds, n.trials))
    tmp <- data.frame(trial.id=all.trial.id, tmp)
    data <- ddply(tmp, .(trial.id), function(d) {
      d$position <- cumsum(d$tmp)
      d
    })
    
    # Correct for the fact that adding repeated noise increases the SD
    # of the distribution, and thus decreases the distance between the signal
    # and noise distributions in terms of SDs.
    noise.vars <- (intervals-min(intervals))*(processing.noise.sd^2)
    actual.sds <- sqrt(1+noise.vars)
    data$position <- data$position/actual.sds
    
    data.frame(data, criterion=criterion, interval=intervals,
               time=time, condition=condition, signal=signal)  
  }
  
  # CONDITION 1: noise
  noise <- generate.noise(condition="noise", signal=0, n.trials=n.noise, id.start=1)

  # CONDITION 2: signal1
  signal1 <- generate.noise(condition="signal1", signal=1, n.trials=n.signal, id.start=n.noise+1)
  signal1$position <- signal1$position + dprime
  data <- rbind(noise, signal1)
  
  # CONDITION 3: signal2
  if(!is.null(fn.satf2)) {
    dprime <- fn.satf2(time)
    signal2 <- generate.noise(condition="signal2", signal=1, n.trials=n.signal, id.start=n.noise+n.signal+1)
    signal2$position <- signal2$position + dprime
    data <- rbind(data, signal2)
  }
  
  data$response <- as.integer(with(data, criterion < position))
  data$last.response <- c(NA, data$response[-nrow(data)])
  data$last.response[data$interval==1] <- NA
  data
}



library(inline)
.sig_sum_trials <- signature(nt="integer", tlen="integer", x="double")
.code_sum_trials <- "
      for (int i=0; i < (*nt)*(*tlen); i++) {
        if (i % *tlen != 0)
         x[i] = x[i]+x[i-1];
      }"
.sum_trials <- cfunction( .sig_sum_trials, .code_sum_trials, convention=".C")
sum_trials <- function(n.trials, trial.len, seq) {
  stopifnot(n.trials*trial.len == length(seq))
  .sum_trials(nt=n.trials, tlen=trial.len, x=seq)$x
}


generate.mrsat <- function(fn.satf1, fn.satf2, fn.bias, time, 
                           n.signal, n.noise, processing.noise.sd=0)
{
  library(plyr)
  # fn.bias specifies the bias between the noise and signal1 distributions.
  #  The location of the criterion is computed from this bias. 
  bias <- fn.bias(time)
  dprime <- fn.satf1(time)
  FAlarms = pnorm(-dprime/2 - bias)
  criterion <- qnorm(1-FAlarms)
  intervals <- 1:length(time)
  n.dpoints <- length(time)
  
  generate.noise <- function(condition, signal, n.trials, id.start=1) {
    all.trial.id <- rep(id.start+1:n.trials, each=n.dpoints)
    sds <- c(1, rep(processing.noise.sd, n.dpoints-1))
    tmp <- rnorm(n.trials*n.dpoints, sd=rep(sds, n.trials))
    position <- sum_trials(n.trials=n.trials, trial.len=n.dpoints, seq=tmp)
    
    # Correct for the fact that adding repeated noise increases the SD
    # of the distribution, and thus decreases the distance between the signal
    # and noise distributions in terms of SDs.
    noise.vars <- (intervals-min(intervals))*(processing.noise.sd^2)
    actual.sds <- sqrt(1+noise.vars)
    position <- position/actual.sds
    
    data.frame(trial.id=all.trial.id, position, criterion=criterion,
               interval=intervals, time=time, condition=condition, 
               signal=signal)
  }
  
  # CONDITION 1: noise
  noise <- generate.noise(condition="noise", signal=0, n.trials=n.noise, id.start=1)
  
  # CONDITION 2: signal1
  signal1 <- generate.noise(condition="signal1", signal=1, n.trials=n.signal, id.start=n.noise+1)
  signal1$position <- signal1$position + dprime
  data <- rbind(noise, signal1)
  
  # CONDITION 3: signal2
  if(!is.null(fn.satf2)) {
    dprime <- fn.satf2(time)
    signal2 <- generate.noise(condition="signal2", signal=1, n.trials=n.signal, id.start=n.noise+n.signal+1)
    signal2$position <- signal2$position + dprime
    data <- rbind(data, signal2)
  }
  
  data$response <- as.integer(with(data, criterion < position))
  data$last.response <- c(NA, data$response[-nrow(data)])
  data$last.response[data$interval==1] <- NA
  data
}

if(FALSE) {
fn.bias <- function(x) -10/(x+2)
fn.satf1 <- function(t) SATF(t, lambda=2, beta=1, delta=.4)
sim.n <- 10^3
data <- generate.mrsat(fn.satf1, fn.satf1, fn.bias, seq(-.5,6,.5),
                    n.noise=sim.n, n.signal=sim.n, processing.noise.sd=0)
data$responseGrammatical <- data$response
data$lastResponseGrammatical <- data$last.response
x.d <- aggregate.dprime(data)
ggplot(x.d, aes(x=time, y=dprime, color=condition))+geom_point()

source("~/Code/satf/R/misc.R")
source("~/Code/satf/R/satf_rawbinomial.R")
source("~/Code/satf/R/satf.R")
(p <- ssatf.binom(dv=list(~responseGrammatical,~lastResponseGrammatical), 
                  signal=~signal, 
                  fit.corr=TRUE, start=c(lambda=2, beta=1, delta=.4),
                  lambda=~1, beta=~1, delta=~1, 
                  data=data, pslope.delta=TRUE, plot=T))

fn <- function(t) logodds2p(RepeatProbPoly(t, p$par$repeat.params))
plot(fn, xlim=c(-1,6))

tail(data)

start.parameters <- function(data.aggregated) {
  time.max.proportion <- .95 # Assume that 95% of the asymptote is reached 
                             # at the last available point in time.
  intercept.proportion <- .1 # Use the earliest time at which 10% of the estimated
                             # asymptote have been reached as the y intercept.
  time.max <- max(data.aggregated$time)
  asymptote <- mean(subset(data.aggregated,time==time.max)$dprime)/time.max.proportion
  intercept.t <- min(with(data.aggregated, time[dprime/asymptote>intercept.proportion]))
  if(is.na(intercept.t)) intercept.t <- 0
  rate = -log(1-time.max.proportion)/(time.max-intercept.t)
  pmax(c(lambda=asymptote, beta=rate, delta=intercept.t), 0.1)
}

run.simulation <- function(n.trials, lambda, beta, delta, time, p.repeat=0) {
  satf1 <- function(t) SATF(t, lambda=lambda[1], beta=beta[1], delta=delta[1])
  satf2 <- function(t) SATF(t, lambda=lambda[2], beta=beta[2], delta=delta[2])
  bias1 <- function(t) -1/exp(t)

  tst <- generate(satf1, satf2, bias1, time, n.noise=n.trials, n.signal=n.trials,
                  p.repeat=p.repeat)
  nmap(colnames(tst), c('response'='responseGrammatical', 'signal'='stimulusGrammatical',
                      'last.response'='lastResponseGrammatical')) -> colnames(tst)
  tst$signaltype <- nmap(tst$condition, c('noise'=0,'signal1'=0,'signal2'=1), as.integer)

  tst.aggregated <- aggregate.dprime(tst)
  tst.aggregated$signaltype <- nmap(asc(tst.aggregated$condition), c('noise'=0,'signal1'=0,'signal2'=1), as.integer)
  par <- start.parameters(tst.aggregated)  
  start.params <- c('lambda'=par[['lambda']], 'beta'=par[['beta']], 'delta'=par[['delta']])
  start.params <- Fit.SATF.Aggregate(log(start.params), metric="adjR2", fixed=c(),
                      data=tst.aggregated, lambda=~1, beta=~1, delta=~1)
  start.params <- start.params$par

  adjR2.111 <- Fit.SATF.Aggregate(start.params, metric="adjR2", fixed=c(),
                      data=tst.aggregated,
                      lambda=~1, beta=~1, delta=~1)
  adjR2.112 <- Fit.SATF.Aggregate(c(start.params, 0), metric="adjR2", fixed=c(),
                      data=tst.aggregated,
                      lambda=~1, beta=~1, delta=~1+signaltype)
  adjR2.121 <- Fit.SATF.Aggregate(c(start.params,0), metric="adjR2", fixed=c(),
                      data=tst.aggregated,
                      lambda=~1, beta=~1+signaltype, delta=~1)
  adjR2.122 <- Fit.SATF.Aggregate(c(start.params,0,0), metric="adjR2", fixed=c(),
                                  data=tst.aggregated,
                      lambda=~1, beta=~1+signaltype, delta=~1+signaltype)
  adjR2.211 <- Fit.SATF.Aggregate(c(start.params,0), metric="adjR2", fixed=c(),
                                  data=tst.aggregated,
                      lambda=~1+signaltype, beta=~1, delta=~1)
  adjR2.212 <- Fit.SATF.Aggregate(c(start.params,0,0), metric="adjR2", fixed=c(),
                                  data=tst.aggregated,
                      lambda=~1+signaltype, beta=~1, delta=~1+signaltype)
  adjR2.221 <- Fit.SATF.Aggregate(c(start.params,0,0), metric="adjR2", fixed=c(),
                                  data=tst.aggregated,
                      lambda=~1+signaltype, beta=~1+signaltype, delta=~1)
  adjR2.222 <- Fit.SATF.Aggregate(c(start.params,0,0,0), metric="adjR2", fixed=c(),
                                  data=tst.aggregated,
                      lambda=~1+signaltype, beta=~1+signaltype, delta=~1+signaltype)

  criterion.fn <- FALSE
  fixed <- c( p.repeat=p2logodds(.1))
  start.params <- c(start.params)
  lik.111 <- Fit.SATF.RawBinomial(start.params, fixed=fixed, data=tst,
                                  criterion.fn=criterion.fn,
                      lambda=~1, beta=~1, delta=~1)
  lik.112 <- Fit.SATF.RawBinomial(c(start.params,0), fixed=fixed, data=tst, 
                                  criterion.fn=criterion.fn,
                      lambda=~1, beta=~1, delta=~1+signaltype)
  lik.121 <- Fit.SATF.RawBinomial(c(start.params,0), fixed=fixed, data=tst, 
                                  criterion.fn=criterion.fn,
                      lambda=~1, beta=~1+signaltype, delta=~1)
  lik.122 <- Fit.SATF.RawBinomial(c(start.params,0,0), fixed=fixed, data=tst, 
                                  criterion.fn=criterion.fn,
                      lambda=~1, beta=~1+signaltype, delta=~1+signaltype)  
  lik.211 <- Fit.SATF.RawBinomial(c(start.params,0), fixed=fixed,  data=tst, 
                                  criterion.fn=criterion.fn,
                      lambda=~1+signaltype, beta=~1, delta=~1)
  lik.212 <- Fit.SATF.RawBinomial(c(start.params,0,0), fixed=fixed, data=tst,
                                  criterion.fn=criterion.fn,
                      lambda=~1+signaltype, beta=~1, delta=~1+signaltype)
  lik.221 <- Fit.SATF.RawBinomial(c(start.params,0,0), fixed=fixed,  data=tst,
                                  criterion.fn=criterion.fn,
                      lambda=~1+signaltype, beta=~1+signaltype, delta=~1)
  lik.222 <- Fit.SATF.RawBinomial(c(start.params,0,0,0), fixed=fixed, data=tst,
                                  criterion.fn=criterion.fn,
                      lambda=~1+signaltype, beta=~1+signaltype, delta=~1+signaltype)

  s <- function(x) c(x$value, length(x$par))
  adjR2 <- c(adjR2.111=s(adjR2.111), adjR2.112=s(adjR2.112), adjR2.121=s(adjR2.121),
             adjR2.122=s(adjR2.122), adjR2.211=s(adjR2.211), adjR2.212=s(adjR2.212),
             adjR2.221=s(adjR2.221), adjR2.222=s(adjR2.222))
  lik <- c(lik.111=s(lik.111$res), lik.112=s(lik.112$res), lik.121=s(lik.121$res),
           lik.122=s(lik.122$res), lik.211=s(lik.211$res), lik.212=s(lik.212$res),
           lik.221=s(lik.221$res), lik.222=s(lik.222$res))
  data.frame(t(c(adjR2, lik)))
}

with(cur.params, run.simulation(n.trials=30, lambda=c(2,2), beta=c(.5,.5),
                                delta=c(1,2),  time=seq(0,6,.4)))

library(reshape)

parameters <- ldply(c(3,2.8,2.6), function(lambda) {
  ldply(c(1,.8,.6), function(beta) {
    ldply(c(1,.8,.6), function(delta) {
      ldply(c(0,.25,.5), function(p.repeat) {
        c(lambda=lambda, beta=beta, delta=delta, p.repeat=p.repeat)
      })
    })
  })
})
parameters


sim.n <- 100
res <- NULL
for(i in 1:nrow(parameters)) {
  print(date())
  cur.params <- parameters[i,]
  print(cur.params)
  cur.res <- ldply(1:sim.n, function(j) {
    #print(date())
    cur.res <- with(cur.params, run.simulation(n.trials=30, time=seq(0,6,.4),
                    lambda=c(3,lambda), beta=c(1,beta), delta=c(1,delta),
                                               p.repeat=p.repeat))
    #print(date())
    cur.res$i <- j + 2*sim.n*(i-1)
    cur.res
  }, .progress = "text")
  cur.res$lambda <- cur.params$lambda
  cur.res$beta <- cur.params$beta
  cur.res$delta <- cur.params$delta
  cur.res$p.repeat <- cur.params$p.repeat
  res <- rbind(res, cur.res)
  save(res, file="./ModelEstimates.rda")
}

run.simulation(n.trials=30,  time=seq(0,6,.4),
                    lambda=c(3.5, lambda), beta=c(.5, beta), delta=c(1,delta))

dprime.var

date()
x <- run.simulation(n.trials=30, lambda=c(3.5,3.5), beta=c(1,1),
                           delta=c(1,1), time=seq(0,6,.4))
date()

}
