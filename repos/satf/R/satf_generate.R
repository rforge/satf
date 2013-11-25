
process.arg <- function(arg, time) {
  if(typeof(arg) == "closure") arg = arg(time)
  if(length(arg) == 1 ) arg = rep(arg, length(time))
  stopifnot(length(arg) == length(time) )
  arg
}

generate.sat.condition.old <- function(criterion, dprime, time, trial.ids, label)
{
  intervals <- 1:length(time)
  # fn.bias specifies the bias between the noise and signal1 distributions.
  #  The location of the criterion is computed from this bias. 
  criterion <- process.arg(criterion, time)
  dprime <- process.arg(dprime, time)

  pYes = 1-pnorm(criterion, mean=dprime)
  pYes = rep(pYes, length(trial.ids))
  response = rbinom(length(pYes), size=1, prob=pYes)
  data.frame(condition=label, interval=intervals, time=time, 
             trial.id=rep(trial.ids, each=length(time)), response=response)
}

generate.sat.condition <- function(criterion, dprime, time, trial.ids, label, rho=0)
{
  intervals <- 1:length(time)
  # fn.bias specifies the bias between the noise and signal1 distributions.
  #  The location of the criterion is computed from this bias. 
  criterion <- process.arg(criterion, time)
  dprime <- process.arg(dprime, time)
  data <-   data.frame( condition=label, interval=intervals, time=time, 
                        trial.id=rep(trial.ids, each=length(time)), 
                        dprime=dprime, criterion=criterion )
  data$noise <- rnorm( nrow(data) )
  data$noise = rcpp_correlate(data$trial.id, data$noise, rho)
  data$dprime.cur <- data$dprime + data$noise
  data$response <- data$dprime.cur > data$criterion
  data
}

generate.sat <- function(criterion, dprime, time, n, rho=0) {
  data0 <- generate.sat.condition(criterion=criterion, dprime=0, time=time, 
                                  trial.ids=1:n, label="condition1", rho=rho)
  data0$signal <- 0
  data1 <- generate.sat.condition(criterion=criterion, dprime=dprime,  time=time, 
                                  trial.ids=1:n+n, label="condition1", rho=rho)
  data1$signal <- 1
  rbind(data0, data1)
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
