
source("./satf_load.R")

  fn.bias <- function(x) 0.1*x #-1/(x+2)
  fn.satf <- function(t) SATF(t, lambda=3, beta=1, delta=.4)
  sim.n <- 10^3

  time = seq(-.5,5,.5)
  data <- generate.sat(criterion=fn.bias, dprime=fn.satf, time=time, n=sim.n)

  data.nyes <- satf.summarize.n.yes(data, c('interval','condition'), dv='response')
  data.nyes

  (p1 <- satf(dv=c(response=~response), signal=~signal,
           start     = c(lambda=2.1, beta=1, delta=.4),
           contrasts = c(lambda=~1, beta=~1, delta=~1), 
           constraints=list(), time=~time, 
           bias=list(bias=~1, poly.degree=10), # trial.id=~trial.id,
           data=data, metric="logLikRaw"))


  (p2 <- satf(dv=c(n.responses.yes=~n.responses.yes, n.responses=~n.responses), 
           signal=~signal,
           start     = c(lambda=2.1, beta=1, delta=.4),
           contrasts = c(lambda=~1, beta=~1, delta=~1), 
           constraints=list(), time=~time, 
           bias=list(bias=~1, poly.degree=5), # trial.id=~trial.id,
           data=data.nyes, metric="logLik"))


x <- ddply(subset(data, subset=signal==0), .(interval), function(d) {c(time=mean(d$time), response=mean(d$response))})
with(x, plot(time, response))

data$time <- jitter(data$time)
with(, plot(response ~ time) )

(m <- lm(response ~ poly(time,5) , data, subset=signal==0))
plot(predict(m, data.frame(time=seq(0,6,.1))))

terms(~low(.0)+high)

x.d <- aggregate.dprime(data)
ggplot(x.d, aes(x=time, y=dprime, color=condition))+geom_point()

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
