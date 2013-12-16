

source("~/CodeSATF/test/satf_load.R")

  fn.bias <- function(t) t*0 # SATF(t, asymptote=1, invrate=1, intercept=.4)
  fn.satf1 <- function(t) SATF(t, asymptote=3, invrate=1, intercept=.4)
  fn.satf2 <- function(t) SATF(t, asymptote=2, invrate=1, intercept=.4)
  sim.n <- 10^3
  time = seq(-.5,5,.5)
  data1 <- satf_generate(criterion=fn.bias, dprime=fn.satf1, time=time, n=sim.n, rho=.72, label="condition1")
  data2 <- satf_generate(criterion=fn.bias, dprime=fn.satf2, time=time, n=sim.n, rho=.72, label="condition2")
  data <- rbind(data1, data2)

  data.nyes <- satf_aggregate_nyes(data, id=c('condition','interval'), time.id=c('interval'), signal='signal')
  data.dprime <- satf_aggregate_dprime(data.nyes, id=c('condition', 'interval'), signal='signal', 
                                       dv=c('n.responses.yes', 'n.responses'))

  # data.nyes <- data.nyes[with(data.nyes, order(condition,signal, interval)),]
  data$c1 <- ifelse(data$condition=="condition1", 1, 0)
  data$c2 <- ifelse(data$condition=="condition2", 1, 0)
  data.nyes$c1 <- ifelse(data.nyes$condition=="condition1", 1, 0)
  data.nyes$c2 <- ifelse(data.nyes$condition=="condition2", 1, 0)

#  with(data.nyes, plot(time, criterion))
#  with(data.dprime, plot(time, c + dprime/2, col=condition))


source("~/CodeSATF/test/satf_load.R")

#data.nyes.cur <- data.nyes
#print(data.nyes.cur)
#save(data, file="/tmp/xxx.rda")
(load("/tmp/xxx.rda"))


date()
(m.raw <- satf(dv=c(response=~response), signal=~signal,
               start     = c(asymptote=2, invrate=1, intercept=.4, bias.min=0),
               contrasts = c(asymptote=~1+c2, invrate=~1+c2, intercept=~1+c2),
               bias = ~1+c2, constraints=list(), time=~time, #trial.id=~trial.id,
               data=data, metric="logLikRaw", method="BFGS", debug=F, stepwise=T))
date()


stop()



(m.raw <- satf(dv=c(response=~response), signal=~signal,
               start     = c(asymptote=2, invrate=1, intercept=.4, bias.min=0),
               contrasts = c(asymptote=~1+c2, invrate=~1+c2, intercept=~1+c2),
               bias = ~1+c2, constraints=list(), time=~time, #trial.id=~trial.id,
               data=data, metric="logLikRaw"))
stop()


(m.nyes <- satf(dv=c(n.responses.yes=~n.responses.yes, n.responses=~n.responses), signal=~signal,
                start     = c(asymptote=2, invrate=1, intercept=.4, bias.min=0),
                contrasts = c(asymptote=~1, invrate=~1, intercept=~1),
                bias = ~1, constraints=list(), time=~time, data=data.nyes.cur, metric="logLik"))



stop()


(m.nyes <- satf(dv=c(n.responses.yes=~n.responses.yes, n.responses=~n.responses), signal=~signal,
                start     = c(asymptote=2, invrate=1, intercept=.4, bias.min=0),
                contrasts = c(asymptote=~1+c2, invrate=~1+c2, intercept=~1+c2),
                bias = ~1+c2, constraints=list(), time=~time, data=data.nyes[1:4,], metric="logLik"))


(m.raw <- satf(dv=c(response=~response), signal=~signal,
               start     = c(asymptote=2, invrate=1, intercept=.4, bias.min=0),
               contrasts = c(asymptote=~1+c2, invrate=~1+c2, intercept=~1+c2),
               bias = ~1+c2, constraints=list(), time=~time, #trial.id=~trial.id,
               data=data, metric="logLikRaw"))

(m.raw <- satf(dv=c(response=~response), signal=~signal,
                start     = c(asymptote=2, invrate=1, intercept=.4,
                              bias.min=-1, bias.max=1, bias.intercept=0, bias.invrate=1),
                contrasts = c(asymptote=~1+c2, invrate=~1+c2, intercept=~1+c2),
                bias = ~1+c2, constraints=list(), time=~time, #trial.id=~trial.id,
                data=data, metric="logLik"))
date()


source("./satf_load.R")

(p1.nyes1 <- satf(dv=c(n.responses.yes=~n.responses.yes, n.responses=~n.responses), signal=~signal,
                start     = c(asymptote=2, invrate=1, intercept=.4, bias.min=0),
                contrasts = c(asymptote=~1+c2, invrate=~1+c2, intercept=~1+c2),
                bias = ~1+c2, constraints=list(), time=~time, data=data.nyes, metric="logLik"))



(p1.nyes2 <- satf(dv=c(n.responses.yes=~n.responses.yes, n.responses=~n.responses), signal=~signal,
                  start     = c(asymptote=2, invrate=1, intercept=.4),
                  contrasts = c(asymptote=~1, invrate=~1, intercept=~1),
                  bias = c(bias.min=~1, bias.max=~1, bias.intercept=~1, bias.invrate=~1), 
                  constraints=list(bias.min=-1.2, bias.max=3, bias.intercept=-1.5, bias.invrate=3.7), time=~time,
                  data=data.nyes, metric="logLik"))


source("./satf_load.R")

date()
(p1.raw <- satf(dv=c(response=~response), signal=~signal,
                start     = c(asymptote=2, invrate=1, intercept=.4,
                              bias.min=0.35, bias.max=0.35, bias.intercept=-27, bias.invrate=291),
                contrasts = c(asymptote=~1, invrate=~1, intercept=~1),
                bias = c(bias.min=~1, bias.max=~1, bias.intercept=~1, bias.invrate=~1), 
                constraints=list(), time=~time, trial.id=~trial.id,
                data=data, metric="logLik"))
date()


x <- list(a=NULL, b=c('a', 'd'), c=NULL)

xx <- c()
for(i in 1:length(x)) {
  cur_name <- names(x)[i]
  xx <- c(xx, cur_name)
  if(length(x[[i]]))
    xx <- c(xx, paste(cur_name, x[[i]], sep='.'))  
}

  data.nyes <- satf.summarize.n.yes(data, c('interval','condition'), dv='response')
data.nyes$dprime <- NULL
  data.dprime <- satf.summarize.dprime(data.nyes, .(interval), signal='signal', dv=c('n.responses.yes', 'n.responses'))
data.dprime

source("./satf_load.R")

  data$lambda3 <- satf(dv=c(response=~response), signal=~signal,
           start     = c(lambda=3, beta=1, delta=.4),
           contrasts = c(lambda=~1, beta=~1, delta=~1), 
           constraints=list(lambda=3, beta=1, delta=.4, corr.mrsat=.99), time=~time, 
           bias=list(bias=~1, poly.degree=10), trial.id=~trial.id,
           data=data, metric="logLikRaw")

  data$lambda1 <- satf(dv=c(response=~response), signal=~signal,
                         start     = c(lambda=3, beta=1, delta=.4),
                         contrasts = c(lambda=~1, beta=~1, delta=~1), 
                         constraints=list(lambda=1, beta=1, delta=.4, corr.mrsat=.99), time=~time, 
                         bias=list(bias=~1, poly.degree=10), trial.id=~trial.id,
                         data=data, metric="logLikRaw")

  data$indep.lambda1 <- satf(dv=c(response=~response), signal=~signal,
                          start     = c(lambda=3, beta=1, delta=.4),
                          contrasts = c(lambda=~1, beta=~1, delta=~1), 
                          constraints=list(lambda=1, beta=1, delta=.4), time=~time, 
                          bias=list(bias=~1, poly.degree=10), #trial.id=~trial.id,
                          data=data, metric="logLikRaw")

  data$indep.lambda3 <- satf(dv=c(response=~response), signal=~signal,
                           start     = c(lambda=3, beta=1, delta=.4),
                           contrasts = c(lambda=~1, beta=~1, delta=~1), 
                           constraints=list(lambda=3, beta=1, delta=.4), time=~time, 
                           bias=list(bias=~1, poly.degree=10), #trial.id=~trial.id,
                           data=data, metric="logLikRaw")

sum(data$lambda1)
sum(data$lambda3)


rcpp_pnorm2d(-2.9, -2.9, .99, TRUE)/pnorm(-2.9)

rcpp_pnorm2d(-2.9, -2.9, .99, TRUE)/pnorm(-2.9)

data.tmp <- subset(data, lambda1 > lambda3)

which.max(data.tmp$lambda1-data.tmp$lambda3)

(data.tmp[2,])

(subset(data, trial.id==1657))[1:10,]

x1 <- with(data.tmp, tapply(paramsCorrect, list(time), sum))
x2 <- with(data.tmp, tapply(paramsIncorr, list(time), sum))
x1-x2

rbind(x1, x2)

source("./satf_load.R")

rcpp_pnorm2d(-0.04, -0.04, .9999, FALSE)

(p1.raw <- satf(dv=c(response=~response), signal=~signal,
                start     = c(lambda=3, beta=1, delta=.4),
                contrasts = c(lambda=~1, beta=~1, delta=~1), 
                constraints=list(lambda=3, beta=1, delta=.4, corr.mrsat=.8), time=~time, 
                bias=list(bias=~1, poly.degree=10), trial.id=~trial.id,
                data=data, metric="logLikRaw"))
[1] -2559

(p1.raw <- satf(dv=c(response=~response), signal=~signal,
                start     = c(lambda=3, beta=1, delta=.4),
                contrasts = c(lambda=~1, beta=~1, delta=~1), 
                constraints=list(), time=~time, 
                bias=list(bias=~1, poly.degree=10), trial.id=~trial.id,
                data=data, metric="logLikRaw"))

(p1.raw <- satf(dv=c(response=~response), signal=~signal,
                start     = c(lambda=3, beta=1, delta=.4),
                contrasts = c(lambda=~1, beta=~1, delta=~1), 
                constraints=list(), time=~time, 
                bias=list(bias=~1, poly.degree=10), #trial.id=~trial.id,
                data=data, metric="logLikRaw"))

  # fit to aggregated data
  (p1.agg <- satf(dv=c(n.responses.yes=~n.responses.yes, n.responses=~n.responses), 
           signal=~signal,
           start     = c(lambda=2.1, beta=1, delta=.4),
           contrasts = c(lambda=~1, beta=~1, delta=~1), 
           constraints=list(), time=~time, 
           bias=list(bias=~1, poly.degree=5), # trial.id=~trial.id,
           data=data.nyes, metric="logLik"))


  data.mr <- generate.mrsat(fn.satf1=fn.satf, fn.satf2=NULL, fn.bias=fn.bias, time=time,
                           n.signal=10^3, n.noise=10^3, processing.noise.sd=0.5)
  data.mr <- data.mr[,c('trial.id','interval','time','condition','signal','response')]
head(data.mr)
  data.mr.nyes <- satf.summarize.n.yes(data.mr, c('interval','condition'), dv='response')
head(data.mr.nyes)

  data.mr.dprime <- satf.summarize.dprime(data.mr.nyes, .(interval), signal='signal', 
                                       dv=c('n.responses.yes', 'n.responses'))
data.mr.dprime

source("./satf_load.R")

  # fit to raw data
  (p.mr.raw1 <- satf(dv=c(response=~response), signal=~signal,
           start     = c(lambda=2.1, beta=1, delta=.4),
           contrasts = c(lambda=~1, beta=~1, delta=~1), 
           constraints=list(corr.mrsat=.7), time=~time, 
           bias=list(bias=~1, poly.degree=10), trial.id=~trial.id,
           data=data.mr, metric="logLikRaw"))
[1] -16280

  # fit to raw data
  (p.mr.raw2 <- satf(dv=c(response=~response), signal=~signal,
                  start     = c(lambda=2.1, beta=1, delta=.4),
                  contrasts = c(lambda=~1, beta=~1, delta=~1), 
                  constraints=list(corr.mrsat=1), time=~time, 
                  bias=list(bias=~1, poly.degree=10), trial.id=~trial.id,
                  data=data1, metric="logLikRaw"))


  (p.mr.agg <- satf(dv=c(n.responses.yes=~n.responses.yes, n.responses=~n.responses),
                  signal=~signal,
                  start     = c(lambda=3, beta=1, delta=.4),
                  contrasts = c(lambda=~1, beta=~1, delta=~1), 
                  constraints=list(), time=~time, 
                  bias=list(bias=~1, poly.degree=10), #trial.id=~trial.id,
                  data=data.mr.nyes, metric="logLikRaw"))


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
