  library(ggplot2)

Eval.SATF.Aggregate <- function(params, fixedparams, metric, action="fit")
{

  data <- fixedparams$data
  n.data <- nrow(data)
  n.params <- length(params); 

  # process parameters
  params <- process.parameters(params, fixedparams$fixed, 
                               fixedparams$params.idx,
                               fixedparams$transform.fn,
                               fixedparams$fn.params.LL)
  # compute d' predictions from the SAT function 
  predicted.dprime <- predict.dprime(params$satf.params, 
                                     fixedparams$design.matrix, data)
  
  # determine criterion function
  if('predicted.criterion' %in% names(fixedparams)) {
    predicted.criterion <- fixedparams$predicted.criterion
  } else {
    fn.crit <- function(t) fixedparams$fn.crit(t, params$crit.params)
    predicted.criterion <- fn.crit(data$time)
  }
  
  if(action=="predict.dprime") {
    return(data.frame(predicted.dprime, predicted.criterion))
  }
  
  stopifnot( !is.na(data$dprime) )
  mean.dprime <- mean( data$dprime )
  
  # compute metric
   if( metric == "RMSD" ) {
      rmsd <- sum( (data$dprime-predicted.dprime)^2 )
      res <- rmsd
    }
    else if( metric == "R2" ) {
      SSresidual   <- sum( (data$dprime-predicted.dprime)^2 )
      SStotal <- sum( (data$dprime-mean.dprime)^2 )
      R2 <- (1-SSresidual/SStotal)
      res <- R2
    }
    else if( metric == "adjR2" ) {
      SSresidual   <- sum( (data$dprime-predicted.dprime)^2 )
      SStotal <- sum( (data$dprime-mean.dprime)^2 )
      adjR2numerator <- SSresidual/ (n.data-n.params)
      adjR2denominator <- SStotal / (n.data-1)
      adjR2 <- (1-adjR2numerator/adjR2denominator)
      res <- adjR2
    }
    else if( metric == "logLik" ) {
      # Binomial likelihood, found by adjusting dprime and bias (c)
      # to maximize likelihood
      fa.defined <- !is.na(data$fas) & !is.na(data$crs)
        
      predicted.hit.prob <- 1-pnorm(predicted.criterion, mean=predicted.dprime)
      predicted.fa.prob  <- 1-pnorm(predicted.criterion, mean=0)
      n.grammatical <- data$hits+data$misses
      logLikHits <- dbinom(data$hits, prob=predicted.hit.prob, size=n.grammatical, log=T)

      n.ungrammatical <- (data$fas+data$crs)[fa.defined]
      fas <- data$fas[fa.defined]
      predicted.fa.prob <- predicted.fa.prob[fa.defined]
      logLikFAs <- dbinom(fas, prob=predicted.fa.prob, size=n.ungrammatical, log=T)
      
      if( 'trial.noise.s' %in% names(params$repeat.params) ) { # Assume a certain amount of dependence between values
        trial.noise.s <- params$repeat.params[['trial.noise.s']]
        dependent <- !is.na(data$last.hits)
        
        # determine differences in d' and criterion
        diff.crit.psi.hits <- predicted.dprime-predicted.criterion
        last.diff.crit.psi.hits <- c(0, diff.crit.psi.hits[-nrow(data)])
        diff.crit.psi.fas <- 0-predicted.criterion
        last.diff.crit.psi.fas <- c(0, diff.crit.psi.fas[-nrow(data)])
        
        probConditionalHits <- rep(NA, nrow(data))
        probConditionalFAs <- rep(NA, nrow(data))
        for(i in (1:nrow(data))[!is.na(data$last.hits)] ) {
          stopifnot(data$hits[i]+data$misses[i] == data$last.hits[i]+data$last.misses[i])
          probYesAfterYes <- prob.dependent.yes(diff.crit.psi.hits[i], trial.noise.s, last.diff.crit.psi.hits[i], last.response.yes=1)
          probYesAfterNo <- prob.dependent.yes(diff.crit.psi.hits[i], trial.noise.s, last.diff.crit.psi.hits[i], last.response.yes=0)
          
          # obtain probabilities for all possible scenarios which would lead to the obtained number of hits
          probYesAfterYes <- dbinom(0:data$hits[i], size=data$last.hits[i], prob=probYesAfterYes, log=T)
          probYesAfterNo <- dbinom(data$hits[i]:0, size=data$last.misses[i], prob=probYesAfterNo, log=T)
          probConditionalHits[i] <- sum(exp(probYesAfterYes+probYesAfterNo))
          
          if(!is.na(data$fas[i]) && !is.na(data$crs[i])) {
            stopifnot(data$fas[i]+data$crs[i] == data$last.fas[i]+data$last.crs[i])
            probYesAfterYes <- prob.dependent.yes(diff.crit.psi.fas[i], trial.noise.s, 
                                              last.diff.crit.psi.fas[i], last.response.yes=1)
            probYesAfterNo <- prob.dependent.yes(diff.crit.psi.fas[i], trial.noise.s, 
                                              last.diff.crit.psi.fas[i], last.response.yes=0)
            probYesAfterYes <- dbinom(0:data$fas[i], size=data$last.fas[i], prob=probYesAfterYes, log=T)
            probYesAfterNo <- dbinom(data$fas[i]:0, size=data$last.crs[i], prob=probYesAfterNo, log=T)
            probConditionalFAs[i] <- sum(exp(probYesAfterYes+probYesAfterNo))
          }
        }
        stopifnot(all(probConditionalHits[dependent] >= 0 & probConditionalHits[dependent] <= 1))
        stopifnot(all(probConditionalFAs[dependent & fa.defined] >= 0 & probConditionalFAs[dependent & fa.defined] <= 1))
        
        logLikHits <- sum(logLikHits[!dependent]) + sum(log(probConditionalHits[dependent]))
        logLikFAs <- sum(logLikFAs[!dependent & fa.defined]) + sum(log(probConditionalFAs[dependent & fa.defined]))
        
        res <- logLikHits + logLikFAs
      } 
      else { 
      # Assume that all values are independent
        res <- sum(logLikHits) + sum(logLikFAs)
      }
      res.bak <- res
      res <- res + params$LL.params
      #print(c(res.bak, res))
    }
  return(res)
}

plot.fit <- function(params, fixedparams, control, predicted=NULL) 
{
  data <- fixedparams$data
  check.columns(c('dprime',control$condition), data)
  corr <- control$correction.fraction
  
#  data[,'condition'] <- data[,control$condition]
  
  # obtain predicted dprime values
  design.matrix <- params$satf.params*fixedparams$design.matrix$design.matrix
  lambda <- col.sums(design.matrix[fixedparams$design.matrix$lambda,])
  beta <- col.sums(design.matrix[fixedparams$design.matrix$beta,])
  delta <- col.sums(design.matrix[fixedparams$design.matrix$delta,])
  plotdata <- data.frame(lambda, beta, delta, condition=data[,c('condition')])
  plotdata <- unique(plotdata)

  newtime <- rep(seq(min(data$time), max(data$time), by=.001), each=nrow(plotdata))
  plotdata <- cbind(time=newtime, plotdata[rep(1:nrow(plotdata), length(newtime)),] )
  plotdata$predicted.dprime <- with(plotdata, SATF(time, lambda, beta, delta))
  plotdata$predicted.criterion <- fixedparams$fn.criterion(plotdata$time)
  plotdata <- plotdata[,c('time','condition','predicted.dprime','predicted.criterion')]

  cr.prop <- with(data, (crs+corr)/(fas+crs+2*corr) )
  data$emp.criterion[!is.na(cr.prop)] <- qnorm(cr.prop[!is.na(cr.prop)])
  if(control$plot.dprime || control$plot.criterion) {
    if(control$plot.dprime) {
      print(head(plotdata))
      p <- ggplot(plotdata, aes(x=time,y=predicted.dprime,linetype=condition))+geom_line()
      p <- p+geom_point(data=data, aes(x=time, y=dprime,linetype=condition))
    }
    if(control$plot.criterion) {
      data.cr <- subset(data, !is.na(emp.criterion))
      p <- p+geom_point(data=data.cr, aes(x=time, y=emp.criterion, color="criterion"))
      p <- p+geom_line(data=data.cr, aes(x=time, y=predicted.criterion, color="criterion")) 
    }
    if(control$plot.action == "return")
      return(p)
    else 
      print(p)
  }
}


# TODO: Change the code in satf_rawbinomial.R to fit the bias as well,
#       like here, instead of the criterion.
# TODO: Implement the plotting function present in  in satf_rawbinomial.R 
# TODO: Add code which checks whether ungrammatical data (fas and crs)
#       has been repeated for every condition.
Fit.SATF.Aggregate <- function(start, lambda, beta, delta, data, control,
                               transform.fn, transform.vec, params.minimum, params.maximum,
                               fixed=c(), metric, fn.params.LL=NULL,
                               optim.control=list(), optim.digits=1, SANN.maxit=0, ...)
{
  control$criterion.intervals <- 0

  # check and inizialize parameters
  logLikCols <- c('hits','fas','crs','misses')
  logLikMRCols <- c(logLikCols, paste('last',logLikCols,sep="."))
  if( metric == "RMSD") check.columns( 'dprime', data)
  else if( metric == "R2" ) check.columns( 'dprime', data)
  else if( metric == "adjR2" ) check.columns( 'dprime', data)
  else if( metric == "logLik" ) {
    if('trial.noise.s' %in% names(start))
      check.columns( logLikMRCols, data)
    else
      check.columns( logLikCols, data)
  }
  
  optim.control <- default(optim.control, 'maxit', 10^6)
  if(metric=="RMSD") optim.control$fnscale = 1
  else optim.control$fnscale = -1
  
  # create the design matrix
  dm <- design.matrix(lambda=lambda, beta=beta, delta=delta, data=data,
                      stimulus.signal=rep(1, nrow(data)))
  
  # init the fixed parameters to Eval.SATF.Aggregate()
  fixedparams <- list(design.matrix=dm, fixed=fixed, data=data,
                      transform.fn=transform.fn, fn.params.LL=fn.params.LL)

  params.idx <- index.params(c(start, fixed))
  fixedparams$params.idx <- params.idx

  fn.criterion <- ComputeEmpiricalBiasFn(data=data, control=control, 
                                         signal=NULL, response=NULL)
  fixedparams$fn.criterion <- fn.criterion
  fixedparams$predicted.criterion <- fn.criterion(data$time)

  # maximize the log-likelihood
  if(metric %in% c("logLik", "adjR2"))
  {
    if(metric == "logLik") {
      response.hits=data$hits; n.grammatical=data$hits+data$misses;
      response.fas=data$fas; n.ungrammatical=data$fas+data$crs;
      response.dprime=NULL;
    } else {
      response.hits=NULL; n.grammatical=NULL;
      response.fas=NULL; n.ungrammatical=NULL;
      response.dprime=data$dprime;
    }
    # set the relevant parameters necessary for optimization
    fp <- initialize.compute.logLik(design.matrix=dm, transform.vec=transform.vec,
                         params.minimum=params.minimum, params.maximum=params.maximum, time=data$time,
                         predicted.criterion=fn.criterion(data$time),
                         response.dprime=response.dprime,
                         response.hits=response.hits, n.grammatical=n.grammatical,
                         response.fas=response.fas, n.ungrammatical=n.ungrammatical,
                         fn.params.LL=fn.params.LL)
    res <- compute.logLik(start,  fixedparams=fp)
  
    if( abs(res) == Inf || is.na(res) ) {
      print(start)
      LL <- compute.logLik(start, fixedparams=fp, return.LL.params=T)
      deinitialize.compute.logLik(fp)
      msg <- sprintf("Function evaluates to %s at starting values. Returning NULL.", res)
      warning(msg)
      return(NULL)
    }
    if(SANN.maxit > 0) {
      res <- optim(start, fn=compute.logLik, fixedparams=fp, method="SANN", control=list(maxit=SANN.maxit))
      start <- res$par
      print(start)
    }
    res <- optim.to.precision(start, optim.digits=optim.digits, fn=compute.logLik, 
                              fixedparams=fp, control=optim.control)
    LL <- compute.logLik(res$par, fixedparams=fp,  coefs.len=length(start),
                         return.LL.params=T)
    oldpar <- res$par
    deinitialize.compute.logLik(fp)
  } else {
   res <- optim.to.precision(start, optim.digits, fn=Eval.SATF.Aggregate, 
                    fixedparams=fixedparams, metric=metric, control=optim.control)
  }
  res$par <- process.parameters(res$par, fixedparams$fixed, fixedparams$params.idx,
                                fixedparams$transform.fn, params.minimum=params.minimum, 
                                params.maximum=params.maximum, fixedparams$fn.params.LL)
  if(control$plot.dprime) {
    plot <- plot.fit(res$par, fixedparams, control)
  }
  res$plot <- plot
  res  
}
