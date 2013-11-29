# TODO: remove the colums 'dprime', 'criterion', 'noise', and 'dprime.cur' from generated data frames,
#       or at least give them better names.

process.arg <- function(arg, time) {
  if(typeof(arg) == "closure") arg = arg(time)
  if(length(arg) == 1 ) arg = rep(arg, length(time))
  stopifnot(length(arg) == length(time) )
  arg
}

satf_generate_condition <- function(criterion, dprime, time, trial.ids, label, rho=0)
{
  intervals <- 1:length(time)
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

satf_generate <- function(criterion, dprime, time, n, rho=0, label="condition1") {
  data0 <- satf_generate_condition(criterion=criterion, dprime=0, time=time, 
                                  trial.ids=1:n, label=label, rho=rho)
  data0$signal <- 0
  data1 <- satf_generate_condition(criterion=criterion, dprime=dprime,  time=time, 
                                  trial.ids=1:n+n, label=label, rho=rho)
  data1$signal <- 1
  rbind(data0, data1)
}
