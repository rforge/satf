
library(sets)
library(plyr)
library(Rcpp)
library(maxLik)
library(nloptr)
source("~/CodeSATF/R/satf_misc.R")
source("~/CodeSATF/R/satf_aggregate.R")
source("~/CodeSATF/R/satf_params.R")
source("~/CodeSATF/R/satf_generate.R")
source("~/CodeSATF/R/satf_resample.R")
source("~/CodeSATF/R/satf.R")
Rcpp::sourceCpp("~/CodeSATF/test/satf_load.cpp", rebuild=T)
