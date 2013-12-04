
library(plyr)
library(Rcpp)
source("~/CodeSATF/R/satf_misc.R")
source("~/CodeSATF/R/satf_aggregate.R")
source("~/CodeSATF/R/satf_params.R")
source("~/CodeSATF/R/satf_generate.R")
source("~/CodeSATF/R/satf.R")
Rcpp::sourceCpp("~/CodeSATF/src/rcpp_exports.cpp", rebuild=T)
