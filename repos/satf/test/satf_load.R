
library(Rcpp)
source("../R/satf_misc.R")
source("../R/satf_aggregate.R")
source("../R/satf_params.R")
source("../R/satf_generate.R")
source("../R/satf.R")
Rcpp::sourceCpp("../src/rcpp_exports.cpp", rebuild=T)
