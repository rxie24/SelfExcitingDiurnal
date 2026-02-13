library(Rcpp)
library(parallel)

#Directory for base function files
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

#Load in needed base function files
source(paste0(wd_files,"base_chp5.R"))
source(paste0(wd_files,"gen_data_chp.R"))
source(paste0(wd_files,"inference_func.R"))

wddata <- "/home/xieryan/Dissertation1/Data/Real_Data_Event_Times_V2_GitHub/"
wdest <- "/home/xieryan/Dissertation1/Data/CHP5_Final/"

#Set parameters needed for all bootstrapping confidence intervals
niter <- 100
ncores <- 60 #60

RNGkind("L'Ecuyer-CMRG")
set.seed(215)

#Look at parameter estimates from first 10 individuals
for(j in 11:20){
  
  #Download estimates for individual j
  times <- readRDS(paste0(wddata,"data",j,".rds"))
  n <- length(times)-1
  
  #Get New CHP Estimates for a given ID
  mles_chp <- readRDS(paste0(wdest,"ests",j,".rds"))
  print(mles_chp)
  
  #Now perform bootstrapping to try and get confidence intervals
  tryCatch({
    mle_bs_ci <- chp_par_bootstrap(mles_chp,n,ncores=ncores,niter=niter)
  },
  error=function(e){
    mle_bs_ci <<- rep(NaN,5)
  })
  print(mle_bs_ci)
}
