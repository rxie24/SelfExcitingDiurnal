library(Rcpp)
library(parallel)

#Directory for base function files
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

#Load in needed base function files
source(paste0(wd_files,"base_chp5.R"))
source(paste0(wd_files,"gen_data_chp.R"))
source(paste0(wd_files,"inference_func.R"))

#Function to generate data from waiting times for Circadian/Hawkes Models
gen_data_chp <- function(n,mu,beta,sigma,eta,gammas,subdiv=100000){
  #n: number of new points you want to generate
  times <- c(0,rep(NA,n)) #WLOG, set times[1] = 0
  
  #Generate waiting times and points
  for(k in 1:n){
    ur <- runif(1,min=0,max=1)
    intstart <- times[k]
    tsubset <- times[1:k]
    et_optim <- uniroot(et_sim_chp,c(0,200),intstart=intstart,times=tsubset,
                        mu=mu,beta=beta,sigma=sigma,eta=eta,gammas=gammas,ur=ur,subdiv=subdiv)
    et_sim <- et_optim$root
    times[k+1] <- times[k] + et_sim
  }
  return(times)
}

wddata <- "/home/xieryan/Dissertation1/Data/Real_Data_Event_Times_V2_GitHub/"
wdest <- "/home/xieryan/Dissertation1/Data/CHP5_Final/"

#Set parameters needed for all bootstrapping confidence intervals
niter <- 100
ncores <- 60 #60
ids <- c(23,41)

RNGkind("L'Ecuyer-CMRG")
set.seed(215)

#Look at parameter estimates from first 10 individuals
for(j in 1:2){
  
  id <- ids[j]
  
  #Download estimates for individual j
  times <- readRDS(paste0(wddata,"data",id,".rds"))
  n <- length(times)-1
  
  #Get New CHP Estimates for a given ID
  mles_chp <- readRDS(paste0(wdest,"ests",id,".rds"))
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

#Problem with 41