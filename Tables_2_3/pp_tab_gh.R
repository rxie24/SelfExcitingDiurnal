library(Rcpp)
library(parallel)

#Directory for base function files
#wd_files <- "/Users/xieryan/Desktop/Dissertation_1/Base Functions/" #Local
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

#Load in needed base function files
source(paste0(wd_files,"base_chp.R"))
source(paste0(wd_files,"base_lse.R"))
source(paste0(wd_files,"pp_gen_data.R"))
source(paste0(wd_files,"chp_complls.R"))
source(paste0(wd_files,"models_five.R"))
source(paste0(wd_files,"filt_outliers.R"))

#Fit five models for data generated from Poisson Process
five_pp <- function(i,mu,thetai_chp,thetai_lse,N){
  
  #i: iteration number (for running tasks in parallel)
  #mu: parameter used to generate test and training data sets
  #thetai_chp: initial parameters used for Circadian Hawkes Process
  #thetai_lse: initial parameters used for Linear Self Exciting Processes
  #N: number of waiting times in training and test sets
  
  #Simulate training data and test data
  times_train <- gen_data_pp(N,mu)
  times_test <- gen_data_pp(N,mu)
  
  mlls <- five_models(times_train,times_test,thetai_chp,thetai_lse)
  return(mlls)
}

#Perform above function over multiple iterations, filter, and find median
five_pp_iter <- function(navg,mu,thetai_chp,thetai_lse,N,ncores=60){
  mod_lls <- mclapply(1:navg,five_pp,mu=mu,thetai_chp=thetai_chp,
                      thetai_lse=thetai_lse,N=N, mc.cores = ncores)
  
  chp_lls <- rep(NA,navg)
  hp_lls <- rep(NA,navg)
  pp_lls <- rep(NA,navg)
  cp_lls <- rep(NA,navg)
  lse_lls <- rep(NA,navg)
  
  for(j in 1:navg){
    all_lls <- mod_lls[[j]]
    chp_lls[j] <- all_lls[1]
    hp_lls[j] <- all_lls[2]
    pp_lls[j] <- all_lls[3]
    cp_lls[j] <- all_lls[4]
    lse_lls[j] <- all_lls[5]
  }
  
  #Filter outliers
  chp_lls <- filt(chp_lls)
  hp_lls <- filt(hp_lls)
  pp_lls <- filt(pp_lls)
  cp_lls <- filt(cp_lls)
  lse_lls <- filt(lse_lls)
  
  #Track number of NA log-likelihoods there are and print them
  no_NAs <- c(sum(is.na(chp_lls)),sum(is.na(hp_lls)),sum(is.na(pp_lls)),sum(is.na(cp_lls)),
              sum(is.na(lse_lls)))
  print(no_NAs)
  
  #Calculate median of loglikelihoods
  chp_avg <- median(chp_lls,na.rm=TRUE)
  hp_avg <- median(hp_lls,na.rm=TRUE)
  pp_avg <- median(pp_lls,na.rm=TRUE)
  cp_avg <- median(cp_lls,na.rm=TRUE)
  lse_avg <- median(lse_lls,na.rm=TRUE)
  
  all_meds <- c(chp_avg,hp_avg,pp_avg,cp_avg,lse_avg)
  return(all_meds)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(215)

#Get set of scaled likelihoods under all five models for different sets of training and test datasets 
#generated from Circadian Hawkes Process

#Set common parameters to test for 
mu <- 0.1

#Set initial parameter estimates (80% of original estimates)
thetai_chp <- rep(0,5)
thetai_lse <- thetai_chp

niter <- 300 #number of iterations to take average over
ncores <- 60

#Set row and column names for tables
cnames <- c("CHP","HP","PP","CP","LSE")

#Table 1: N = 50
N <- 50

ll_meds <- five_pp_iter(niter,mu,thetai_chp,thetai_lse,N,ncores=ncores)
tab1_results <- cbind(ll_meds)
row.names(tab1_results) <- cnames

#Table 2: N = 500
N <- 500

ll_meds <- five_pp_iter(niter,mu,thetai_chp,thetai_lse,N,ncores=ncores)
tab2_results <- cbind(ll_meds)
row.names(tab2_results) <- cnames

print("N=50")
print(tab1_results)

print("N=500")
print(tab2_results)
