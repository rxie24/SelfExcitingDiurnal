library(Rcpp)
library(parallel)

#Directory for base function files (MODIFY TO BE YOUR OWN DIRECTORY)
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

#Load in needed base function files
source(paste0(wd_files,"base_chp.R"))
source(paste0(wd_files,"base_lse.R"))
source(paste0(wd_files,"gen_data_lse2.R"))
source(paste0(wd_files,"chp_complls.R"))
source(paste0(wd_files,"models_five.R"))
source(paste0(wd_files,"filt_outliers.R"))

#Fit five models for data generated from Circadian Linear Self-Exciting Model
five_lse <- function(i,theta_lse,thetai_chp,thetai_lse,N,b=2){
  
  #i: iteration number (for running tasks in parallel)
  #theta_lse: parameters used to generate test and training data sets
  #thetai_chp: initial parameters used for Circadian Hawkes Process
  #thetai_lse: initial parameters used for Circadian Linear Self Exciting Processes
  #N: number of waiting times in training and test sets
  #b: order of linear self-exciting model
  
  #Extract parameters to generate data
  e_ind <- 2 + b
  
  ltheta <- length(theta_lse)
  mu <- theta_lse[1]
  beta <- theta_lse[2:(e_ind-1)]
  
  ncoef <- length(theta_lse[e_ind:ltheta])/2
  eta <- theta_lse[e_ind:(e_ind+ncoef-1)]
  gammas <- theta_lse[(e_ind+ncoef):ltheta]
  
  #Simulate training data and test data
  times_train <- gen_data_lse(N,mu,beta,eta,gammas)
  times_test <- gen_data_lse(N,mu,beta,eta,gammas)
  
  mlls <- five_models(times_train,times_test,thetai_chp,thetai_lse)
  return(mlls)
}

#Perform above function over multiple iterations, filter, and find median
five_lse_iter <- function(navg,theta_lse,thetai_chp,thetai_lse,N,ncores=60){
  mod_lls <- mclapply(1:navg,five_lse,theta_lse=theta_lse,thetai_chp=thetai_chp,
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
#generated from Circadian Linear Self-Exciting Process

#Set common parameters to test for 
mu <- 0.1
beta_seq <- list(c(0.2,0.1),c(0.1,0.2))
eta <- -0.3
gammas <- 0.3

thetai_chp <- rep(0,5)
thetai_lse <- rep(0,5)

niter <- 300 #number of iterations to take average over
bl <- length(beta_seq)
ncores <- 60

#Set row and column names for tables
cnames <- c("CHP","HP","PP","CP","LSE")
rnames <- c("(0.2,0.1)","(0.1,0.2)")

#Table 1: N = 50
N <- 50
chp_testlls <- rep(NA,bl)
hp_testlls <- rep(NA,bl)
pp_testlls <- rep(NA,bl)
cp_testlls <- rep(NA,bl)
lse_testlls <- rep(NA,bl)

for(k in 1:bl){
  beta <- beta_seq[[k]]
  theta_lse <- c(mu,beta,eta,gammas)
  
  ll_meds <- five_lse_iter(niter,theta_lse,thetai_chp,thetai_lse,N,ncores=ncores)
  
  chp_testlls[k] <- ll_meds[1]
  hp_testlls[k] <- ll_meds[2]
  pp_testlls[k] <- ll_meds[3]
  cp_testlls[k] <- ll_meds[4]
  lse_testlls[k] <- ll_meds[5]
}

tab1_results <- cbind(chp_testlls,hp_testlls,pp_testlls,cp_testlls,lse_testlls)
colnames(tab1_results) <- cnames
rownames(tab1_results) <- rnames

#Table 2: N = 500
N <- 500
chp_testlls <- rep(NA,bl)
hp_testlls <- rep(NA,bl)
pp_testlls <- rep(NA,bl)
cp_testlls <- rep(NA,bl)
lse_testlls <- rep(NA,bl)

for(k in 1:bl){
  beta <- beta_seq[[k]]
  theta_lse <- c(mu,beta,eta,gammas)
  
  ll_meds <- five_lse_iter(niter,theta_lse,thetai_chp,thetai_lse,N,ncores=ncores)
  
  chp_testlls[k] <- ll_meds[1]
  hp_testlls[k] <- ll_meds[2]
  pp_testlls[k] <- ll_meds[3]
  cp_testlls[k] <- ll_meds[4]
  lse_testlls[k] <- ll_meds[5]
}

tab2_results <- cbind(chp_testlls,hp_testlls,pp_testlls,cp_testlls,lse_testlls)
colnames(tab2_results) <- cnames
rownames(tab2_results) <- rnames

print("N=50")
print(tab1_results)

print("N=500")
print(tab2_results)
