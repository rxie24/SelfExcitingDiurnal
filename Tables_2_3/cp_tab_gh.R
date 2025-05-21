library(Rcpp)
library(parallel)

#Directory for base function files (MODIFY TO BE YOUR SPECIFIC DIRECTORY)
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

#Load in needed base function files
source(paste0(wd_files,"base_chp.R"))
source(paste0(wd_files,"base_lse.R"))
source(paste0(wd_files,"gen_data_chp.R"))
source(paste0(wd_files,"chp_complls.R"))
source(paste0(wd_files,"models_five.R"))
source(paste0(wd_files,"chp_five.R"))
source(paste0(wd_files,"filt_outliers.R"))

#Get set of scaled likelihoods under all five models for different sets of training and test datasets 
#generated from Circadian Process

RNGkind("L'Ecuyer-CMRG")
set.seed(215)

#Set common parameters to test for 
mu <- 0.1
beta <- 0
sigma <- 0
etagam_seq <- list(c(-0.3,0.3),c(-0.3,0.3,0.2,-0.2)) #Consider 2 and 4 Fourier coefficients
gl <- length(etagam_seq)

niter <- 300 #number of iterations to take average over
ncores <- 60

#Set row and column names for tables
cnames <- c("CHP","HP","PP","CP","LSE")
rnames <- c("(-0.3,0.3)","(-0.3,0.3,0.2,-0.2)")

#Table 1: N = 50
N <- 50
chp_testlls <- rep(NA,gl)
hp_testlls <- rep(NA,gl)
pp_testlls <- rep(NA,gl)
cp_testlls <- rep(NA,gl)
lse_testlls <- rep(NA,gl)

for(k in 1:gl){
  etagammas <- etagam_seq[[k]]
  theta_chp <- c(mu,beta,sigma,etagammas)
  thetai_chp <- rep(0,length(theta_chp))
  thetai_lse <- thetai_chp
  
  ll_meds <- five_chp_iter(niter,theta_chp,thetai_chp,thetai_lse,N,ncores=ncores)

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
chp_testlls <- rep(NA,gl)
hp_testlls <- rep(NA,gl)
pp_testlls <- rep(NA,gl)
cp_testlls <- rep(NA,gl)
lse_testlls <- rep(NA,gl)

for(k in 1:gl){
  etagammas <- etagam_seq[[k]]
  theta_chp <- c(mu,beta,sigma,etagammas)
  thetai_chp <- rep(0,length(theta_chp))
  thetai_lse <- thetai_chp
  
  ll_meds <- five_chp_iter(niter,theta_chp,thetai_chp,thetai_lse,N,ncores=ncores)
  
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


