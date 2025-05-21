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
source(paste0(wd_files,"chp_five_sim.R"))
source(paste0(wd_files,"filt_outliers.R"))

#Get set of scaled likelihoods under all five models for different sets of training and test datasets 
#generated from Hawkes Process

RNGkind("L'Ecuyer-CMRG")
set.seed(215)

#Set common parameters to test for 
mu <- 0.1
beta <- 0.25
eta <- 0
sigma_seq <- c(0.2,-2)
gammas <- 0

theta_wb <- c(mu,beta,eta,gammas)

#Set initial parameter estimates
thetai_chp <- rep(0,5)
thetai_lse <- rep(0,5)

niter <- 300 #number of iterations to take average over

#Set row and column names for tables
rnames <- c("Log(Sigma)=0.2","Log(Sigma)=-2")
cnames <- c("CHP","HP","PP","CP","LSE")

#Table 1: N = 50
N <- 50
five_testlls <- test_lls(sigma_seq,niter,theta_wb,3,thetai_chp,thetai_lse,N,ncores=60)
chp_testlls <- five_testlls[[1]]
hp_testlls <- five_testlls[[2]]
pp_testlls <- five_testlls[[3]]
cp_testlls <- five_testlls[[4]]
lse_testlls <- five_testlls[[5]]

tab1_results <- cbind(chp_testlls,hp_testlls,pp_testlls,cp_testlls,lse_testlls)
colnames(tab1_results) <- cnames
rownames(tab1_results) <- rnames

#Table 2: N = 500
N <- 500
five_testlls <- test_lls(sigma_seq,niter,theta_wb,3,thetai_chp,thetai_lse,N,ncores=60)
chp_testlls <- five_testlls[[1]]
hp_testlls <- five_testlls[[2]]
pp_testlls <- five_testlls[[3]]
cp_testlls <- five_testlls[[4]]
lse_testlls <- five_testlls[[5]]

tab2_results <- cbind(chp_testlls,hp_testlls,pp_testlls,cp_testlls,lse_testlls)
colnames(tab2_results) <- cnames
rownames(tab2_results) <- rnames

print("N=50")
print(tab1_results)

print("N=500")
print(tab2_results)


