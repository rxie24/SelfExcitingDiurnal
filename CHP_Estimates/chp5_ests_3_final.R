library(Rcpp)
library(parallel)

#Directory for base function files
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

source(paste0(wd_files,"base_chp5.R"))
source(paste0(wd_files,"chp_complls.R"))

chp_llest <- function(k,wd,wddata,nfc){
  
  #Load in the data
  times <- readRDS(paste0(wddata,"data",k,".rds"))
  
  #Fit Hawkes Model
  thetai_hp <- rep(0,3)
  mles_hp <- hp_mles(thetai_hp,times)
  ll_hp <- -1*hp_ll(mles_hp,times)
  
  #Fit Circadian Process
  nfc_k <- nfc[k]
  thetai_cp <- rep(0,(1+nfc_k))
  mles_cp <- cp_mles(thetai_cp,times)
  ll_cp <- -1*cp_ll(mles_cp,times)
  
  #Fit Circadian Hawkes Model w/ Hawkes and Circadian parameters as initial values
  cp_init <- mles_cp[2:length(mles_cp)]
  thetai_chp <- c(mles_hp,cp_init)
  mles_chp <- chp_mles(thetai_chp,times)
  
  #Save mles of full Circadian Hawkes model
  saveRDS(object = mles_chp,file=paste0(wd,"ests",k,".rds"))  
  ll_chp <- -1*chp_ll(mles_chp,times)
  
  #Fit Poisson Process
  ll_pp <- max_logpp(times)
  
  #Save 
  ll_est <- c(k,ll_pp,ll_cp,ll_hp,ll_chp) #Poisson, Circadian, Hawkes, CHP
  
  return(ll_est)
}

#Establish working directories to save parameter estimates and where data is located 
wd <- "/home/xieryan/Dissertation1/Data/CHP5_Final/"
wddata <- "/home/xieryan/Dissertation1/Data/Real_Data_Event_Times_V2_GitHub/"

#Load in number of Fourier coefficients to be used for all individuals
nfc <- readRDS("/home/xieryan/Dissertation1/Data/CHP4/chp_lengths.rds")

#Fit CHP Model for individual
id <- 34
id_lls <- chp_llest(id,wd=wd,wddata=wddata,nfc=nfc) #have no roughness penalty
print(id_lls)

#Save matrix of log-likelihoods when fitting four models
saveRDS(object = id_lls,file=paste0(wd,"CHP_ID_LLs3.rds")) 

