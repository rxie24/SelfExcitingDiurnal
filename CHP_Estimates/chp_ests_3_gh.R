library(Rcpp)
library(parallel)

#Directory for base function files (MODIFY TO BE YOUR SPECIFIC DIRECTORY)
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

source(paste0(wd_files,"base_chp.R"))
source(paste0(wd_files,"chp_complls.R"))

chp_llest <- function(k,wd,wddata){
  
  #Load in the data
  times <- readRDS(paste0(wddata,"data",k,".rds"))
  
  #Fit Hawkes Model
  thetai_hp <- rep(0,3)
  mles_hp <- hp_mles(thetai_hp,times)
  ll_hp <- -1*hp_ll(mles_hp,times)
  
  #Fit Circadian Hawkes Model w/ fixed set of coefficients (6)
  thetai_chp <- rep(0,9)
  mles_chp <- chp_mles(thetai_chp,times)
  saveRDS(object = mles_chp,file=paste0(wd,"ests",k,".rds")) #Save mles of full Circadian Hawkes model
  ll_chp <- -1*chp_ll(mles_chp,times)
  
  #Fit Circadian Process
  lfc <- length(mles_chp[4:length(mles_chp)])
  thetai_cp <- rep(0,(1+lfc))
  mles_cp <- cp_mles(thetai_cp,times)
  ll_cp <- -1*cp_ll(mles_cp,times)
  
  #Fit Poisson Process
  ll_pp <- max_logpp(times)
  
  #Save 
  ll_est <- c(k,ll_pp,ll_cp,ll_hp,ll_chp) #Poisson, Circadian, Hawkes, CHP
  
  return(ll_est)
}

ids <- c(10,13,17)

#Establish working directories to saved CHP parameter estimates and where data is located (MODIFY TO BE YOUR SPECIFIC DIRECTORY)
wd <- "/home/xieryan/Dissertation1/Data/CHP_GH/"
wddata <- "/home/xieryan/Dissertation1/Data/Real_Data_Event_Times_V2_GitHub/"

all_lls <- mclapply(ids,chp_llest,wd=wd,wddata=wddata,mc.cores=3) #have no roughness penalty
print(all_lls)

#Create matrix of all likelihoods for each individual
id_lls <- matrix(NA,nrow=3,ncol=5)
for(j in 1:3){
  id_lls[j,] <- all_lls[[j]]
}

#Save matrix of log-likelihoods when fitting four models
saveRDS(object = id_lls,file=paste0(wd,"CHP_ID_LLs_Final3.rds")) 
