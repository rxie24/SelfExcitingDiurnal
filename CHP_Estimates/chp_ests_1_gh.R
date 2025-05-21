library(Rcpp)
library(parallel)

#Directory for base function files (MODIFY TO FIT YOUR SPECIFIC DIRECTORY)
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

source(paste0(wd_files,"base_chp.R"))
source(paste0(wd_files,"chp_complls.R"))
source(paste0(wd_files,"rope_functions.R"))

chp_llest <- function(k,wd,wddata,rp=10,nfc=4){
  
  #Load in the data
  times <- readRDS(paste0(wddata,"data",k,".rds"))
  
  #Fit Circadian Hawkes Model
  thetai_chp <- rep(0,5)
  mles_chp <- chp_mles(thetai_chp,times)

  #Calculate total sum of squares (across cumulative and max)
  ml_fs0 <- mles_chp[4:5]
  sscum <- cum_ss(times,ml_fs0)
  ropen <- rpen(ml_fs0)
  ropen <- rp*ropen
  
  ro_pe0cum <- sscum + ropen
  
  for(j in 2:nfc){
    thetai_chp <- rep(0,(3+(2*j)))
    par_estj <- chp_mles(thetai_chp,times)
    
    ml_fsj <- par_estj[4:length(par_estj)]
    sscumj <- cum_ss(times,ml_fsj)
    ropenj <- rpen(ml_fsj)
    ropenj <- rp*ropenj
    
    ro_pejcum <- sscumj + ropenj
    
    #If current penalized SS is greater than before, stop loop; otherwise, update new estimate
    if(ro_pejcum >= ro_pe0cum){
      break
    }else{
      ro_pe0cum <- ro_pejcum
      mles_chp <- par_estj
    }
  }
  
  #Save mles of full Circadian Hawkes model
  saveRDS(object = mles_chp,file=paste0(wd,"ests",k,".rds")) 
  ll_chp <- -1*chp_ll(mles_chp,times)
  
  #Fit Circadian Process
  lfc <- length(mles_chp[4:length(mles_chp)])
  thetai_cp <- rep(0,(1+lfc))
  mles_cp <- cp_mles(thetai_cp,times)
  ll_cp <- -1*cp_ll(mles_cp,times)
  
  #Fit Hawkes Model
  thetai_hp <- rep(0,3)
  mles_hp <- hp_mles(thetai_hp,times)
  ll_hp <- -1*hp_ll(mles_hp,times)
  
  #Fit Poisson Process
  ll_pp <- max_logpp(times)
  
  #Save 
  ll_est <- c(k,ll_pp,ll_cp,ll_hp,ll_chp) #Poisson, Circadian, Hawkes, CHP
  
  return(ll_est)
}

#Establish working directories to saved CHP parameter estimates and where data is located (MODIFY TO FIT YOUR SPECIFIC DIRECTORY)
wd <- "/home/xieryan/Dissertation1/Data/CHP_GH/"
wddata <- "/home/xieryan/Dissertation1/Data/Real_Data_Event_Times_V2_GitHub/"
ncores <- 33

ids <- c(1,2,3,4,5,6,8,9,11,12,14,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,36,37,39,40)
all_lls <- mclapply(ids,chp_llest,wd=wd,wddata=wddata,rp=0.005,nfc=4,mc.cores=ncores)

#Create matrix of all likelihoods for each individual (41 rows for all individuals)
id_lls <- matrix(NA,nrow=33,ncol=5)
for(j in 1:33){
  id_lls[j,] <- all_lls[[j]]
}

#Save matrix of log-likelihoods when fitting four models
saveRDS(object = id_lls,file=paste0(wd,"CHP_ID_LLs_Final.rds")) 
