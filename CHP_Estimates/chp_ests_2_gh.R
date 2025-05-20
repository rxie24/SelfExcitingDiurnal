library(Rcpp)
library(parallel)

#Directory for base function files
#wd_files <- "/Users/xieryan/Desktop/Dissertation_1/Base Functions/" #Local
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

source(paste0(wd_files,"base_chp.R"))
source(paste0(wd_files,"chp_complls.R"))
source(paste0(wd_files,"rope_functions.R"))

chp_llest <- function(k,wd,wddata,rp=10,nfc=4){
  
  #Load in the data
  times <- readRDS(paste0(wddata,"data",k,".rds"))
  
  #Fit Hawkes Model
  thetai_hp <- rep(0,3)
  mles_hp <- hp_mles(thetai_hp,times)
  ll_hp <- -1*hp_ll(mles_hp,times)
  
  #Fit Circadian Hawkes Model w/ two sets of initial values, and find the one with the better fit
  thetai_chp1 <- rep(0,5)
  mles_chp1 <- chp_mles(thetai_chp1,times)

  thetai_chp2 <- c(mles_hp,rep(0,2))
  mles_chp2 <- chp_mles(thetai_chp2,times)
  
  #Calculate total sum of squares (across cumulative and max)
  #1st initial values
  ml_fs01 <- mles_chp1[4:5]
  sscum1 <- cum_ss(times,ml_fs01)
  ropen1 <- rpen(ml_fs01)
  ropen1 <- rp*ropen1
  ro_pe0cum1 <- sscum1 + ropen1
  
  #2nd initial values
  ml_fs02 <- mles_chp2[4:5]
  sscum2 <- cum_ss(times,ml_fs02)
  ropen2 <- rpen(ml_fs02)
  ropen2 <- rp*ropen2
  ro_pe0cum2 <- sscum2 + ropen2
  
  if(ro_pe0cum2 < ro_pe0cum1){
    ro_pe0cum <- ro_pe0cum2
    mles_chp <- mles_chp2
  }else{
    ro_pe0cum <- ro_pe0cum1
    mles_chp <- mles_chp1
  }
  
  for(j in 2:nfc){
    #First initial values
    thetai_chp1 <- rep(0,(3+(2*j)))
    par_estj1 <- chp_mles(thetai_chp1,times)
    
    ml_fsj1 <- par_estj1[4:length(par_estj1)]
    sscumj1 <- cum_ss(times,ml_fsj1)
    ropenj1 <- rpen(ml_fsj1)
    ropenj1 <- rp*ropenj1
    ro_pejcum1 <- sscumj1 + ropenj1
    
    #Second initial values
    thetai_chp2 <- c(mles_hp,rep(0,2*j))
    par_estj2 <- chp_mles(thetai_chp2,times)
    
    ml_fsj2 <- par_estj2[4:length(par_estj2)]
    sscumj2 <- cum_ss(times,ml_fsj2)
    ropenj2 <- rpen(ml_fsj2)
    ropenj2 <- rp*ropenj2
    ro_pejcum2 <- sscumj2 + ropenj2
    
    if(ro_pejcum2 < ro_pejcum1){
      ro_pejcum <- ro_pejcum2
      par_estj <- par_estj2
    }else{
      ro_pejcum <- ro_pejcum1
      par_estj <- par_estj1
    }
    
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
  
  #Fit Poisson Process
  ll_pp <- max_logpp(times)
  
  #Save 
  ll_est <- c(k,ll_pp,ll_cp,ll_hp,ll_chp) #Poisson, Circadian, Hawkes, CHP
  
  return(ll_est)
}

#Establish working directories to save parameter estimates and where data is located 
wd <- "/home/xieryan/Dissertation1/Data/CHP_GH/"
wddata <- "/home/xieryan/Dissertation1/Data/Real_Data_Event_Times_V2_GitHub/"

ids <- c(7,34,35,38,41)
all_lls <- mclapply(ids,chp_llest,wd=wd,wddata=wddata,rp=0,nfc=4,mc.cores=5) #have no roughness penalty

#Create matrix of all likelihoods for each individual
id_lls <- matrix(NA,nrow=5,ncol=5)
for(j in 1:5){
  id_lls[j,] <- all_lls[[j]]
}

#Save matrix of log-likelihoods when fitting four models
saveRDS(object = id_lls,file=paste0(wd,"CHP_ID_LLs_Final2.rds")) 

