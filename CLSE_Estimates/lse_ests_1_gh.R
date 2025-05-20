library(Rcpp)
library(parallel)

#Directory for base function files
#wd_files <- "/Users/xieryan/Desktop/Dissertation_1/Base Functions/" #Local
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

#Load in base function files
source(paste0(wd_files,"base_lse.R"))
source(paste0(wd_files,"lse_complls.R"))
source(paste0(wd_files,"rope_functions.R"))

lse_llest <- function(k,wd,wddata,rp=10,nfc=4,b=2){
  
  #Load in the data
  times <- readRDS(paste0(wddata,"data",k,".rds"))
  
  thetai_lse <- rep(0,5)
  mles_lse <- lse_mles(thetai_lse,times,b=b)
  
  #Calculate total sum of squares (across cumulative and max)
  ml_fs0 <- mles_lse[4:5]
  sscum <- cum_ss(times,ml_fs0)
  ropen <- rpen(ml_fs0)
  ropen <- rp*ropen
  
  ro_pe0cum <- sscum + ropen
  
  for(j in 2:nfc){
    thetai_lse <- rep(0,(3+(2*j)))
    par_estj <- lse_mles(thetai_lse,times,b=b)
    
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
      mles_lse <- par_estj
    }
  }
  
  #Save mles of full linear self-exciting model
  saveRDS(object = mles_lse,file=paste0(wd,"ests",k,".rds"))
  ll_lse <- -1*lse_ll(mles_lse,times,b=b)
  
  #Fit Circadian Process
  lfc <- length(mles_lse[4:length(mles_lse)])
  thetai_cp <- rep(0,(1+lfc))
  mles_cp <- cp_mles2(thetai_cp,times,b=b)
  ll_cp <- -1*cp_ll2(mles_cp,times,b=b)
  
  #Fit Linear Self-Exciting Model w/o Circadian Component
  thetai_nc <- rep(0,3)
  mles_nc <- lse_nc_mles(thetai_nc,times,b=b)
  ll_nc <- -1*lse_nc_ll(mles_nc,times,b=b)
  
  #Fit Poisson Process
  ll_pp <- max_logpp2(times,b=b)
  
  #Save 
  ll_est <- c(k,ll_pp,ll_cp,ll_nc,ll_lse) #Poisson, Circadian, LSE w/o FS, LSE
  
  return(ll_est)
}

#Establish working directories to save parameter estimates and where data is located 
wd <- "/home/xieryan/Dissertation1/Data/LSE_GH/"
wddata <- "/home/xieryan/Dissertation1/Data/Real_Data_Event_Times_V2_GitHub/"

ids <- c(2,3,4,5,8,9,11,12,13,14,15,16,18,19,20,21,22,23,24,25,27,28,29,30,31,32,33,36,37,38,40)
lids <- length(ids)
ncores <- lids

all_lls <- mclapply(ids,lse_llest,wd=wd,wddata=wddata,rp=0.005,nfc=4,b=2,mc.cores=ncores)

#Create matrix of all likelihoods for each individual
id_lls <- matrix(NA,nrow=lids,ncol=5)
for(j in 1:lids){
  id_lls[j,] <- all_lls[[j]]
}

#Save matrix of log-likelihoods when fitting four models
saveRDS(object = id_lls,file=paste0(wd,"LSE_ID_LLs_Final.rds")) 

