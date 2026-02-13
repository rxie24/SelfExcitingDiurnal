library(Rcpp)
library(parallel)

#Directory for base function files
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

#Load in base function files
source(paste0(wd_files,"base_lse4.R"))
source(paste0(wd_files,"lse_complls.R"))

lse_llest <- function(k,wd,wddata,nfc,b=2){
  
  #Load in the data
  times <- readRDS(paste0(wddata,"data",k,".rds"))
  
  #Fit Circadian Process
  nfc_k <- nfc[k]
  thetai_cp <- rep(0,(1+nfc_k))
  mles_cp <- cp_mles2(thetai_cp,times,b=b)
  ll_cp <- -1*cp_ll2(mles_cp,times,b=b)
  mu_init <- mles_cp[1]
  cp_init <- mles_cp[2:length(mles_cp)]
  
  #Fit CLSE model
  thetai_lse <- c(mu_init,rep(0,2),cp_init)
  mles_lse <- lse_mles(thetai_lse,times,b=b)
  
  #Save mles of full linear self-exciting model
  saveRDS(object = mles_lse,file=paste0(wd,"ests",k,".rds"))
  ll_lse <- -1*lse_ll(mles_lse,times,b=b)
  
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
wd <- "/home/xieryan/Dissertation1/Data/CLSE2_Final/"
wddata <- "/home/xieryan/Dissertation1/Data/Real_Data_Event_Times_V2_GitHub/"
nfc <- readRDS("/home/xieryan/Dissertation1/Data/CLSE2/clse_lengths.rds")

ids <- seq(1,41,1) #23
lids <- length(ids)
ncores <- lids

all_lls <- mclapply(ids,lse_llest,wd=wd,wddata=wddata,nfc=nfc,b=2,mc.cores=ncores)
print(all_lls)

#Create matrix of all likelihoods for each individual
id_lls <- matrix(NA,nrow=lids,ncol=5)
for(j in 1:lids){
  id_lls[j,] <- all_lls[[j]]
}

#Save matrix of log-likelihoods when fitting four models
saveRDS(object = id_lls,file=paste0(wd,"LSE_ID_LLs0.rds")) 
