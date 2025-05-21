library(Rcpp)
library(parallel)

#Directory for base function files (MODIFY TO BE YOUR SPECIFIC DIRECTORY)
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

#Load in base function files
source(paste0(wd_files,"base_lse.R"))
source(paste0(wd_files,"lse_complls.R"))
source(paste0(wd_files,"rope_functions.R"))

lse_llest <- function(k,wd,wddata,b=2){
  
  #Load in the data
  times <- readRDS(paste0(wddata,"data",k,".rds"))
  
  #Fit Linear Self-Exciting Model w/o Circadian Component
  thetai_nc <- rep(0,3)
  mles_nc <- lse_nc_mles(thetai_nc,times,b=b)
  ll_nc <- -1*lse_nc_ll(mles_nc,times,b=b)
  
  #Fit full LSE Model to the data w/ two different sets of initial values
  thetai_lse <- rep(0,9) #For these individuals, fix the number of coefficients to estimate on
  mles_lse <- lse_mles(thetai_lse,times,b=b)
  saveRDS(object = mles_lse,file=paste0(wd,"ests",k,".rds"))  ##Save mles of full linear self-exciting model
  ll_lse <- -1*lse_ll(mles_lse,times,b=b)
  
  #Fit Circadian Process
  lfc <- length(mles_lse[4:length(mles_lse)])
  thetai_cp <- rep(0,(1+lfc))
  mles_cp <- cp_mles2(thetai_cp,times,b=b)
  ll_cp <- -1*cp_ll2(mles_cp,times,b=b)
  
  #Fit Poisson Process
  ll_pp <- max_logpp2(times,b=b)
  
  #Save 
  ll_est <- c(k,ll_pp,ll_cp,ll_nc,ll_lse) #Poisson, Circadian, LSE w/o FS, LSE
  
  return(ll_est)
}

#Establish working directories to saved CLSE parameter estimates and where data is located (MODIFY TO BE YOUR SPECIFIC DIRECTORY)
wd <- "/home/xieryan/Dissertation1/Data/LSE_GH/"
wddata <- "/home/xieryan/Dissertation1/Data/Real_Data_Event_Times_V2_GitHub/"

ids <- c(10,34)

all_lls <- mclapply(ids,lse_llest,wd=wd,wddata=wddata,b=2,mc.cores=2)

#Create matrix of all likelihoods for each individual
id_lls <- matrix(NA,nrow=2,ncol=5)
for(j in 1:2){
  id_lls[j,] <- all_lls[[j]]
}

#Save matrix of log-likelihoods when fitting four models
saveRDS(object = id_lls,file=paste0(wd,"LSE_ID_LLs_Final6.rds")) 

