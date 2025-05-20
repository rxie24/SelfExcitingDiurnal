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
  
  #Fit Linear Self-Exciting Model w/o Circadian Component
  thetai_nc <- rep(0,3)
  mles_nc <- lse_nc_mles(thetai_nc,times,b=b)
  ll_nc <- -1*lse_nc_ll(mles_nc,times,b=b)
  
  #Fit full LSE Model to the data w/ two different sets of initial values
  thetai_lse1 <- rep(0,5)
  mles_lse1 <- lse_mles(thetai_lse1,times,b=b)

  thetai_lse2 <- c(mles_nc,rep(0,2))
  mles_lse2 <- lse_mles(thetai_lse2,times,b=b)
  
  #Calculate total sum of squares (across cumulative and max) for two different initial values
  ml_fs01 <- mles_lse1[4:5]
  sscum1 <- cum_ss(times,ml_fs01)
  ropen1 <- rpen(ml_fs01)
  ropen1 <- rp*ropen1
  ro_pe0cum1 <- sscum1 + ropen1

  ml_fs02 <- mles_lse2[4:5]
  sscum2 <- cum_ss(times,ml_fs02)
  ropen2 <- rpen(ml_fs02)
  ropen2 <- rp*ropen2
  ro_pe0cum2 <- sscum2 + ropen2
  
  if(ro_pe0cum2 < ro_pe0cum1){
    ro_pe0cum <- ro_pe0cum2
    mles_lse <- mles_lse2
  }else{
    ro_pe0cum <- ro_pe0cum1
    mles_lse <- mles_lse1
  }
  
  for(j in 2:nfc){
    #1st initial values
    thetai_lse1 <- rep(0,(3+(2*j)))
    par_estj1 <- lse_mles(thetai_lse1,times,b=b)
    
    ml_fsj1 <- par_estj1[4:length(par_estj1)]
    sscumj1 <- cum_ss(times,ml_fsj1)
    ropenj1 <- rpen(ml_fsj1)
    ropenj1 <- rp*ropenj1
    ro_pejcum1 <- sscumj1 + ropenj1
    
    #2nd initial values
    thetai_lse2 <- c(mles_nc,rep(0,(2*j)))
    par_estj2 <- lse_mles(thetai_lse2,times,b=b)
    
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
  
  #Fit Poisson Process
  ll_pp <- max_logpp2(times,b=b)
  
  #Save 
  ll_est <- c(k,ll_pp,ll_cp,ll_nc,ll_lse) #Poisson, Circadian, LSE w/o FS, LSE
  
  return(ll_est)
}

#Establish working directories to save parameter estimates and where data is located 
wd <- "/home/xieryan/Dissertation1/Data/LSE_GH/"
wddata <- "/home/xieryan/Dissertation1/Data/Real_Data_Event_Times_V2_GitHub/"

ids <- c(7,17,35)

all_lls <- mclapply(ids,lse_llest,wd=wd,wddata=wddata,rp=0,nfc=4,b=2,mc.cores=3)

#Create matrix of all likelihoods for each individual
id_lls <- matrix(NA,nrow=3,ncol=5)
for(j in 1:3){
  id_lls[j,] <- all_lls[[j]]
}

#Save matrix of log-likelihoods when fitting four models
saveRDS(object = id_lls,file=paste0(wd,"LSE_ID_LLs_Final5.rds")) 

