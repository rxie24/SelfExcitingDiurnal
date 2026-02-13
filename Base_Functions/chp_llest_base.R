
chp_llest <- function(k,wd,wddata,nfc,ip){
  
  #Load in the data
  times <- readRDS(paste0(wddata,"data",k,".rds"))
  
  #Fit Hawkes Model
  thetai_hp <- rep(0,3)
  mles_hp <- hp_mles(thetai_hp,times)
  ll_hp <- -1*hp_ll(mles_hp,times)
  
  #Fit Circadian Hawkes Model w/ Hawkes parameters as initial values
  nfc_k <- nfc[k]
  thetai_chp <- ip*c(mles_hp,rep(0,nfc_k))
  mles_chp <- chp_mles(thetai_chp,times)
  
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
