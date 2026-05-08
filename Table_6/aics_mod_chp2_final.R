
#Create function to calculate AICs
aic_lls <- function(j,wdest,ll_chp){
  #Load in estimates generated from fitting respective model
  mles <- readRDS(paste0(wdest,"ests",j,".rds"))
  k <- length(mles)

  #Calculate an output AIC
  aic_chp <- 2*k - 2*ll_chp
  return(aic_chp)
}

#Set working directory for all of the estimates
wdest <- "/home/xieryan/Dissertation1/Data/Mod_CHP2/"

#Load in log-likelihoods
ll_ids <- readRDS(paste0(wdest,"Mod_CHP_LLs0.rds"))

#Calculate AIC values
aic_results <- rep(NA,41)

for(j in 1:41){
  ll_chp <- ll_ids[j]
  aic_results[j] <- aic_lls(j,wdest,ll_chp)
}

#Print AIC results (rounded to nearest tenth)
print(round(aic_results,1))
