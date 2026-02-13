

#Create function to calculate AICs
aic_lls <- function(j,wdest,id_lls){
  #Load in estimates generated from fitting respective model
  mles <- readRDS(paste0(wdest,"ests",j,".rds"))
  k <- length(mles)
  nfoco <- length(mles[4:k])
  cp_k <- nfoco + 1
  
  ll_pp <- id_lls[1]
  ll_cp <- id_lls[2]
  ll_hp <- id_lls[3]
  ll_chp <- id_lls[4]
  
  aic_pp <- 2 - 2*ll_pp
  aic_cp <- 2*cp_k - 2*ll_cp
  aic_hp <- 2*3 - 2*ll_hp
  aic_chp <- 2*k - 2*ll_chp
  
  #Output test statistics with comparison
  aic_ids <- c(aic_pp,aic_cp,aic_hp,aic_chp)
  return(aic_ids)
}

#Set working directory for all of the estimates
wdest <- "/home/xieryan/Dissertation1/Data/CLSE2_Final/"

#Aggregrate likelihoods together
all_lls <- matrix(NA,nrow=41,ncol=4)

ll_ids <- readRDS(paste0(wdest,"LSE_ID_LLs0.rds"))

all_lls[,1] <- as.numeric(ll_ids[,2])
all_lls[,2] <- as.numeric(ll_ids[,3])
all_lls[,3] <- as.numeric(ll_ids[,4])
all_lls[,4] <- as.numeric(ll_ids[,5])

#Calculate AICs
aic_results <- matrix(NA,nrow=41,ncol=4)
colnames(aic_results) <- c("PP","CP","LSE","CLSE")

for(j in 1:41){
  id_lls <- all_lls[j,]
  aic_results[j,] <- aic_lls(j,wdest,id_lls)
}

#Print AIC results (rounded to nearest tenth)
print(round(aic_results,1))
