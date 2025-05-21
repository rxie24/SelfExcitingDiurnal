
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

#Set working directory for all of the CHP likelihoods (MODIFY TO BE YOUR SPECIFIC DIRECTORY)
wdest <- "/home/xieryan/Dissertation1/Data/CHP_GH/"

#Aggregrate saved likelihoods together into matrix

all_lls <- matrix(NA,nrow=41,ncol=4)

#CHP_ESTS_1
ll_ids1 <- readRDS(paste0(wdest,"CHP_ID_LLs_Final.rds"))
id_ind1 <- c(1,2,3,4,5,6,8,9,11,12,14,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,36,37,39,40)

all_lls[id_ind1,1] <- as.numeric(ll_ids1[,2])
all_lls[id_ind1,2] <- as.numeric(ll_ids1[,3])
all_lls[id_ind1,3] <- as.numeric(ll_ids1[,4])
all_lls[id_ind1,4] <- as.numeric(ll_ids1[,5])

#CHP_ESTS_2
ll_ids2 <- readRDS(paste0(wdest,"CHP_ID_LLs_Final2.rds"))
id_ind2 <- c(7,34,35,38,41)
all_lls[id_ind2,1] <- as.numeric(ll_ids2[,2])
all_lls[id_ind2,2] <- as.numeric(ll_ids2[,3])
all_lls[id_ind2,3] <- as.numeric(ll_ids2[,4])
all_lls[id_ind2,4] <- as.numeric(ll_ids2[,5])

#CHP_ESTS_3
ll_ids3 <- readRDS(paste0(wdest,"CHP_ID_LLs_Final3.rds"))
id_ind3 <- c(10,13,17)
all_lls[id_ind3,1] <- as.numeric(ll_ids3[,2])
all_lls[id_ind3,2] <- as.numeric(ll_ids3[,3])
all_lls[id_ind3,3] <- as.numeric(ll_ids3[,4])
all_lls[id_ind3,4] <- as.numeric(ll_ids3[,5])

#Calculate AIC values
aic_results <- matrix(NA,nrow=41,ncol=4)
colnames(aic_results) <- c("PP","CP","HP","CHP")

for(j in 1:41){
  id_lls <- all_lls[j,]
  aic_results[j,] <- aic_lls(j,wdest,id_lls)
}

#Print AIC results (rounded to nearest tenth)
print(round(aic_results,1))
