#Establish directory of saved CLSE likelihoods (MODIFY TO BE YOUR SPECIFIC DIRECTORY)
wdest <- "/home/xieryan/Dissertation1/Data/LSE_GH/"

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

#Aggregrate saved likelihoods together into matrix

all_lls <- matrix(NA,nrow=41,ncol=4)

#LSE_ESTS_1
ll_ids1 <- readRDS(paste0(wdest,"LSE_ID_LLs_Final.rds"))
id_ind1 <- c(2,3,4,5,8,9,11,12,13,14,15,16,18,19,20,21,22,23,24,25,27,28,29,30,31,32,33,36,37,38,40)

all_lls[id_ind1,1] <- as.numeric(ll_ids1[,2])
all_lls[id_ind1,2] <- as.numeric(ll_ids1[,3])
all_lls[id_ind1,3] <- as.numeric(ll_ids1[,4])
all_lls[id_ind1,4] <- as.numeric(ll_ids1[,5])

#LSE_ESTS_2
ll_ids2 <- readRDS(paste0(wdest,"LSE_ID_LLs_Final2.rds"))
id_ind2 <- c(1,26,39)
all_lls[id_ind2,1] <- as.numeric(ll_ids2[,2])
all_lls[id_ind2,2] <- as.numeric(ll_ids2[,3])
all_lls[id_ind2,3] <- as.numeric(ll_ids2[,4])
all_lls[id_ind2,4] <- as.numeric(ll_ids2[,5])

#LSE_ESTS_3
ll_ids3 <- readRDS(paste0(wdest,"LSE_ID_LLs_Final3.rds"))
all_lls[41,] <- ll_ids3[2:5]

#LSE_ESTS_4
ll_ids4 <- readRDS(paste0(wdest,"LSE_ID_LLs_Final4.rds"))
all_lls[6,] <- ll_ids4[2:5]

#LSE_ESTS_5
ll_ids5 <- readRDS(paste0(wdest,"LSE_ID_LLs_Final5.rds"))
id_ind5 <- c(7,17,35)
all_lls[id_ind5,1] <- as.numeric(ll_ids5[,2])
all_lls[id_ind5,2] <- as.numeric(ll_ids5[,3])
all_lls[id_ind5,3] <- as.numeric(ll_ids5[,4])
all_lls[id_ind5,4] <- as.numeric(ll_ids5[,5])

#LSE_ESTS_6
ll_ids6 <- readRDS(paste0(wdest,"LSE_ID_LLs_Final6.rds"))
id_ind6 <- c(10,34)
all_lls[id_ind6,1] <- as.numeric(ll_ids6[,2])
all_lls[id_ind6,2] <- as.numeric(ll_ids6[,3])
all_lls[id_ind6,3] <- as.numeric(ll_ids6[,4])
all_lls[id_ind6,4] <- as.numeric(ll_ids6[,5])

#Calculate AICs
aic_results <- matrix(NA,nrow=41,ncol=4)
colnames(aic_results) <- c("PP","CP","LSE","CLSE")

for(j in 1:41){
  id_lls <- all_lls[j,]
  aic_results[j,] <- aic_lls(j,wdest,id_lls)
}

#Print AIC results (rounded to nearest tenth)
print(round(aic_results,1))
