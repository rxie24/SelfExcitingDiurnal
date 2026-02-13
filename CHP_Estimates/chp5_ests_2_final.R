library(Rcpp)
library(parallel)

#Directory for base function files
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

source(paste0(wd_files,"base_chp5.R"))
source(paste0(wd_files,"chp_complls.R"))
source(paste0(wd_files,"chp_llest_base.R"))

#Establish working directories to save parameter estimates and where data is located 
wd <- "/home/xieryan/Dissertation1/Data/CHP5_Final/"
wddata <- "/home/xieryan/Dissertation1/Data/Real_Data_Event_Times_V2_GitHub/"

#Load in number of Fourier coefficients used for each individual
nfc <- readRDS("/home/xieryan/Dissertation1/Data/CHP4/chp_lengths.rds")

#Subset individuals to be used
ids <- c(13, 14, 26, 41)
ncores <- length(ids)

#Fit CHP Models for chosen subset of individuals
all_lls <- mclapply(ids,chp_llest,wd=wd,wddata=wddata,nfc=nfc,ip=0.8,mc.cores=ncores) #have no roughness penalty
print(all_lls)

#Create matrix of all likelihoods for each individual
id_lls <- matrix(NA,nrow=ncores,ncol=5)
for(j in 1:ncores){
  id_lls[j,] <- all_lls[[j]]
}

#Save matrix of log-likelihoods when fitting four models
saveRDS(object = id_lls,file=paste0(wd,"CHP_ID_LLs2.rds")) 

