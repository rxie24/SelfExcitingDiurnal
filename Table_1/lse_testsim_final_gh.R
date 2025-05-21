library(Rcpp)
library(parallel)

#Directory for base function files (MODIFY TO BE YOUR SPECIFIC DIRECTORY)
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

#Load in needed base function files
source(paste0(wd_files,"base_lse.R"))
source(paste0(wd_files,"gen_data_lse2.R"))

#Create function to get CLSE estimates from CLSE simulated data w/ two initial values: all 0 and ninety percent of true
lse_est <- function(i,theta,b=2,n=500){
  #b: order of beta (number of points to go back to)
  e_ind <- 2 + b
  ltheta <- length(theta)
  thetai <- rep(0,ltheta) #first initial values (all 0)
  thetai2 <- 0.9*theta #second initial values (ninety percent of true)
  
  mu <- theta[1]
  beta <- theta[2:(e_ind-1)]
  
  ncoef <- length(theta[e_ind:ltheta])/2
  eta <- theta[e_ind:(e_ind+ncoef-1)]
  gammas <- theta[(e_ind+ncoef):ltheta]
  
  times <- gen_data_lse(n,mu,beta,eta,gammas) #simulate CLSE data
  mles_lse1 <- lse_mles(thetai,times) #get CLSE estimates under first initial values (all 0)
  mles_lse2 <- lse_mles(thetai2,times) #get CLSE estimates under second initial values (ninety percent of true)
  
  mles_lse <- rbind(mles_lse1,mles_lse2)
  return(mles_lse)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(215)

theta <- c(0.5,0.2,0.1,-0.3,0.3)
est_mat <- matrix(NA,nrow=100,ncol=5)
est_mat2 <- matrix(NA,nrow=100,ncol=5)

est_parallel <- mclapply(1:100,lse_est,theta=theta,mc.cores=10) #10 cores locally

for(j in 1:100){
  estj <- est_parallel[[j]]
  est_mat[j,] <- estj[1,]
  est_mat2[j,] <- estj[2,]
}


#COMPARE 0 INTIALIZATION VERSUS 90% TRUE
print("Initialize to 0")
colMeans(est_mat,na.rm=TRUE)
apply(est_mat,2,median,na.rm=TRUE)
apply(est_mat,2,var,na.rm=TRUE)
apply(est_mat,2,IQR,na.rm=TRUE)
cor(est_mat[,2],est_mat[,3])
#SDs = square root of variances

#Look at results when values initialized to 90% of true
print("Initialize to Ninety Percent of True")
colMeans(est_mat2,na.rm=TRUE)
apply(est_mat2,2,median,na.rm=TRUE)
apply(est_mat2,2,var,na.rm=TRUE)
apply(est_mat2,2,IQR,na.rm=TRUE) 
cor(est_mat2[,2],est_mat2[,3]) 
#SDs = square root of variances
