library(Rcpp)
library(parallel)

#Directory for base function files
#wd_files <- "/Users/xieryan/Desktop/Dissertation_1/Base Functions/" #Local
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

#Load in needed base function files
source(paste0(wd_files,"base_lse.R"))
source(paste0(wd_files,"gen_data_lse2.R"))

lse_est <- function(i,theta,b=2,n=500){
  #b: order of beta (number of points to go back to)
  e_ind <- 2 + b
  ltheta <- length(theta)
  thetai <- rep(0,ltheta)
  thetai2 <- 0.9*theta
  
  mu <- theta[1]
  beta <- theta[2:(e_ind-1)]
  
  ncoef <- length(theta[e_ind:ltheta])/2
  eta <- theta[e_ind:(e_ind+ncoef-1)]
  gammas <- theta[(e_ind+ncoef):ltheta]
  
  times <- gen_data_lse(n,mu,beta,eta,gammas)
  mles_lse1 <- lse_mles(thetai,times)
  mles_lse2 <- lse_mles(thetai2,times)
  
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


## USING RESULTS BELOW ##

#COMPARE 0 INTIALIZATION VERSUS 90% TRUE
print("Initialize to 0")
colMeans(est_mat,na.rm=TRUE) #0.4895437  0.2090635  0.1009286 -0.3086907  0.2961485
apply(est_mat,2,median,na.rm=TRUE) #0.4937076  0.2131867  0.1042848 -0.3048419  0.2946589
apply(est_mat,2,var,na.rm=TRUE) #0.008287688 0.010119920 0.009161813 0.005941252 0.005243471
apply(est_mat,2,IQR,na.rm=TRUE) #
cor(est_mat[,2],est_mat[,3]) #0.1251023
#SDs: 0.09103674 0.10059781 0.09571736 0.07707952 0.07241182 (square root of variances)

#Look at results when values initialized to 90% of true
print("Initialize to Ninety Percent of True")
colMeans(est_mat2,na.rm=TRUE) #0.4764045  0.2225372  0.1128101 -0.3086346  0.3021806
apply(est_mat2,2,median,na.rm=TRUE) #0.4860190  0.2210120  0.1108207 -0.3050475  0.3048771
apply(est_mat2,2,var,na.rm=TRUE) #0.007707801 0.008814560 0.007406082 0.005735614 0.005501466
apply(est_mat2,2,IQR,na.rm=TRUE) #
cor(est_mat2[,2],est_mat2[,3]) #0.1415362
#SDs: 0.08779408 0.09388589 0.08605860 0.07573384 0.07417187 (square root of variances)
