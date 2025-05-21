library(Rcpp)
library(parallel)

#Directory for base function files (MODIFY TO BE YOUR SPECIFIC DIRECTORY)
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

#Load in needed base function files
source(paste0(wd_files,"base_chp.R"))
source(paste0(wd_files,"gen_data_chp.R"))

#Create function to get CHP estimates from CHP simulated data w/ two sets of initial values
chp_est <- function(i,theta,n=500){
  ltheta <- length(theta)
  thetai <- rep(0,ltheta) #first initial values: all 0
  thetai2 <- 0.9*theta #second initial values: ninety percent of true
  
  mu <- theta[1]
  beta <- theta[2]
  sigma <- theta[3]
  
  ncoef <- length(theta[4:ltheta])/2
  
  eta <- theta[4:(4+ncoef-1)]
  gammas <- theta[(4+ncoef):ltheta]
  
  times <- gen_data_chp(n,mu,beta,sigma,eta,gammas) #simulate CHP data
  mles_chp1 <- chp_mles(thetai,times) #get estimates under first initial values (all 0)
  mles_chp2 <- chp_mles(thetai2,times) #get estimates under second initial values (ninety percent of true)
  
  mles_chp <- rbind(mles_chp1,mles_chp2)
  return(mles_chp)
}


RNGkind("L'Ecuyer-CMRG")
set.seed(215)

theta <- c(0.5,0.15,0.1,-0.3,0.3)
est_mat <- matrix(NA,nrow=100,ncol=5)
est_mat2 <- matrix(NA,nrow=100,ncol=5)

est_parallel <- mclapply(1:100,chp_est,theta=theta,mc.cores=10)

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
apply(est_mat,2,sd,na.rm=TRUE)
apply(est_mat,2,IQR,na.rm=TRUE)
cor(est_mat[,2],est_mat[,3])

print("Initialize to Ninety Percent of True")
colMeans(est_mat2,na.rm=TRUE)
apply(est_mat2,2,median,na.rm=TRUE)
apply(est_mat2,2,var,na.rm=TRUE)
apply(est_mat2,2,sd,na.rm=TRUE)
apply(est_mat2,2,IQR,na.rm=TRUE)
cor(est_mat2[,2],est_mat2[,3])
