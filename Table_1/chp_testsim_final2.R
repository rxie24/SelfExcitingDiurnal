library(Rcpp)
library(parallel)

#Local directory
wd_files <- "/Users/xieryan/Desktop/Dissertation_1/Base Functions/"

source(paste0(wd_files,"base_chp5.R"))
source(paste0(wd_files,"gen_data_chp.R"))

chp_est <- function(i,theta,n=500){
  ltheta <- length(theta)
  thetai <- rep(0,ltheta)
  thetai2 <- 0.9*theta
  
  mu <- theta[1]
  beta <- theta[2]
  sigma <- theta[3]
  
  ncoef <- length(theta[4:ltheta])/2
  
  eta <- theta[4:(4+ncoef-1)]
  gammas <- theta[(4+ncoef):ltheta]
  
  times <- gen_data_chp(n,mu,beta,sigma,eta,gammas)
  mles_chp1 <- chp_mles(thetai,times)
  mles_chp2 <- chp_mles(thetai2,times)
  
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

## USING RESULTS BELOW ##

#COMPARE 0 INTIALIZATION VERSUS 90% TRUE
print("Initialize to 0")
colMeans(est_mat,na.rm=TRUE) #0.5335315  0.1460096 -0.3420262 -0.3192684  0.3078155
apply(est_mat,2,median,na.rm=TRUE) #0.5403891  0.1376150 -0.3443379 -0.3162256  0.3054392
apply(est_mat,2,var,na.rm=TRUE) #0.005848560 0.011770486 0.080850104 0.006027273 0.005174769
apply(est_mat,2,sd,na.rm=TRUE) #0.07647588 0.10849187 0.28434153 0.07763551 0.07193587 
apply(est_mat,2,IQR,na.rm=TRUE) #0.09583277 0.16044432 0.29324209 0.09577362 0.10404845
cor(est_mat[,2],est_mat[,3]) #-0.3190628

print("Initialize to Ninety Percent of True")
colMeans(est_mat2,na.rm=TRUE) #0.49737548  0.14111962  0.09521601 -0.30872575  0.30752896
apply(est_mat2,2,median,na.rm=TRUE) #0.4972288  0.1367605  0.1084518 -0.3015399  0.3059159 
apply(est_mat2,2,var,na.rm=TRUE) #0.009146880 0.009887825 0.067827133 0.006464886 0.004569927
apply(est_mat2,2,sd,na.rm=TRUE) #0.09563932 0.09943754 0.26043643 0.08040452 0.06760124
apply(est_mat2,2,IQR,na.rm=TRUE) #0.12035890 0.13329551 0.09791799 0.10597130 0.09006373
cor(est_mat2[,2],est_mat2[,3]) #-0.1552322
