library(Rcpp)
library(parallel)

#Local directory
wd_files <- "/Users/xieryan/Desktop/Dissertation_1/Base Functions/"

source(paste0(wd_files,"base_lse4.R"))
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

theta <- c(0.2,0.01,0.01,-0.3,0.3)
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
colMeans(est_mat,na.rm=TRUE) #0.15514266  0.01024336  0.01051211 -0.33008705  0.28502750
apply(est_mat,2,median,na.rm=TRUE) #0.14669305  0.01042750  0.01051953 -0.32593226  0.27592848
apply(est_mat,2,var,na.rm=TRUE) #2.104990e-03 4.126814e-06 5.054978e-06 5.087406e-03 5.837201e-03
apply(est_mat,2,sd,na.rm=TRUE) #0.045880174 0.002031456 0.002248328 0.071326057 0.076401576
apply(est_mat,2,IQR,na.rm=TRUE) #0.057753473 0.002906906 0.002373687 0.089465443 0.097191687
cor(est_mat[,2],est_mat[,3]) #-0.02115539

#Look at results when values initialized to 90% of true
print("Initialize to Ninety Percent of True")
colMeans(est_mat2,na.rm=TRUE) #0.201780170  0.009744320  0.009973239 -0.299464840  0.291955269
apply(est_mat2,2,median,na.rm=TRUE) #0.204625807  0.009550994  0.010062085 -0.293370942  0.289559705
apply(est_mat2,2,var,na.rm=TRUE) #1.315499e-03 4.534030e-06 3.447100e-06 3.072600e-03 2.823884e-03
apply(est_mat2,2,sd,na.rm=TRUE) #0.036269809 0.002129326 0.001856637 0.055431039 0.053140234
apply(est_mat2,2,IQR,na.rm=TRUE) #0.039488656 0.002687011 0.002241733 0.057547461 0.060051928
cor(est_mat2[,2],est_mat2[,3]) #-0.2002988
