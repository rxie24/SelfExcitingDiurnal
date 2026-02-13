library(Rcpp)
library(parallel)

#Directories for base function files and to save results
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/"
wd <- "/home/xieryan/Dissertation1/Simulations2/"

#Load in needed base function files
source(paste0(wd_files,"base_chp5.R"))
source(paste0(wd_files,"base_lse4.R"))
source(paste0(wd_files,"gen_data_chp.R"))
source(paste0(wd_files,"chp_complls.R"))
source(paste0(wd_files,"models_five.R"))
source(paste0(wd_files,"chp_five.R"))
source(paste0(wd_files,"chp_five_sim.R"))
source(paste0(wd_files,"filt_outliers.R"))

RNGkind("L'Ecuyer-CMRG")
set.seed(215)

mu <- 0.1
beta_seq <- seq(0.15,0.5,0.025)
eta <- -0.15
gammas <- 0.15

thetai_lse <- rep(0,5)
thetai_chp <- rep(0,5)

navg <- 300 #number of iterations to take average over
nocores <- 60

#Subplot 1
print("Subplot 1")

#Log(Sigma) = 0.2, N = 50
sigma <- 0.2 #Vary sigma between different subfigures
theta_wb <- c(mu,sigma,eta,gammas)
N <- 50 #Number of waiting times

five_testlls1 <- test_lls(beta_seq,navg,theta_wb,2,thetai_chp,thetai_lse,N,ncores=nocores)
saveRDS(five_testlls1,file=paste0(wd,"five_testlls1_2.rds"))

#Subplot 2
print("Subplot 2")

#Log(Sigma) = 0.2, N = 500
sigma <- 0.2 #Vary sigma between different subfigures
theta_wb <- c(mu,sigma,eta,gammas)
N <- 500

five_testlls2 <- test_lls(beta_seq,navg,theta_wb,2,thetai_chp,thetai_lse,N,ncores=nocores)
saveRDS(five_testlls2,file=paste0(wd,"five_testlls2_2.rds"))

#Subplot 3
print("Subplot 3")

#Log(Sigma) = -2, N = 50
sigma <- -2 #Vary sigma between different subfigures
theta_wb <- c(mu,sigma,eta,gammas)
N <- 50

five_testlls3 <- test_lls(beta_seq,navg,theta_wb,2,thetai_chp,thetai_lse,N,ncores=nocores)
saveRDS(five_testlls3,file=paste0(wd,"five_testlls3_2.rds"))

#Subplot 4
print("Subplot 4")

#Log(Sigma) = -2, N = 500
sigma <- -2 #Vary sigma between different subfigures
theta_wb <- c(mu,sigma,eta,gammas)
N <- 500

five_testlls4 <- test_lls(beta_seq,navg,theta_wb,2,thetai_chp,thetai_lse,N,ncores=nocores)
saveRDS(five_testlls4,file=paste0(wd,"five_testlls4_2.rds"))
