library(Rcpp)
library(parallel)

#Directory for base function files (MODIFY TO BE YOUR SPECIFIC DIRECTORY)
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster

#Load in needed base function files
source(paste0(wd_files,"base_chp.R"))
source(paste0(wd_files,"base_lse.R"))
source(paste0(wd_files,"gen_data_chp.R"))
source(paste0(wd_files,"chp_complls.R"))
source(paste0(wd_files,"models_five.R"))
source(paste0(wd_files,"chp_five.R"))
source(paste0(wd_files,"chp_five_sim.R"))
source(paste0(wd_files,"filt_outliers.R"))

RNGkind("L'Ecuyer-CMRG")
set.seed(215)

mu <- 0.1
beta_seq <- seq(0.15,0.25,0.025)
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
chp_testlls1 <- five_testlls1[[1]]
hp_testlls1 <- five_testlls1[[2]]
pp_testlls1 <- five_testlls1[[3]]
cp_testlls1 <- five_testlls1[[4]]
lse_testlls1 <- five_testlls1[[5]]

#Subplot 2
print("Subplot 2")

#Log(Sigma) = 0.2, N = 500
sigma <- 0.2 #Vary sigma between different subfigures
theta_wb <- c(mu,sigma,eta,gammas)
N <- 500

five_testlls2 <- test_lls(beta_seq,navg,theta_wb,2,thetai_chp,thetai_lse,N,ncores=nocores)
chp_testlls2 <- five_testlls2[[1]]
hp_testlls2 <- five_testlls2[[2]]
pp_testlls2 <- five_testlls2[[3]]
cp_testlls2 <- five_testlls2[[4]]
lse_testlls2 <- five_testlls2[[5]]

#Subplot 3
print("Subplot 3")

#Log(Sigma) = -2, N = 50
sigma <- -2 #Vary sigma between different subfigures
theta_wb <- c(mu,sigma,eta,gammas)
N <- 50

five_testlls3 <- test_lls(beta_seq,navg,theta_wb,2,thetai_chp,thetai_lse,N,ncores=nocores)
chp_testlls3 <- five_testlls3[[1]]
hp_testlls3 <- five_testlls3[[2]]
pp_testlls3 <- five_testlls3[[3]]
cp_testlls3 <- five_testlls3[[4]]
lse_testlls3 <- five_testlls3[[5]]

#Subplot 4
print("Subplot 4")

#Log(Sigma) = -2, N = 500
sigma <- -2 #Vary sigma between different subfigures
theta_wb <- c(mu,sigma,eta,gammas)
N <- 500

five_testlls4 <- test_lls(beta_seq,navg,theta_wb,2,thetai_chp,thetai_lse,N,ncores=nocores)
chp_testlls4 <- five_testlls4[[1]]
hp_testlls4 <- five_testlls4[[2]]
pp_testlls4 <- five_testlls4[[3]]
cp_testlls4 <- five_testlls4[[4]]
lse_testlls4 <- five_testlls4[[5]]

#Find limits for plots (test)
lims1 <- unlist(five_testlls1)
lims2 <- unlist(five_testlls2)
lims3 <- unlist(five_testlls3)
lims4 <- unlist(five_testlls4)

miny1 <- min(lims1,lims2)
maxy1 <- max(lims1,lims2)

miny2 <- min(lims3,lims4)
maxy2 <- max(lims3,lims4)

#Create plots
#Establish working directory to save resulting figure (MODIFY TO BE YOUR SPECIFIC DIRECTORY) 
wd <- "/home/xieryan/Dissertation1/Simulations/"

#Plot 1: Trajectories of all five model comparisons
pdf(file=paste0(wd,"five_comp_final.pdf"), width = 5, height = 5)
par(mfrow = c(2,2), mar = c(3,3,1.5,0.5), oma = c(2, 2, 2.5, 0))

plot(beta_seq,chp_testlls1,type="l",col="red",
     xlab="",ylab="",main=expression(bold(paste("N=50, ", phi, " = 0.2"))),
     ylim=c(miny1,maxy1))
lines(beta_seq,hp_testlls1,col="blue")
lines(beta_seq,pp_testlls1,col="purple")
lines(beta_seq,cp_testlls1,col="orange")
lines(beta_seq,lse_testlls1,col="green")
axis(2, at = seq(0, 0.2, by = 0.1), labels = seq(0, 0.2, by = 0.1))  # Custom Y-ticks

plot(beta_seq,chp_testlls2,type="l",col="red",
     xlab="",ylab="",main=expression(bold(paste("N=500, ", phi, " = 0.2"))),
     ylim=c(miny1,maxy1))
lines(beta_seq,hp_testlls2,col="blue")
lines(beta_seq,pp_testlls2,col="purple")
lines(beta_seq,cp_testlls2,col="orange")
lines(beta_seq,lse_testlls2,col="green")

plot(beta_seq,chp_testlls3,type="l",col="red",
     xlab="",ylab="",main=expression(bold(paste("N=50, ", phi, " = -2"))),
     ylim=c(miny2,maxy2))
lines(beta_seq,hp_testlls3,col="blue")
lines(beta_seq,pp_testlls3,col="purple")
lines(beta_seq,cp_testlls3,col="orange")
lines(beta_seq,lse_testlls3,col="green")

plot(beta_seq,chp_testlls4,type="l",col="red",
     xlab="",ylab="",main=expression(bold(paste("N=500, ", phi, " = -2"))),
     ylim=c(miny2,maxy2))
lines(beta_seq,hp_testlls4,col="blue")
lines(beta_seq,pp_testlls4,col="purple")
lines(beta_seq,cp_testlls4,col="orange")
lines(beta_seq,lse_testlls4,col="green")

legend("bottomright",legend=c("CHP","Hawkes","Poisson","Circadian","CLSE"),
       fill=c("red","blue","purple","orange","green"),cex = 0.7)


mtext("Comparison of Test Log-Likelihoods", side = 3, outer = TRUE, line = 0.5, cex = 1.4, font = 2)
mtext(bquote(beta), side = 1, outer = TRUE, line = 0.2, cex = 1)
mtext("Test Log-Likelihoods", side = 2, outer = TRUE, line = 0.2, cex = 1)
dev.off()

