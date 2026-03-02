library(Rcpp)

#Set directories for PDF file, data, and base functions
wd_files <- "/Users/xieryan/Desktop/Dissertation_1/Base Functions/"
wdest <- "/Users/xieryan/Desktop/Dissertation_1/Real_Data2/Baseline_Hawkes/"
wddata <- '/Users/xieryan/Desktop/Dissertation_1/Real Data/Real_Data_Event_Times_V2_GitHub/'

#Local directory
source(paste0(wd_files,"base_chp5.R"))
source(paste0(wd_files,"chp_complls.R"))
source(paste0(wd_files,"gen_data_chp.R"))

set.seed(100)

#LOAD IN EXAMPLE DATASET
times <- readRDS(paste0(wddata,"data11.rds"))

ltimes <- length(times)-1

#Fit Hawkes Model
thetai_hp <- rep(0,3)
mles_hp <- hp_mles(thetai_hp,times)
 
#Simulate data according to Hawkes Process
mu <- mles_hp[1]
beta <- mles_hp[2]
sigma <- mles_hp[3]
eta <- 0
gammas <- 0

hp_data <- gen_data_chp(ltimes,mu,beta,sigma,eta,gammas)

#Plot Figure 1 (to prove Hawkes Process is insufficient)
hist_breaks <- seq(0,24,1)

pdf(file=paste0(wdest,"base_hawkes_example3.pdf"),width=5,height=5) #or 4x4
par(mgp = c(2.5, 1, 0), mar = c(4, 3.5, 0.1, 0.25))
hist(times%%24,prob=TRUE,
     #breaks=hist_breaks, main="Individual Data Compared to Hawkes", 
     breaks=hist_breaks, main="",
     axes=FALSE, ylim=c(0,0.1),xlim=c(0,24),
     xlab="Time of Day (Hrs)")
lines(density(hp_data%%24),col="blue")

axis(2, at = seq(0, 0.1, by = 0.05), labels = seq(0, 0.1, by = 0.05))
axis(1, at = seq(0, 24, by = 12), labels = seq(0, 24, by = 12))
dev.off()
