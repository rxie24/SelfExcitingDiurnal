#Script encompassing code to generate histogram figures

#Directory for base function files
#wd_files <- "/Users/xieryan/Desktop/Dissertation_1/Base Functions/" #Local
wd_files <- "/home/xieryan/Dissertation1/Base_Functions/" #Cluster
source(paste0(wd_files,"rope_functions.R"))

cum_hist <- function(j,wdest,wddata){
  #id: individual
  #wdest: working directory of the estimates
  
  #Load in the data
  times <- readRDS(paste0(wddata,"data",j,".rds"))
  tmod24 <- times %% 24
  
  #Specific bin midpoints and breaks for histograms
  hist_mids <- seq(0.5,23.5,1)
  hist_breaks <- seq(0,24,1)
  
  #Load in estimates generated from fitting respective model
  mles <- readRDS(paste0(wdest,"ests",j,".rds"))
  lmles <- length(mles)
  ncoef <- length(mles[4:lmles])/2
  ml_eta <- as.double(mles[4:(4+ncoef-1)])
  ml_gam <- as.double(mles[(4+ncoef):lmles])
  ml_fs <- c(ml_eta,ml_gam)
  
  #Calculate scaled Fourier series for plots
  x <- seq(0,24,0.1)
  fs <- fourier(x,ml_gam,ml_eta)
  scal_fs <- (fs - min(fs))/(max(fs) - min(fs))
  
  #CREATE HISTOGRAM LOOKING AT ALL EVENT TIMES MOD 24
  
  #Extract counts from histogram across 24 bins
  tmod24hist <- hist(tmod24,breaks=hist_breaks,plot=FALSE)
  histcounts <- tmod24hist$density
  
  #Find scaling of Fourier series that minimizes the distance between the transformed
  #Fourier series and the histogram counts with
  FS_transform <- optimize(cp_transf,c(0,10000),
                           ml_gam=ml_gam,ml_eta=ml_eta,
                           histcounts=histcounts,hist_mids=hist_mids)
  scal <- FS_transform$minimum
  
  #Calculate Fourier series to get smoother function for Fourier series
  fs_plot1 <- scal*scal_fs
  
  #Plot histogram of event times mod 24 with modified Fourier series laid over it
  hist(tmod24,prob=TRUE,
       breaks=hist_breaks, main=NULL, axes=FALSE, ylim=c(0,0.2),xlim=c(0,24))
  lines(x,fs_plot1,col="red")
  
  if(j %% 6 == 1){
    axis(2, at = seq(0, 0.2, by = 0.1), labels = seq(0, 0.2, by = 0.1))  # Custom Y-ticks
  }
  
  #if(j >= 13){
    axis(1, at = seq(0, 24, by = 12), labels = seq(0, 24, by = 12))  # Custom X-ticks
  #}
  box()
  return()
}

#HISTOGRAMS OF ALL INDIVIDUALS

#Load in directories for data and estimates
wdest <- "/home/xieryan/Dissertation1/Data/LSE_GH/"
wddata <- "/home/xieryan/Dissertation1/Data/Real_Data_Event_Times_V2_GitHub/"

#CLSE HISTOGRAMS
pdf(file=paste0(wdest,"lse_hists_v2_final.pdf"),width=6,height=6) #width=4,height=4
par(mfrow=c(7,6),mar = c(1,1,0.25,0.5), oma = c(3, 3, 2.5, 0.5), mgp = c(2, 0.2, 0), tck = -0.05) #mar = c(1,1,1,0.5)
for(j in 1:41){
  cum_hist(j,wdest,wddata)
}
mtext("Time of Day (hrs)", side = 1, outer = TRUE, line = 1.5, cex = 1)
mtext("Density", side = 2, outer = TRUE, line = 1.5, cex = 1)
dev.off()


