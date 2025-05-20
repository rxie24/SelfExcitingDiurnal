#FUNCTIONS TO PERFORM ROUGHNESS PENALTY

#Create function to plot fourier series in 24 hour interval
fourier <- function(x,gammas,eta){
  m <- length(x)
  M <- length(gammas)
  integers <- c(1:M)
  fouriers <- rep(NA,m)
  
  for(k in 1:m){
    #CALCULATE FOURIER SERIES
    xm <- x[k]
    f_terms <- gammas*sin(integers*xm*pi/12) + eta*cos(integers*xm*pi/12)
    fouriers[k] <- sum(f_terms)
  }
  return(fouriers)
}

#Create function to calculate the square of the second derivative of the Fourier series with respect to time
sqsec <- function(x,theta){
  
  #x: needs to be able to take in a vector of times (for integrand to work)
  #theta: (Eta,Gammas)
  
  ltheta <- length(theta)
  nfc <- ltheta/2
  
  eta <- theta[1:(1+nfc-1)]
  gammas <- theta[(1+nfc):ltheta]
  xl <- length(x)
  
  ints <- seq(1,nfc,1)
  
  sec_der <- rep(NA,xl)
  
  for(j in 1:xl){
    xj <- x[j]
    fseries <- (-1*gammas*((ints*pi)^2)/144)*sin(ints*xj*pi/12) + (-1*eta*((ints*pi)^2)/144)*cos(ints*xj*pi/12)
    sec_der[j] <- sum(fseries)^2
  }
  
  return(sec_der)
}

#Create function to calculate integral of the square of the second derivative of the Fourier series with respect to time to calculate roughness penalty
rpen <- function(theta){
  #theta: (Eta,Gammas)
  pen_int <- integrate(sqsec,lower=0,upper=24,theta=theta)
  
  return(pen_int[[1]])
}

#Get parameter estimates for an individual person
cp_transf = function(beta,ml_gam,ml_eta,histcounts,hist_mids) {
  
  #beta: scale and shift factors
  #mles: Fourier series parameters
  #histcounts: histogram counts of each of predetermined bins
  #hist_mids: midpoints of histogram bins
  
  #Calculate Fourier series and calculate least squares
  fs <- fourier(hist_mids,ml_gam,ml_eta)
  scaled_fs <- beta*(fs - min(fs))/(max(fs) - min(fs))
  
  #Calculate sum of squared error between optimized scaled/shifted Fourier series 
  #estimate and histogram counts
  sqerror <- sum((histcounts - scaled_fs)^2)
  
  #Return squared error
  return(sqerror)
}

cum_ss <- function(times,ml_fs){
  #id: individual
  #wdest: working directory of the estimates
  
  tmod24 <- times %% 24
  
  #Specific bin midpoints and breaks for histograms
  hist_mids <- seq(0.5,23.5,1)
  hist_breaks <- seq(0,24,1)
  
  #Extract eta and gammas from estimates
  lmles <- length(ml_fs)
  ncoef <- lmles/2
  ml_eta <- ml_fs[1:(1+ncoef-1)]
  ml_gam <- ml_fs[(1+ncoef):lmles]
  
  #Extract counts from histogram across 24 bins
  tmod24hist <- hist(tmod24,breaks=hist_breaks,plot=FALSE)
  histcounts <- tmod24hist$density
  
  #Find minimum SS of estimated Fourier series scaled to cumulative counts
  FS_transform <- optimize(cp_transf,c(0,10000),ml_gam=ml_gam,ml_eta=ml_eta,
                           histcounts=histcounts,hist_mids=hist_mids)
  ss_cum <- FS_transform$objective
  
  return(ss_cum)
}
