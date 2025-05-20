
#Fit five models for data generated from some form of Circadian +/- Hawkes Processes
five_chp <- function(i,theta_chp,thetai_chp,thetai_lse,N){
  
  #i: iteration number (for running tasks in parallel)
  #theta_chp: parameters used to generate test and training data sets
  #thetai_chp: initial parameters used for Circadian Hawkes Process
  #thetai_lse: initial parameters used for Linear Self Exciting Processes
  #N: number of waiting times in testing and training data sets
  
  #Extract parameters to generate data
  ltheta <- length(theta_chp)
  mu <- theta_chp[1]
  beta <- theta_chp[2]
  sigma <- theta_chp[3]
  
  ncoef <- length(theta_chp[4:ltheta])/2
  eta <- theta_chp[4:(4+ncoef-1)]
  gammas <- theta_chp[(4+ncoef):ltheta]
  
  #Simulate training data and test data
  times_train <- gen_data_chp(N,mu,beta,sigma,eta,gammas)
  times_test <- gen_data_chp(N,mu,beta,sigma,eta,gammas)
  
  mlls <- five_models(times_train,times_test,thetai_chp,thetai_lse)
  return(mlls)
}

#Perform above function over multiple iterations, filter, and find median
five_chp_iter <- function(navg,theta_chp,thetai_chp,thetai_lse,N,ncores=60){
  mod_lls <- mclapply(1:navg,five_chp,theta_chp=theta_chp,thetai_chp=thetai_chp,
                      thetai_lse=thetai_lse,N=N, mc.cores = ncores)
  
  chp_lls <- rep(NA,navg)
  hp_lls <- rep(NA,navg)
  pp_lls <- rep(NA,navg)
  cp_lls <- rep(NA,navg)
  lse_lls <- rep(NA,navg)
  
  for(j in 1:navg){
    all_lls <- mod_lls[[j]]
    chp_lls[j] <- all_lls[1]
    hp_lls[j] <- all_lls[2]
    pp_lls[j] <- all_lls[3]
    cp_lls[j] <- all_lls[4]
    lse_lls[j] <- all_lls[5]
  }
  
  #Filter outliers
  chp_lls <- filt(chp_lls)
  hp_lls <- filt(hp_lls)
  pp_lls <- filt(pp_lls)
  cp_lls <- filt(cp_lls)
  lse_lls <- filt(lse_lls)
  
  #Track number of NA log-likelihoods there are and print them
  no_NAs <- c(sum(is.na(chp_lls)),sum(is.na(hp_lls)),sum(is.na(pp_lls)),sum(is.na(cp_lls)),
              sum(is.na(lse_lls)))
  print(no_NAs)
  
  #Calculate median of loglikelihoods
  chp_avg <- median(chp_lls,na.rm=TRUE)
  hp_avg <- median(hp_lls,na.rm=TRUE)
  pp_avg <- median(pp_lls,na.rm=TRUE)
  cp_avg <- median(cp_lls,na.rm=TRUE)
  lse_avg <- median(lse_lls,na.rm=TRUE)
  
  all_meds <- c(chp_avg,hp_avg,pp_avg,cp_avg,lse_avg)
  return(all_meds)
}
