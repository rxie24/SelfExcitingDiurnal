
#Fit five models to some given set of training data and get test data log-likelihoods
five_models <- function(times_train,times_test,thetai_chp,thetai_lse){
  #times_train: training data set
  #times_test: test data set
  #thetai_chp: initial parameters used for Circadian Hawkes Process
  #thetai_lse: initial parameters used for Linear Self Exciting Processes
  
  #Set initial parameter estimates for all models
  ltheta <- length(thetai_chp)
  mu1 <- thetai_chp[1]
  beta1 <- thetai_chp[2]
  sigma1 <- thetai_chp[3]
  
  ncoef <- length(thetai_chp[4:ltheta])/2
  eta1 <- thetai_chp[4:(4+ncoef-1)]
  gammas1 <- thetai_chp[(4+ncoef):ltheta]

  thetai_hp <- c(mu1,beta1,sigma1)
  thetai_cp <- c(mu1,eta1,gammas1)
  
  #Calculate average waiting time of training dataset
  lt <- length(times_train)
  wts <- times_train[2:lt] - times_train[1:(lt-1)]
  mtrain <- mean(wts)
  
  #Extract number of waiting times
  n <- length(wts)
  
  #Calculate waiting times of testing dataset
  lt2 <- length(times_test)
  wts2 <- times_test[2:lt2] - times_test[1:(lt2-1)]
  
  #Perform parameter estimations
  
  #Circadian Hawkes
  tryCatch({
    mles_chp <- chp_mles(thetai_chp,times_train)
    chp_testll <- -1*chp_ll(mles_chp,times_test)
  },
  error=function(e){
    chp_testll <<- NaN
  })
  
  #Hawkes
  tryCatch({
    mles_hp <- hp_mles(thetai_hp,times_train)
    hp_testll <- -1*hp_ll(mles_hp,times_test)
  },
  error=function(e){
    hp_testll <<- NaN
  })
  
  #Poisson
  #Get MLE of Lambda as 1/average waiting time
  tryCatch({
    lam <- 1/mtrain
    pp_testll <- log_pp(lam,wts2)
  },
  error=function(e){
    pp_testll <<- NaN
  })

  
  #Circadian
  tryCatch({
    mles_cp <- cp_mles(thetai_cp,times_train)
    cp_testll <- -1*cp_ll(mles_cp,times_test)
  },
  error=function(e){
    cp_testll <<- NaN
  })
  
  #Linear Self-Exciting Model
  tryCatch({
    mles_lse <- lse_mles(thetai_lse,times_train,b=2)
    lse_testll <- -1*lse_ll(mles_lse,times_test,b=2)
  },
  error=function(e){
    lse_testll <<- NaN
  })
  
  #Save test LLs into vector
  results <- rep(NA,5)
  results[1] <- chp_testll
  results[2] <- hp_testll
  results[3] <- pp_testll
  results[4] <- cp_testll
  results[5] <- lse_testll
  
  #Scale log-likelihoods by dividing by number of waiting times(n)
  results <- results/n
  return(results)
}
