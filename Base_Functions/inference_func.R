
#Create function to generate a dataset given parameter value, then calculate estimate (CHP Model)
chp_par_bootstrap_i <- function(i,theta,n,subdiv=100000){
  
  #Extract corresponding parameter values
  ltheta <- length(theta)
  mu <- theta[1]
  beta <- theta[2]
  sigma <- theta[3]
  
  ncoef <- length(theta[4:ltheta])/2
  eta <- theta[4:(4+ncoef-1)]
  gammas <- theta[(4+ncoef):ltheta]
  
  #Simulate data based on parameter values
  times <- gen_data_chp(n,mu,beta,sigma,eta,gammas)
  
  #Get parameter estimates from new simulated data
  thetai <- rep(0,ltheta)
  
  tryCatch({
    mchp <- optim(thetai,chp_ll,times=times,subdiv=subdiv,method="Nelder-Mead")
    chp_mles <- mchp$par
  },
  error=function(e){
    chp_mles <<- rep(NaN,ltheta)
  })
  
  return(chp_mles)
}

#Create function to perform parametric bootstrapping given parameter estimates from dataset (CHP Model)
#and extract 95% confidence interval
chp_par_bootstrap <- function(theta,n,ncores=10,niter=100,subdiv=100000){
  
  #Generate niter (default=100) datasets and extract parameter estimates from them
  bs_samples <- mclapply(1:niter,chp_par_bootstrap_i,theta=theta,n=n,subdiv=subdiv,mc.cores=ncores)
  
  #Create matrix of all parameter estimates
  bs_ests <- matrix(NA,nrow=niter,ncol=length(theta))
  for(j in 1:niter){
    bs_ests[j,] <- bs_samples[[j]]
  }
  
  #Extract 2.5th and 97.5th percentiles of all estimates to get confidence intervals
  mle_ci <- apply(bs_ests,2,quantile,probs=c(0.025,0.975),na.rm=TRUE)
  return(mle_ci)
}

#Create function to generate a dataset given parameter value, then calculate estimate (CLSE Model)
lse_par_bootstrap_i <- function(i,theta,n,b=2,subdiv=100000){
  
  #Extract corresponding parameter values
  e_ind <- 2 + b
  ltheta <- length(theta)
  
  mu <- theta[1]
  beta <- theta[2:(e_ind-1)]
  
  ncoef <- length(theta[e_ind:ltheta])/2
  eta <- theta[e_ind:(e_ind+ncoef-1)]
  gammas <- theta[(e_ind+ncoef):ltheta]
  
  #Simulate data based on parameter values
  times <- gen_data_lse(n,mu,beta,eta,gammas)
  
  #Get parameter estimates from new simulated data
  thetai <- theta
  
  tryCatch({
    mlse <- optim(thetai,lse_ll,times=times,b=b,subdiv=subdiv,method="Nelder-Mead")
    lse_mles <- mlse$par
  },
  error=function(e){
    lse_mles <<- rep(NaN,ltheta)
  })
  
  return(lse_mles)
}


#Create function to perform parametric bootstrapping given parameter estimates from dataset (CLSE Model)
#and extract 95% confidence interval
lse_par_bootstrap <- function(theta,n,b=2,ncores=10,niter=100,subdiv=100000){
  
  #Generate niter (default=100) datasets and extract parameter estimates from them
  bs_samples <- mclapply(1:niter,lse_par_bootstrap_i,theta=theta,n=n,b=b,subdiv=subdiv,mc.cores=ncores)
  
  #Create matrix of all parameter estimates
  bs_ests <- matrix(NA,nrow=niter,ncol=length(theta))
  for(j in 1:niter){
    bs_ests[j,] <- bs_samples[[j]]
  }
  
  #Extract 2.5th and 97.5th percentiles of all estimates to get confidence intervals
  mle_ci <- apply(bs_ests,2,quantile,probs=c(0.025,0.975),na.rm=TRUE)
  return(mle_ci)
}

