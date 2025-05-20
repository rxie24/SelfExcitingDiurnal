#Functions for calculating and optimizing on Log-likelihoods for Poisson, LSE (SE Only), and Circadian 
#(looking only at subset of time that full LSE model uses to calculate log-likelihood)

#Function to calculate negative log-likelihood of LSE Model w/o Circadian component
lse_nc_ll <- function(theta,times,b=2,subdiv=100000){
  #b: order of beta (number of points to go back to)
  e_ind <- 2 + b
  ltheta <- length(theta)
  
  mu <- theta[1]
  beta <- theta[2:(e_ind-1)]
  eta <- 0
  gammas <- 0
  
  lt <- length(times)
  start <- b+1
  
  tsubset <- times[start:lt]
  ltsub <- length(tsubset)
  ltsub2 <- ltsub - 1
  
  first_terms <- compute_lambda_lse(tsubset[2:ltsub],times,mu,beta,eta,gammas)
  log_first <- log(first_terms)
  first <- sum(log_first)
  
  intstart <- times[start]
  intend <- times[lt]
  second <- second_lse_ll(intend,intstart,times,mu,beta,eta,gammas,subdiv=subdiv)
  
  ll <- first - second
  return(-1*ll)
}

#Get MLEs under full LSE Model
lse_nc_mles <- function(theta,times,b=2,subdiv=100000){
  mlse_nc <- optim(theta,lse_nc_ll,times=times,b=b,subdiv=subdiv,method="Nelder-Mead")
  lse_nc_mles <- mlse_nc$par
  return(lse_nc_mles)
}

#Function to calculate negative log-likelihood of full Circadian only Model
cp_ll2 <- function(theta,times,b=2,subdiv=100000){
  #b: order of beta (number of points to go back to)
  e_ind <- 2 + b
  ltheta <- length(theta)
  
  mu <- theta[1]
  beta <- c(0,0)
  
  ncoef <- length(theta[2:ltheta])/2
  eta <- theta[2:(2+ncoef-1)]
  gammas <- theta[(2+ncoef):ltheta]
  
  lt <- length(times)
  start <- b+1
  
  tsubset <- times[start:lt]
  ltsub <- length(tsubset)
  ltsub2 <- ltsub - 1
  
  first_terms <- compute_lambda_lse(tsubset[2:ltsub],times,mu,beta,eta,gammas)
  log_first <- log(first_terms)
  first <- sum(log_first)
  
  intstart <- times[start]
  intend <- times[lt]
  second <- second_lse_ll(intend,intstart,times,mu,beta,eta,gammas,subdiv=subdiv)
  
  ll <- first - second
  return(-1*ll)
}

#Get MLEs under Circadian Only Model
cp_mles2 <- function(theta,times,b=2,subdiv=100000){
  mcp2 <- optim(theta,cp_ll2,times=times,b=b,subdiv=subdiv,method="Nelder-Mead")
  cp2_mles <- mcp2$par
  return(cp2_mles)
}

#Poisson Process Log-likelihood
log_pp2 <- function(lambda,diffs){
  #lambda: value of parameter
  #diffs: waiting times
  
  lt2 <- length(diffs)
  first <- lt2*log(lambda)
  second <- lambda*sum(diffs)
  
  ll <- first  - second
  return(ll)
}

#Create function to find maximum likelihood under Poisson Process
max_logpp2 <- function(times,b=2){
  #diffs: waiting times
  
  #Find waiting times for respective subset of time
  lt <- length(times)
  start <- b+1
  
  tsubset <- times[start:lt]
  ltsub <- length(tsubset)
  ltsub2 <- ltsub - 1
  
  diffs <- tsubset[2:ltsub] - tsubset[1:ltsub2]
  
  #Calculate MLE of exp(mu) = inverse of average of waiting times
  inv_lam <- mean(diffs)
  lam_mle <- 1/inv_lam
  #Calculate maximum log-likelihood at MLE
  max_pp <- log_pp2(lam_mle,diffs)
  
  return(max_pp)
}

