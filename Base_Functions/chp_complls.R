#Functions for calculating log-likelihood under Hawkes, Circadian, and Poisson models

#Function to calculate negative log-likelihood of Hawkes Model
hp_ll <- function(theta,times,subdiv=100000){
  ltheta <- length(theta)
  mu <- theta[1]
  beta <- theta[2]
  sigma <- theta[3]
  eta <- 0
  gammas <- 0
  
  lt <- length(times)
  lt2 <- lt - 1
  
  first_terms <- compute_lambda_hp(times[2:lt],times,mu,beta,sigma,eta,gammas)
  log_first <- log(first_terms)
  first <- sum(log_first)
  
  intstart <- times[1]
  intend <- times[lt]
  second <- second_chp_ll(intend,intstart,times,mu,beta,sigma,eta,gammas,subdiv)
  
  ll <- first - second
  return(-1*ll)
}

#Get MlEs under Hawkes Model
hp_mles <- function(theta,times,subdiv=100000){
  mhp <- optim(theta,hp_ll,times=times,subdiv=subdiv,method="Nelder-Mead")
  hp_mles <- mhp$par
  return(hp_mles)
}

#Function to calculate negative log-likelihood of Circadian Model
cp_ll <- function(theta,times,subdiv=100000){
  ltheta <- length(theta)
  mu <- theta[1]
  beta <- 0
  sigma <- 0
  
  ncoef <- length(theta[2:ltheta])/2
  eta <- theta[2:(2+ncoef-1)]
  gammas <- theta[(2+ncoef):ltheta]
  
  lt <- length(times)
  lt2 <- lt - 1
  
  first_terms <- compute_lambda_hp(times[2:lt],times,mu,beta,sigma,eta,gammas)
  log_first <- log(first_terms)
  first <- sum(log_first)
  
  intstart <- times[1]
  intend <- times[lt]
  second <- second_chp_ll(intend,intstart,times,mu,beta,sigma,eta,gammas,subdiv)
  
  ll <- first - second
  return(-1*ll)
}

#Get MlEs under Circadian Model
cp_mles <- function(theta,times,subdiv=100000){
  mcp <- optim(theta,cp_ll,times=times,subdiv=subdiv,method="Nelder-Mead")
  cp_mles <- mcp$par
  return(cp_mles)
}

#Poisson Process Log-likelihood
log_pp <- function(lambda,diffs){
  #lambda: value of parameter
  #diffs: waiting times
  
  lt2 <- length(diffs)
  first <- lt2*log(lambda)
  second <- lambda*sum(diffs)
  
  ll <- first  - second
  return(ll)
}

#Create function to find maximum likelihood under Poisson Process
max_logpp <- function(times){
  #diffs: waiting times
  
  #Find waiting times
  lt <- length(times)
  lt2 <- lt - 1
  diffs <- times[2:lt] - times[1:lt2]
  
  #Calculate MLE of exp(mu) = inverse of weighted average of waiting times
  inv_lam <- mean(diffs)
  lam_mle <- 1/inv_lam
  #Calculate maximum log-likelihood at MLE
  max_pp <- log_pp(lam_mle,diffs)
  
  return(max_pp)
}
