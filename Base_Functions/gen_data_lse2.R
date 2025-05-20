
#Function to generate data from waiting times for Circadian/Hawkes Models
et_sim_lse <- function(wt,intstart,times,mu,beta,eta,gammas,ur,subdiv=100000){
  #ur: random variable generate from a Unif(0,1) distribution
  first <- log(ur)
  intend <- intstart + wt
  second <- second_lse_ll(intend,intstart,times,mu,beta,eta,gammas,subdiv)
  final <- first + second
  return(final)
}

#Function to generate data from waiting times for Circadian/Hawkes Models
gen_data_lse <- function(n,mu,beta,eta,gammas,subdiv=100000){
  #n: number of new points you want to generate
  times <- c(0,rep(NA,n)) #WLOG, set times[1] = 0
  times[2] <- 0 + rexp(1,rate=1)
  
  #Generate waiting times and points
  for(k in 2:n){
    ur <- runif(1,min=0,max=1)
    intstart <- times[k]
    tsubset <- times[1:k]
    et_optim <- uniroot(et_sim_lse,c(0,100),intstart=intstart,times=tsubset,
                        mu=mu,beta=beta,eta=eta,gammas=gammas,ur=ur,subdiv=subdiv)
    et_sim <- et_optim$root
    times[k+1] <- times[k] + et_sim
  }
  return(times)
}
