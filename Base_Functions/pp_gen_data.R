
#Function to generate data from waiting times from Poisson processes
gen_data_pp <- function(n,mu){
  #n: number of new points you want to generate
  #mu: log(lambda) parameter
  
  times <- c(0,rep(NA,n)) #WLOG, set times[1] = 0
  wts <- rexp(n,rate=exp(mu))
  
  #Generate waiting times and points
  for(k in 1:n){
    times[k+1] <- times[k] + wts[k]
  }
  return(times)
}
