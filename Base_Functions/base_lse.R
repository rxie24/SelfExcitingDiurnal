#Base Functions for Linear Self-Exciting Log-Likelihoods and Full LSE Model Optimization

#Function to calculate lambda at multiple time points
cppFunction("
std::vector<double> compute_lambda_lse(const std::vector<double>& t, const std::vector<double>& times, double mu, const std::vector<double>& beta, const std::vector<double>& eta, const std::vector<double>& gammas) {
    int m = t.size();
    std::vector<double> lam(m);

    for (int k = 0; k < m; ++k) {
        lam[k] = lambda(t[k], times, mu, beta, eta, gammas);
    }

    return lam;
}",
            includes = {"double lambda(double t, const std::vector<double>& times, double mu, const std::vector<double>& beta, const std::vector<double>& eta, const std::vector<double>& gammas) {
    int M = gammas.size();
    int p = beta.size();
    
    double result = mu;
    
    // Calculate the Fourier series component
    double fourier = 0.0;
    double pi_over_12 = M_PI / 12.0;
    for (int i = 0; i < M; ++i) {
        double angle = (i + 1) * t * pi_over_12;
        fourier += (gammas[i] * std::sin(angle)) + (eta[i] * std::cos(angle));
    }
    result += fourier;

    // Calculate the self-exciting component
    double se_sum = 0.0;
    int nt2 = std::count_if(times.begin(), times.end(), [t](double val) { return val < t; });
    se_sum += beta[0] * (t - times[nt2 - 1]);
    
    for (int i = 1; i < p; ++i) {
        double diff = times[nt2-i] - times[nt2-i-1];
        se_sum += diff*beta[i];
    } 

    result += se_sum;

    // Calculate lambda(t)
    return std::exp(result);
};
  "})

#Create function to integrate lambda from starting time point to ending time point
second_lse_ll <- function(intend,intstart,times,mu,beta,eta,gammas,subdiv=100000){
  second <- integrate(compute_lambda_lse,intstart,intend,times=times,
                      mu=mu,beta=beta,eta=eta,gammas=gammas,subdivisions=subdiv,
                      stop.on.error=FALSE)
  second <- second$value
  return(second)
}

#Function to calculate negative log-likelihood of full LSE Model
lse_ll <- function(theta,times,b=2,subdiv=100000){
  #b: order of beta (number of points to go back to)
  e_ind <- 2 + b
  ltheta <- length(theta)
  
  mu <- theta[1]
  beta <- theta[2:(e_ind-1)]
  
  ncoef <- length(theta[e_ind:ltheta])/2
  eta <- theta[e_ind:(e_ind+ncoef-1)]
  gammas <- theta[(e_ind+ncoef):ltheta]
  
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
lse_mles <- function(theta,times,b=2,subdiv=100000){
  mlse <- optim(theta,lse_ll,times=times,b=b,subdiv=subdiv,method="Nelder-Mead")
  lse_mles <- mlse$par
  return(lse_mles)
}

