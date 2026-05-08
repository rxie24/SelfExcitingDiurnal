library(Rcpp)
library(parallel)

#Function to compute lambda to calculate lambda across multiple time points
cppFunction("
std::vector<double> compute_lambda_hp(const std::vector<double>& t, const std::vector<double>& times, double mu, double sigma, const std::vector<double>& eta, const std::vector<double>& gammas) {
    int m = t.size();
    std::vector<double> lam(m);

    for (int k = 0; k < m; ++k) {
        lam[k] = lambda(t[k], times, mu, sigma, eta, gammas);
    }

    return lam;
}",
            includes = {"double lambda(double t, const std::vector<double>& times, double mu, double sigma, const std::vector<double>& eta, const std::vector<double>& gammas,int ne = 5) {

    int M = gammas.size();
    double result = mu;

    // Calculate the Fourier series component
    double fourier = 0.0;
    double pi_over_12 = M_PI / 12.0;
    for (int i = 0; i < M; ++i) {
        double angle = (i + 1) * t * pi_over_12;
        fourier += (gammas[i] * std::sin(angle)) + (eta[i] * std::cos(angle));
    }
    // result += fourier;

    // Calculate the integral component
    double sum_int = 0.0;
    int nt2 = std::count_if(times.begin(), times.end(), [t](double val) { return val < t; });
    int n2 = std::min(nt2,ne);
    for (int i = 0; i < n2; ++i) {
        double quant = t - times[nt2 - i - 1];
        double norm_cdf = 0.5 * (1.0 + erf((quant - 0.0) / (exp(sigma) * sqrt(2.0))));
        double norm_cdf2 = 1.0 - norm_cdf;
        sum_int += norm_cdf2;
          
        if (norm_cdf2 < 1e-08) {
            break;
        }
      }
  
// Calculate lambda(t)
    result += (fourier * sum_int);
    return std::exp(result);
};
  "})

#Function to calculate integrate of lambda from some start time to some end time
second_chp_ll <- function(intend,intstart,times,mu,sigma,eta,gammas,subdiv=100000){
  second <- integrate(compute_lambda_hp,intstart,intend,times=times,
                      mu=mu,sigma=sigma,eta=eta,gammas=gammas,subdivisions=subdiv,
                      stop.on.error=FALSE)
  second <- second$value
  return(second)
}

#Function to calculate negative log-likelihood of CHP Model
chp_ll <- function(theta,times,subdiv=100000){
  ltheta <- length(theta)
  mu <- theta[1]
  sigma <- theta[2]
  
  ncoef <- length(theta[3:ltheta])/2
  
  eta <- theta[3:(3+ncoef-1)]
  gammas <- theta[(3+ncoef):ltheta]
  
  lt <- length(times)
  lt2 <- lt - 1
  
  first_terms <- compute_lambda_hp(times[2:lt],times,mu,sigma,eta,gammas)
  log_first <- log(first_terms)
  first <- sum(log_first)
  
  intstart <- times[1]
  intend <- times[lt]
  second <- second_chp_ll(intend,intstart,times,mu,sigma,eta,gammas,subdiv)
  
  ll <- first - second
  return(-1*ll)
}

#Get MlEs under Circadian Hawkes Model
chp_mles <- function(theta,times,subdiv=100000){
  mchp <- optim(theta,chp_ll,times=times,subdiv=subdiv,method="Nelder-Mead")
  chp_mles <- mchp$par
  return(chp_mles)
}

chp_llest <- function(k,wd,wddata,nfc){
  
  #Load in the data
  times <- readRDS(paste0(wddata,"data",k,".rds"))
  
  #Fit Hawkes Process Modulated by Fourier Series
  nfc_k <- nfc[k]
  thetai_chp <- rep(0,(2+nfc_k))
  mles_chp <- chp_mles(thetai_chp,times)
  #print(mles_chp) #Print estimate
  
  #Save mles of full Circadian Hawkes model
  saveRDS(object = mles_chp,file=paste0(wd,"ests",k,".rds"))  
  ll_chp <- -1*chp_ll(mles_chp,times)
  
  #Return log-likelihood
  return(ll_chp)
}

#Establish working directories to save parameter estimates and where data is located 
wd <- "/home/xieryan/Dissertation1/Data/Mod_CHP2/"
wddata <- "/home/xieryan/Dissertation1/Data/Real_Data_Event_Times_V2_GitHub/"

#Subset individuals to perform estimation
ids <- seq(1,41,1)
ncores <- length(ids)
nfc <- readRDS("/home/xieryan/Dissertation1/Data/CHP/chp_lengths.rds")

#Fit Mod CHP Models for individuals
all_lls <- mclapply(ids,chp_llest,wd=wd,wddata=wddata,nfc=nfc,mc.cores=ncores)
print(all_lls)

id_lls <- unlist(all_lls)
#Save vector of log-likelihoods when fitting four models
saveRDS(object = id_lls,file=paste0(wd,"Mod_CHP_LLs0.rds"))
