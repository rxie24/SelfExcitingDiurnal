
#Create function to run chp_five function over multiple values of parameter value to get list of vectors
#of test-loglikelihoods across all five models

test_lls <- function(param_seq,navg,theta_wb,mp,thetai_chp,thetai_lse,N,ncores=60){
  pl <- length(param_seq)
  mind <- mp-1
  chp_testlls <- rep(NA,pl)
  hp_testlls <- rep(NA,pl)
  pp_testlls <- rep(NA,pl)
  cp_testlls <- rep(NA,pl)
  lse_testlls <- rep(NA,pl)
  
  for(k in 1:pl){
    param_val <- param_seq[k]
    theta_chp <- append(theta_wb,param_val,after=mind)

    ll_meds <- five_chp_iter(navg,theta_chp,thetai_chp,thetai_lse,N,ncores=ncores)
    
    chp_testlls[k] <- ll_meds[1]
    hp_testlls[k] <- ll_meds[2]
    pp_testlls[k] <- ll_meds[3]
    cp_testlls[k] <- ll_meds[4]
    lse_testlls[k] <- ll_meds[5]
  }
  
  five_lls <- list(chp_testlls,hp_testlls,pp_testlls,cp_testlls,lse_testlls)
  return(five_lls)
}

