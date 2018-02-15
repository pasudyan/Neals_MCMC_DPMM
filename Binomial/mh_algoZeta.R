## Function to do Metropolis-Hastings algorithm     ##
## in the gibbs sampling for                        ##
## Binomial with Beta Prior                         ##
######################################################

######################################################
## - Prior: Beta(1,1)                               ##
## - dataMat contains Y_i, N_i, binPar, and cluster ##
## assignments                                      ##
## - Sampling for zeta                              ##
######################################################

mh_algoZeta <- function(dataMat){
  zetaCurr = dataMat$zeta[1]
  alpha_base = 1
  beta_base = 1
  mhloop = 5
  
  for(t in 1:mhloop){
    zetaCand = runif(1)

    #Calculating the acceptance probability 
    n_vec = sum(dbinom(dataMat$mixDat,
                       dataMat$trials, 
                       dataMat$theta*zetaCand, 
                       log = TRUE), 
                dbeta(zetaCand, 
                      alpha_base,
                      beta_base,
                      log = TRUE))
    
    d_vec = sum(dbinom(dataMat$mixDat,
                       dataMat$trials,
                       dataMat$theta*zetaCurr,
                       log = TRUE), 
                dbeta(zetaCurr,
                      alpha_base,
                      beta_base,
                      log = TRUE)) 
    
    ratio = n_vec-d_vec
    
    if (runif(1) < exp(ratio))
      params = zetaCand
    else 
      params = zetaCurr
    
  }
  return(params)
}
