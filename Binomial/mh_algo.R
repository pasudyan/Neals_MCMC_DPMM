#Metropolis-Hastings algorithm for Binomial-Beta
mh_algo <- function(dataMat, phi.8MH){
  nuCurr = phi.8MH[1:K]
  alpha_base = 1
  eta_base = 1
  
  ratio = numeric(length(nuCurr))
  params = numeric(length(nuCurr))
  
  for(t in 1:5){
    nuCand = rbeta(length(nuCurr),1,1)
    
    #Calculating the acceptance probability 
    for (i in 1:length(nuCurr)){
      mydatta  = dataMat[dataMat$cluster==i,]
      
      n_vec = sum(dbinom(mydatta$mixDat,
                         mydatta$trials, 
                         mydatta$zeta*nuCand[i], 
                         log = TRUE), 
                  dbeta(nuCand[i], 
                        alpha_base,
                        beta_base,
                        log = TRUE))
      
      d_vec = sum(dbinom(mydatta$mixDat,
                         mydatta$trials,
                         mydatta$zeta*nuCurr[i],
                         log = TRUE), 
                  dbeta(nuCurr[i],
                        alpha_base,
                        beta_base,
                        log = TRUE)) 
      
      
      ratio = n_vec-d_vec
      
      if (runif(1) < exp(ratio))
        params[i] = nuCand[i]
      else  
        params[i]= nuCurr[i]
    }
    nuCurr = params
  }
  return(nuCurr)
}