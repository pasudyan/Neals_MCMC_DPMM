#########################################################  
## MCMC methods for parameter estimation               ##
## Algorithm 3 of Neal's 2000                          ##
##                                                     ## 
## Observations: Y_i ~ Bin(N, theta_i)                 ##
## where theta_i has prior Beta(1,1)                   ## 
#########################################################

test_phi.3 = numeric(numIter)

#initialize number of clusters
K = floor(alpha*log(num_cust)) 

#initialize parameters for each cluster
phi.3 = numeric(num_cust) #for algo 3

#storing the cluster assignments for each iteration
cluster_store3 = numeric(numIter)

#initialize cluster assignment
sampling = c(1:K, sample(1:K, num_cust-K, replace=TRUE))
cluster  = sample(sampling, num_cust, replace=FALSE)

#calculate number of people in each table
counts = numeric(num_cust)
counts[1:K] = sapply(1:K, function(r) sum(cluster == r))

#conditional prob for each cluster
prob_clust = numeric(num_cust) 

#sum of observations in each cluster
sum_cluster = numeric(num_cust)

for (ii in 1:numIter){
  
  rand_samples = sample(1:num_cust,
                        num_cust,
                        replace=FALSE)
  
  #iteration for each observation
  for (b in 1:num_cust){
    
    j = rand_samples[b]
    
    #removing the observation from the cluster
    counts[cluster[j]] = counts[cluster[j]]-1
    
    #reevaluating cluster assignments
    if (counts[cluster[j]]==0){
      
      counts[cluster[j]]  = counts[K]
      cluster[cluster==K] = cluster[j]
      counts[K] = 0
      K = K-1
      
    }
    
    cluster[j] = 0
    
    #iteration for each table
    for (iter in 1:K){
      
      #prob of choosing existing table 
      a = sum(mixDat[cluster==iter]) + alpha_base 
      b = sum(n - mixDat[cluster==iter]) + beta_base
      
      prob_clust[iter] = log(counts[iter]) +
        log(choose(n,mixDat[j])) +
        log(beta(mixDat[j]+a, n-mixDat[j]+b)/beta(a,b))
    }
    
    #prob of choosing new table
    prob_clust[K+1] = log(alpha) + 
      log(choose(n, mixDat[j])) +
      log(beta(mixDat[j]+alpha_base, 
               n-mixDat[j]+beta_base)) -
      log(beta(alpha_base, beta_base))
    
    #normalizing constant
    prob_norm = prob_clust[1:(K+1)] -
      logSumExp(prob_clust[1:(K+1)])
    
    #sampling new cluster assignments
    new_table = sample(1:(K+1),
                       1,
                       replace=FALSE,
                       exp(prob_norm))
    
    #new table addition 
    cluster[j] <- new_table
    if (new_table == (K+1)){
      counts[K+1] = 1
      K = K+1
    } 
    else {
      counts[new_table] = counts[new_table]+1
    }
  }
  
  #sampling the parameter for each cluster
  sum_cluster = sapply(1:K, function(r) 
    sum(mixDat[cluster==r]))
  
  a = sum_cluster + alpha_base 
  
  b = counts[1:K]*n - sum_cluster + beta_base
  
  phi.3[1:K] = rbeta(K,a,b)
  
  #storing params value for observation 
  test_phi.3[ii] = phi.3[cluster[test_obs]]
  cluster_store3[ii] = cluster[test_obs]
  
#   if ( ii%%telliter == 0)
#     print(sprintf("Iteration: %d", ii))
  
}



