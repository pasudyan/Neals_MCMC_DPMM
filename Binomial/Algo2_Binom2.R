#########################################################  
## MCMC methods for parameter estimation               ##
## Algorithm 2 of Neal's 2000                          ##
##                                                     ## 
## Observations: Y_i ~ Bin(N_i, theta_i)               ##
## where theta_i has prior Beta(1,1)                   ## 
#########################################################

test_phi.2 = numeric(numIter)

#initialize number of clusters
K = floor(alpha*log(num_cust)) 

#initialize parameters for each cluster
phi.2 = numeric(num_cust)
phi.2[1:K] = rbeta(K, alpha_base, beta_base) 

#initialize cluster assignment
sampling = c(1:K, sample(1:K, num_cust-K, replace=TRUE))
cluster = sample(sampling, num_cust, replace=FALSE) 

#calculate number of people in each table
counts = numeric(num_cust)
counts[1:K] = sapply(1:K, function(r) sum(cluster == r))

##Conditional prob for each cluster
prob_clust = numeric(num_cust)

#sum of observations in each cluster
sum_cluster = numeric(num_cust)
sum_trials = numeric(num_cust)

for (ii in 1:numIter){
  
  rand_samples = sample(1:num_cust, 
                        num_cust, 
                        replace=FALSE)
  
  #iteration for each observation
  for (b in 1:num_cust){
    
    j = rand_samples[b]
    
    #removing the observation from the cluster
    counts[cluster[j]] <- counts[cluster[j]]-1
    
    #reevaluating cluster assignments
    if (counts[cluster[j]]==0){
      
      #renumbering clusters
      counts[cluster[j]] = counts[K]
      cluster[cluster==K] = cluster[j]
      counts[K] = 0
      
      #remove and renumber the parameter phi.2
      phi.2[cluster[j]] = phi.2[K]
      phi.2[K] = 0
      K = K-1
    } 
    
    #prob of choosing existing table 
    prob_clust[1:K] = log(counts[1:K]) + 
      dbinom(mixDat[j],
             trials[j],
             phi.2[1:K],
             log = TRUE)
    
    #prob of choosing new table
    prob_clust[K+1] = log(alpha) + 
      log(choose(trials[j],mixDat[j])) +
      log(beta(
        mixDat[j]+alpha_base,
        trials[j]-mixDat[j]+beta_base)) -
      log(beta(alpha_base, beta_base))
    
    #normalizing constant
    prob_norm = prob_clust[1:(K+1)] - 
      logSumExp(prob_clust[1:(K+1)])
    
    #sampling new cluster assignments for j'th obs
    new_table = sample(1:(K+1),
                       1,
                       replace=FALSE,
                       exp(prob_norm))
    
    #new table addition 
    cluster[j] = new_table
    if (new_table == (K+1))
    {
      counts[K+1] = 1
      
      #sample from the posterior
      a = mixDat[j] + alpha_base
      b = trials[j] - mixDat[j] + beta_base
      phi.2[K+1] = rbeta(1, a, b)
      K = K+1
    } 
    else {
      counts[new_table] = counts[new_table]+1
    }
  }
  
  #sampling the parameter for the ii'th iteration
  sum_cluster = sapply(1:K, function(r) 
    sum(mixDat[cluster==r]))
  
  a = sum_cluster + alpha_base 
  
  sum_trials = sapply(1:K, function(r) 
    sum(trials[cluster==r]))
  
  b = sum_trials - sum_cluster + beta_base
  
  phi.2[1:K] = rbeta(K,a,b)
  
  #storing params value for observation 
  test_phi.2[ii] = phi.2[cluster[test_obs]]

#     if ( ii%%telliter == 0)
#       print(sprintf("Iteration: %d", ii))
  
} 






















