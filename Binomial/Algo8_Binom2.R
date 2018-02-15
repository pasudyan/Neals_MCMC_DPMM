#########################################################  
## MCMC methods for parameter estimation               ##
## Algorithm 8 of Neal's 2000                          ##
##                                                     ## 
## Observations: Y_i ~ Bin(N_i, theta_i)               ##
## where theta_i has prior Beta(1,1)                   ## 
#########################################################

test_phi.8 = numeric(numIter)

#initialize number of clusters
K = floor(alpha*log(num_cust))

#initialize cluster assignment
phi.8 = numeric(xtraSpace)
phi.8[1:K] = rbeta(K, alpha_base, beta_base) #for algo 2

#initialize cluster assignment
sampling = c(1:K, sample(1:K, num_cust-K, replace=TRUE))
cluster  = sample(sampling, num_cust, replace=FALSE)

#calculate number of people in each table
counts = numeric(xtraSpace)
counts[1:K] = sapply(1:K, function(r) sum(cluster == r))

#conditional prob for each cluster
prob_clust = numeric(xtraSpace) 

#sum of observations in each cluster
sum_cluster = numeric(xtraSpace)
sum_trials = numeric(xtraSpace)

for (ii in 1:numIter){
  
  rand_samples = sample(1:num_cust,
                        num_cust,
                        replace=FALSE)
  
  #iteration for each observation
  for (b in 1:num_cust){
    
    j = rand_samples[b]
    
    h = K + m
    
    #removing the observation from the cluster
    counts[cluster[j]] = counts[cluster[j]]-1
    
    #reevaluating cluster assignments
    if (counts[cluster[j]]!=0){
      
      #drawing values for the aux params from G_0
      phi.8[(K+1):h] = rbeta(m, alpha_base, beta_base) 
      
    } else if (counts[cluster[j]]==0 & cluster[j]!=K){
      
      #renumbering the clusters
      counts[cluster[j]] = counts[K]
      cluster[cluster==K] = cluster[j]
      counts[K] = 0
      
      #renumbering the phi.8
      temp = phi.8[cluster[j]]
      phi.8[cluster[j]] = phi.8[K]
      phi.8[K] = temp
      
      #reducing the number of clusters and aux clusters
      K = K-1
      h = h-1
      
      #drawing values for the aux params from G_0
      phi.8[(K+2):h] = rbeta(m-1, alpha_base, beta_base)
      
    } else if (counts[cluster[j]]==0 & cluster[j]==K) {
      
      #reducing the number of clusters and aux clusters
      K = K-1
      h = h-1
      
      #drawing values for the aux params from G_0
      phi.8[(K+2):h] = rbeta((m-1), alpha_base, beta_base)
      
    }
    
    #prob of choosing existing cluster 
    prob_clust[1:K] = log(counts[1:K]) +
      dbinom(mixDat[j], 
             trials[j],
             phi.8[1:K],
             log=TRUE)
    
    #prob of choosing new cluster
    prob_clust[(K+1):h] = log(alpha/m) + 
      dbinom(mixDat[j],
             trials[j],
             phi.8[(K+1):h], 
             log=TRUE)
    
    #normalizing constant
    prob_norm = prob_clust[1:h] -
      logSumExp(prob_clust[1:h])
    
    #sampling new cluster assignments
    new_table = sample(1:h,
                       1,
                       replace=FALSE,
                       exp(prob_norm))
    
    #new table addition
    if (new_table > K){
      
      cluster[j] = K+1
      phi.8[K+1] = phi.8[new_table]
      counts[K+1] = 1
      K = K+1
      
    } else {
      
      cluster[j] = new_table
      counts[new_table] = counts[new_table]+1
      
    }
    
    phi.8[(K+1):xtraSpace]  = rep(0,xtraSpace-K)
    counts[(K+1):xtraSpace] = rep(0,xtraSpace-K)
    
  } 
  
  #sampling the new parameters
  #this one here could be done using the MH algorithm
  sum_cluster = sapply(1:K, function(r) 
    sum(mixDat[cluster==r]))
  
  a = sum_cluster + alpha_base 
  
  sum_trials = sapply(1:K, function(r) 
    sum(trials[cluster==r]))
  
  b = sum_trials - sum_cluster + beta_base
  
  phi.8[1:K] = rbeta(K,a,b)
  
  #storing params value for observation 
  test_phi.8[ii] = phi.8[cluster[test_obs]]

}
