##MCMC methods for parameter estimation (Algorithm 8)##
#######################################################

test_phi.8MH3 = numeric(numIter)

#initialize number of clusters
K = floor(alpha*log(num_cust))

#initialize cluster assignment
phi.8MH3      = numeric(xtraSpace)
phi.8MH3[1:K] = rbeta(K, alpha_base, beta_base) 

#initialize cluster assignment
sampling = c(1:K, sample(1:K, num_cust-K, replace=TRUE))
cluster  = sample(sampling, num_cust, replace=FALSE)

dataMat$cluster = cluster

#calculate number of people in each table
counts      = numeric(xtraSpace)
counts[1:K] = sapply(1:K, function(r) sum(dataMat$cluster == r))

#conditional prob for each cluster
prob_clust = numeric(xtraSpace) 

#sum of observations in each cluster
sum_cluster = numeric(xtraSpace)
sum_trials  = numeric(xtraSpace)

for (ii in 1:numIter){
  
  rand_samples = sample(1:num_cust,
                        num_cust,
                        replace=FALSE)
  
  #iteration for each observation
  for (b in 1:num_cust){
    
    j = rand_samples[b]
    
    h = K + m
    
    #removing the observation from the cluster
    counts[dataMat$cluster[j]] = counts[dataMat$cluster[j]]-1
    
    #reevaluating cluster assignments
    if (counts[dataMat$cluster[j]]!=0){
      
      #drawing values for the aux params from G_0
      phi.8MH3[(K+1):h] = rbeta(m, alpha_base, beta_base) 
      
    } else if (counts[dataMat$cluster[j]]==0 & dataMat$cluster[j]!=K){
      
      #renumbering the clusters
      counts[dataMat$cluster[j]] = counts[K]
      dataMat$cluster[dataMat$cluster == K] = dataMat$cluster[j]
      counts[K] = 0
      
      #renumbering the phi.8MH3
      temp = phi.8MH3[dataMat$cluster[j]]
      phi.8MH3[dataMat$cluster[j]] = phi.8MH3[K]
      phi.8MH3[K] = temp
      
      #reducing the number of clusters and aux clusters
      K = K-1
      h = h-1
      
      #drawing values for the aux params from G_0
      phi.8MH3[(K+2):h] = rbeta(m-1, alpha_base, beta_base)
      
    } else if (counts[dataMat$cluster[j]]==0 & dataMat$cluster[j]==K) {
      
      #reducing the number of clusters and aux clusters
      K = K-1
      h = h-1
      
      #drawing values for the aux params from G_0
      phi.8MH3[(K+2):h] = rbeta((m-1), alpha_base, beta_base)
      
    }
    
    #prob of choosing existing cluster 
    prob_clust[1:K] = log(counts[1:K]) +
      dbinom(dataMat$mixDat[j], 
             dataMat$trials[j],
             dataMat$zeta[j]*phi.8MH3[1:K],
             log=TRUE)
    
    #prob of choosing new cluster
    prob_clust[(K+1):h] = log(alpha/m) + 
      dbinom(dataMat$mixDat[j],
             dataMat$trials[j],
             dataMat$zeta[j]*phi.8MH3[(K+1):h], 
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
      
      dataMat$cluster[j] = K+1
      phi.8MH3[K+1]      = phi.8MH3[new_table]
      counts[K+1] = 1
      K = K+1
      
    } else {
      
      dataMat$cluster[j] = new_table
      counts[new_table]  = counts[new_table]+1
      
    }
    
    phi.8MH3[(K+1):xtraSpace]  = rep(0,xtraSpace-K)
    counts[(K+1):xtraSpace]    = rep(0,xtraSpace-K)
    
  } 
  
  #sampling the new parameters using the MH Algorithm 
  phi.8MH3[1:K] = mh_algo(dataMat, 
                          phi.8MH3[1:K])
  
  #storing params value for observation 
  test_phi.8MH3[ii] = phi.8MH3[dataMat$cluster[test_obs]]
  
  if (ii%%telliter == 0)
    print(sprintf("Iteration: %d", ii))
  
}
