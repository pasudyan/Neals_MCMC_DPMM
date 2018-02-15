##MCMC methods for parameter estimation (Algorithm 2)##
#######################################################

test_phi.2 = numeric(numIter)

#initialize number of clusters
K = floor(alpha*log(num_cust)) 

#initialize parameters for each cluster
phi.2 = numeric(xtraSpace)
phi.2[1:K] = rnorm(K, mu_base, sd_base) #for algo 2

#storing the cluster assignments for each iteration
cluster_store2 = numeric(numIter)

#initialize cluster assignment
sampling = c(1:K, sample(1:K, num_cust-K, replace=TRUE))
cluster  = sample(sampling, num_cust, replace=FALSE) 

#calculate number of people in each table
counts = numeric(xtraSpace)
counts[1:K] = sapply(1:K, function(r) sum(cluster == r))

#Conditional prob for each cluster
prob_clust <- numeric(xtraSpace) 

#sum of observations in each cluster
sum_cluster = numeric(xtraSpace)

for (ii in 1:numIter){
  
  rand_samples = sample(1:num_cust, 
                        num_cust,
                        replace=FALSE
                        )
  
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
      dnorm(mixDat[j],
            phi.2[1:K],
            sd_obs,
            log = TRUE)
    
    #prob of choosing new table
    prob_clust[K+1] = log(alpha) + 
    dnorm(mixDat[j],
         mu_base,
         sqrt(var_obs + var_base),
         log = TRUE
      )
    
    #normalizing constant
    prob_norm = prob_clust[1:(K+1)]-
      logSumExp(prob_clust[1:(K+1)])
    
    #sampling new cluster assignments for j'th obs
    new_table = sample(1:(K+1),
                       1,
                       replace = FALSE,
                       exp(prob_norm))
    
    #new table addition 
    cluster[j] = new_table
    if (new_table == (K+1))
    {
      counts[K+1] = 1
      
      #sample the parameter from the posterior
      var_phi.2 = (var_base*var_obs)/
        (var_base+var_obs)
      mean_phi.2 = var_phi.2*
        (mu_base/var_base + 
           mixDat[cluster==(K+1)]/var_obs)
      phi.2[K+1] = rnorm(1, mean_phi.2, sqrt(var_phi.2))
      K = K+1
      
    } 
    else {
      counts[new_table] = counts[new_table]+1
    }
  }
  
  #sampling the parameter for the ii'th iteration
  sum_cluster = sapply(1:K, function(r) 
    sum(mixDat[cluster==r]))
  
  var_phi.2 = (var_base*var_obs)/
    (counts[1:K]*var_base+var_obs)
  
  mean_phi.2 = var_phi.2*
    (mu_base/var_base + 
       sum_cluster/var_obs)
  
  phi.2[1:K] = rnorm(K, mean_phi.2, sqrt(var_phi.2))
  
  #storing params value for observation 
  test_phi.2[ii] = phi.2[cluster[test_obs]]

} 


























