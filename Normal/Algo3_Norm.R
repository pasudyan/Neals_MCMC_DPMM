##MCMC methods for parameter estimation (Algorithm 3)##
#######################################################

test_phi.3 = numeric(numIter)

#initialize number of clusters
K = floor(alpha*log(num_cust)) 

#initialize parameters for each cluster
phi.3 = numeric(xtraSpace) #for algo 3

#storing the cluster assignments for each iteration
cluster_store3 = numeric(numIter)

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
    
    #calculating parameters for table probabilities
    sum_cluster = sapply(1:K, function(r) 
      sum(mixDat[cluster==r]))
    
    var_post = (var_base*var_obs)/
      (counts[1:K]*var_base+var_obs)
    
    mean_post = var_post*(mu_base/var_base +
                            sum_cluster/var_obs)
    
    #prob of choosing existing table
    prob_clust[1:K] = log(counts[1:K]) + 
      dnorm(mixDat[j],
            mean_post, 
            sqrt(var_post + var_obs),
            log = TRUE
            )
    
    #prob of choosing new table
    prob_clust[K+1] = log(alpha) + 
      dnorm(mixDat[j], 
            mu_base, 
            sqrt(var_base + var_obs),
            log = TRUE
            )
    
    #normalizing constant
    prob_norm = prob_clust[1:(K+1)]-
      logSumExp(prob_clust[1:(K+1)])
    
    #sampling new cluster assignments
    new_table = sample(1:(K+1),
                       1,
                       replace = FALSE,
                       exp(prob_norm)
                       )
    
    #new table addition 
    cluster[j] = new_table
    if (new_table == (K+1)){
      counts[K+1] = 1
      K = K+1
    } else {
      counts[new_table] = counts[new_table]+1
    }
  }
  
  #sampling the parameter for each cluster 
  #from the posterior
  sum_cluster = sapply(1:K, function(r) 
    sum(mixDat[cluster==r]))
  
  var_phi.3  = (var_base*var_obs)/
    (counts[1:K]*var_base+var_obs)
  
  mean_phi.3 = var_phi.3*
    (mu_base/var_base + 
       sum_cluster/var_obs)
  
  phi.3[1:K] = rnorm(K, mean_phi.3, sqrt(var_phi.3))
  
  #storing params value for observation 
  test_phi.3[ii] = phi.3[cluster[test_obs]]
  
}




