##MCMC methods for parameter estimation (Algorithm 8)##
#######################################################

test_phi.8MH = numeric(numIter)

#initialize number of clusters
K = floor(alpha*log(num_cust))

#initialize cluster assignment
phi.8MH = numeric(xtraSpace)
phi.8MH[1:K] = rnorm(K, mu_base, sd_base) #for algo 8

#storing the cluster assignments for each iteration
cluster_store8 = numeric(numIter)

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

mh_algo <- function(datta, phi.8MH, cluster){
  nuCurr  = phi.8MH
  ratio   = numeric(length(nuCurr))
  params  = numeric(length(nuCurr))
  mu_base = 0 
  sd_obs  = 1 
  sd_base = 10 
  var_obs = sd_obs^2
  var_base  = sd_base^2
  sd_prop   = 1
  
  counts = sapply(1:length(nuCurr), function(r)
    sum(cluster == r))

  for(t in 1:5){
    nuCand = rnorm(length(nuCurr),nuCurr,sd_prop)
    
    #Calculating the acceptance probability 
    for (i in 1:length(nuCurr)){
      mydatta = datta[cluster==i]
      
      var_post = var_base*var_obs/
        (counts[i]*var_base + var_obs)
      
      mu_post = var_post*(mu_base/var_base + 
                           sum(mydatta)/var_obs)
      
      n_vec = sum(dnorm(mydatta, nuCand[i], sqrt(var_obs), log=TRUE), 
                  dnorm(nuCand[i], mu_base, sd_base, log=TRUE))
      
      d_vec = sum(dnorm(mydatta, nuCurr[i], sqrt(var_obs), log=TRUE), 
                  dnorm(nuCurr[i],  mu_base, sd_base, log=TRUE))

      ratio[i] = sum(n_vec-d_vec)
      
      if (exp(ratio[i]) > 1){
        params[i] = nuCand[i]
      } else if (runif(1) < exp(ratio[i])){
        params[i] = nuCand[i]
      } else {
        params[i]= nuCurr[i]
      }
    }
    nuCurr = params
  }
  return(nuCurr)
}

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
      phi.8MH[(K+1):h] = rnorm(m, mu_base, sd_base) 
      
    } else if (counts[cluster[j]]==0 & cluster[j]!=K){
      
      #renumbering the clusters
      counts[cluster[j]]  = counts[K]
      cluster[cluster==K] = cluster[j]
      counts[K] = 0
      
      #renumbering the phi.8MH
      temp = phi.8MH[cluster[j]]
      phi.8MH[cluster[j]] = phi.8MH[K]
      phi.8MH[K] = temp
      
      #reducing the number of clusters and aux clusters
      K = K-1
      h = h-1
      
      #drawing values for the aux params from G_0
      phi.8MH[(K+2):h] = rnorm(m-1, mu_base, sd_base)
      
    } else if (counts[cluster[j]]==0 & cluster[j]==K){
      
      #reducing the number of clusters and aux clusters
      K = K-1
      h = h-1
      
      #drawing values for the aux params from G_0
      phi.8MH[(K+2):h] = rnorm((m-1), mu_base, sd_base)
      
    }
    
    #prob of choosing existing cluster 
    prob_clust[1:K] = log(counts[1:K]) +
      dnorm(mixDat[j], phi.8MH[1:K], sd_obs, log=TRUE)
    
    #prob of choosing new cluster
    prob_clust[(K+1):h] = log(alpha/m) + 
      dnorm(mixDat[j], phi.8MH[(K+1):h], sd_obs, log=TRUE)
    
    #normalizing constant
    prob_norm = prob_clust[1:h] - 
      logSumExp(prob_clust[1:h])
    
    #sampling new cluster assignments
    new_table = sample(1:h, 
                       1, 
                       replace = FALSE,
                       exp(prob_norm)
    )
    
    #new table addition
    if (new_table > K){
      
      cluster[j] = K+1
      phi.8MH[K+1] = phi.8MH[new_table]
      counts[K+1] = 1
      K = K+1
      
    } else {
      
      cluster[j] = new_table
      counts[new_table] = counts[new_table]+1
      
    }
    
    phi.8MH[(K+1):xtraSpace]  = rep(0,xtraSpace-K)
    counts[(K+1):xtraSpace] = rep(0,xtraSpace-K)
    
  } 
  
  #sampling the new parameters
  #this one here could be done using the MH algorithm
  
  phi.8MH[1:K] = mh_algo(mixDat, 
                         phi.8MH[1:K],
                         cluster)
  
  #storing params value for observation 
  test_phi.8MH[ii] = phi.8MH[cluster[test_obs]]
  
}

