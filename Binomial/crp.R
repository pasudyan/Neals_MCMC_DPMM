crp <- function(X, sz, alp = 1, G0 = c(1,1), num_iter = 1000)
{

  # CRP sampler for DP mixture model
  #   based on Algorithm 8 of Neal 2000
  #
  # Vinayak Rao
  #
  # X is a dataframe with columns X, eta, c and p_bin
  # Returns a list of two matrices giving cluster assignment and parameter of each observation

  empt      <- 5

  num_cls   <- max(X$c)
  num_obs   <- nrow(X)

  ret_c     <- matrix(rep(0,num_obs*num_iter), nrow=num_iter)
  ret_p_bin <- ret_c
  denom     <- num_obs + alp - 1

  n_c       <- rep(0, 100)  # Hopefully, 100 > typical num of clusters
  p_bin     <- rep(0, 100)
  p_vec     <- rep(0, 100)

  for(i in 1:num_cls) {
    n_c[i]   <- sum(X$c == i)
    p_bin[i] <- X$p_bin[X$c == i][1]     # cluster parameter
  }

  for(iter in 1:num_iter) {
    if (iter %% 100 == 0) print(iter)
    for(i in 1:num_obs) {

      X_curr      <- X$X[i]
      c_curr      <- X$c[i]
      eta_curr    <- X$eta[i]

      n_c[c_curr] <- n_c[c_curr] - 1

      if(n_c[c_curr] == 0) {
        flag          <- 1
        n_c[c_curr]   <- n_c[num_cls]
        n_c[num_cls]  <- 0
        tmp           <- p_bin[num_cls]; 
        p_bin[num_cls] <- p_bin[c_curr]; 
        p_bin[num_cls] <- tmp
        
        X$c[X$c == num_cls] <- c_curr
        num_cls       <- num_cls - 1 
        
      } else  {
        flag          <- 0 }

      c_empty      <- num_cls + seq_len(empt)
      
      p_vec[1:num_cls] <- log(n_c[1:num_cls])
      p_vec[c_empty]   <- log(alp) - log(empt)

      if(flag == 1) 
        p_bin[c_empty[-1]]   <- rbeta(empt-1, G0[1], G0[2]) 
      else
        p_bin[c_empty]       <- rbeta(empt, G0[1], G0[2])

      indx <- 1:(num_cls+empt)
      p_vec[indx]  <- p_vec[indx] + 
        dbinom(X_curr, sz, p_bin[indx]*eta_curr, log=T)
      p_vec[indx]  <- p_vec[indx] - max(p_vec[indx])
      p_vec[indx]  <- exp(p_vec[indx])
      p_vec[indx]  <- p_vec[indx] / sum(p_vec[indx])

      c_new        <- sample(num_cls+empt, 1, FALSE, p_vec[indx])

      if(c_new > num_cls) {
        p_bin[num_cls+1] <- p_bin[c_new]
        c_new            <- num_cls+1
        num_cls          <- num_cls + 1
      }
      X$c[i]      <- c_new
      n_c[c_new]  <- n_c[c_new] + 1
   }

   for(j in 1:num_cls) {
     indx          <- X$c == j
     
     llik          <- function(p_bin) {
       sum(dbinom(X$X[indx], sz, p_bin*X$eta[indx], log=T)) + 
         dbeta(p_bin, G0[1], G0[2], log=T)
     }
     
     p_bin[j]      <- as.vector(uni.slice(p_bin[j], 
                                          llik, 
                                          w=1, 
                                          m=Inf, 
                                          lower=0, 
                                          upper=1))

     ret_c[iter,indx]     <- j
     ret_p_bin[iter,indx] <- p_bin[j] 
   }
  }
  return(list(c=ret_c, p_bin=ret_p_bin))
}
