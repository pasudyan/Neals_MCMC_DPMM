##Chinese restaurant process: Beta base distribution#
#########################################################

crpBeta <- function(num_cust, alpha, trials, alpha_base, beta_base){
  if (num_cust == 0){
    
    print("Number of customers must be greater than 0")
    
  }  else if (alpha == 0){
    
    print("Alpha must be positive integers")
    
  } else{
    
    table = numeric(num_cust)
    table[1] = 1
    next_table = 2
    
    for (i in 2:num_cust){
      if (runif(1) < alpha/(i-1+alpha)){
        
        table[i] = next_table
        next_table = next_table + 1
        
      } else {
        
        size = sapply(1:(next_table-1), 
                      function(r) sum(table == r))
        which_table = sample(1:(next_table-1),
                             1, 
                             replace = FALSE, 
                             (1/(i-1+alpha))*size)
        table[i] = which_table
      }
    }
  }
  
  params = rbeta((next_table-1), alpha_base, beta_base)

  val = list(
    parameters = params,
    values = params[table], 
    tables = table
  )
  
  return(val)
}