#########################################################  
## Script to estimate parameters using                 ##
## Algorithm 1, 2, 8, and 8-MH of Neal's 2000 paper    ##
##                                                     ## 
## Observations: Y_i ~ Bin(N_i, theta_i)               ##
## where theta_i has prior Beta(1,1)                   ## 
#########################################################

rm(list=ls())

##Loading Library##
library(ggplot2)
library(matrixStats)

##Generating data from the Gaussian Mixtures##
source("./crpBinom.R")

##Parameter initialization##
num_cust  = 20
alpha     = 2
n         = c(5, 25, 50, 100, 150)

alpha_base = 1 #parameters of base distribution
beta_base  = 1
m          = 2 #number of auxiliary params

test_obs  = sample((num_cust/2):num_cust,1,replace=F)

#Generating Mixed Data
# dat     = crpBinom(num_cust, alpha, n, alpha_base, beta_base)
# mixDat  = dat$dirichBinom
mixDat  = c(rbinom(num_cust/2, n, .1), rbinom(num_cust/2, n, .90))

trials = numeric(num_cust)
trials = rep(n, length.out=num_cust)

numIter = 7000 #number of observations
burnIn = 2000 #number of burnIn
telliter = 500 #number of iteration to print
xtraSpace = 2*num_cust

print(sprintf("Number of Samples: %d",num_cust))
print("Starting Algorithm 2")

ptm = proc.time()
source("./Algo2_Binom2.R")
print(proc.time() - ptm) 

print("Starting Algorithm 3")
ptm = proc.time()
source("./Algo3_Binom2.R")
print(proc.time() - ptm)

print("Starting Algorithm 8")
ptm = proc.time()
source("./Algo8_Binom2.R")
print(proc.time() - ptm)

print("Starting Algorithm 8MH")
ptm = proc.time()
source("./Algo8MH_Binom2.R")
print(proc.time() - ptm)

print("Starting Algorithm 8MH3")
ptm = proc.time()
source("./Algo8MH_Binom3.R")
print(proc.time() - ptm)

#test for equal distribution
ksTest <- function(x, y){
  ks.test(x,y,alternative="two.sided")
}

test_phi.2=test_phi.2[(burnIn+1):numIter]
test_phi.3=test_phi.3[(burnIn+1):numIter]
test_phi.8=test_phi.8[(burnIn+1):numIter]
test_phi.8MH=test_phi.8MH[(burnIn+1):numIter]
test_phi.8MH3=test_phi.8MH3[(burnIn+1):numIter]

print(ksTest(test_phi.2, test_phi.3))
print(ksTest(test_phi.3, test_phi.8))
print(ksTest(test_phi.8, test_phi.8MH))
print(ksTest(test_phi.8MH3, test_phi.8MH))

# Plotting the density of the parameters
test_params3 = data.frame(alg2 = test_phi.2, 
                          alg3 = test_phi.3,
                          alg8 = test_phi.8)  

## Plotting for all 3 Algorithms 
p = ggplot(test_params3) + 
  geom_density(aes(x=alg2, color="2")) + 
  geom_density(aes(x=alg3, color="3")) + 
  geom_density(aes(x=alg8, color="8")) + 
  ggtitle("Posterior Density of Binomial Parameter") + 
  scale_y_continuous("Density values") + xlab("Parameter values") + labs(color="Algorithm")
print(p)

## Plotting for all 4 Algorithms 
test_params4 = data.frame(alg2 = test_phi.2, 
                          alg3 = test_phi.3,
                          alg8 = test_phi.8,
                          alg8MH = test_phi.8MH,
                          alg8MH3 = test_phi.8MH3)  

## Plotting for all 3 Algorithms 
p = ggplot(test_params4) + 
  geom_density(aes(x=alg2, color="2")) + 
  geom_density(aes(x=alg3, color="3")) + 
  geom_density(aes(x=alg8, color="8")) +
  geom_density(aes(x=alg8MH, color="8MH")) +
  geom_density((aes(x=alg8MH3, color="8MH3"))) + 
  ggtitle("Posterior Density of Binomial Parameter") + 
  scale_y_continuous("Density values") + xlab("Parameter values") + labs(color="Algorithm")
print(p)

print(summary(test_params4$alg2))
print(summary(test_params4$alg3))
print(summary(test_params4$alg8))
print(summary(test_params4$alg8MH))
print(summary(test_params4$alg8MH3))

test_params = data.frame(alg8 = test_phi.8,
                         alg8MH = test_phi.8MH,
                         alg8MH3 = test_phi.8MH3)

q = ggplot(test_params) + 
  geom_density(aes(x=alg8, color="8")) +
  geom_density(aes(x=alg8MH, color="8MH")) +
  geom_density((aes(x=alg8MH3, color="8MH3"))) + 
  ggtitle("Posterior Density of Binomial Parameter") + 
  scale_y_continuous("Density values") + xlab("Parameter values") + labs(color="Algorithm")
print(q)

print(summary(test_params$alg8))
print(summary(test_params$alg8MH))
print(summary(test_params$alg8MH3))


