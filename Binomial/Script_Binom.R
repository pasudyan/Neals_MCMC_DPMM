#########################################################  
## Script to estimate parameters using	               ##
## Algorithm 1, 2, and 8 of Neal's 2000 paper          ##
##                                                     ## 
## Observations: Y_i ~ Bin(N, theta_i)                 ##
## where theta_i has prior Beta(1,1)                   ## 
#########################################################

rm(list=ls())

##Loading Library##
library(ggplot2)
library(matrixStats)

##Generating data from the Gaussian Mixtures##
source("./crpBinom.R")

##Parameter initialization##
num_cust  = 50
alpha     = 2
n         = 10

alpha_base = 1 #parameters of base distribution
beta_base  = 1
m          = 2 #number of auxiliary params

test_obs  = sample((num_cust/2):num_cust,1,replace=F)

#Generating Mixed Data
# dat     = crpBinom(num_cust, alpha, n, alpha_base, beta_base)
# mixDat  = dat$dirichBinom
mixDat  = c(rbinom(num_cust/2, n, .1), rbinom(num_cust/2,n, .90))

numIter = 5000 #number of observations
burnIn = 2000 #number of burnIn
telliter = 500 #number of iteration to print
xtraSpace = 2*num_cust

print("Starting Algorithm 2")
source("./Algo2_Binom.R")
print("Starting Algorithm 3")
source("./Algo3_Binom.R")
print("Starting Algorithm 8")
source("./Algo8_Binom.R")

#test for equal distribution
ksTest <- function(x, y){
  ks.test(x,y,alternative="two.sided")
}

test_phi.2=test_phi.2[(burnIn+1):numIter]
test_phi.3=test_phi.3[(burnIn+1):numIter]
test_phi.8=test_phi.8[(burnIn+1):numIter]

ksTest(test_phi.2, test_phi.3)
ksTest(test_phi.2, test_phi.8)
ksTest(test_phi.3, test_phi.8)

## Plotting the density of the parameters
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

print(summary(test_params3$alg2))
print(summary(test_params3$alg3))
print(summary(test_params3$alg8))








