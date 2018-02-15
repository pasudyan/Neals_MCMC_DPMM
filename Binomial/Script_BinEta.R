### Script to run the Algorithms ###
####################################

rm(list=ls())

##Loading Library##
library(ggplot2)
library(matrixStats)

##Parameter initialization##
num_cust  = 20
alpha     = 2
n         = c(5, 25, 50, 100, 150)

alpha_base = 1 #parameters of base distribution
beta_base  = 1
m          = 2 #number of auxiliary params

#test_obs  = sample((num_cust/2):num_cust,1,replace=F)
test_obs   = 2

#Generating Mixed Data
zeta      = .9
theta     = c(rep(1, num_cust/2), rep(1, num_cust/2))
trials    = numeric(num_cust)
trials    = rep(n, length.out=num_cust)

mixDat    = rbinom(num_cust, trials, zeta*theta)


dataMat = data.frame(mixDat = mixDat, 
                     zeta   = zeta,
                     theta  = theta,
                     pi     = zeta*theta,
                     trials = trials)

numIter = 5000 #number of observations
burnIn = 2000 #number of burnIn
telliter = 500 #number of iteration to print
xtraSpace = 2*num_cust
source("./mh_algo.R")
source("./mh_algoZeta.R")

# source("./Algo8MH3_BinEta.R")

source("./Algo8_BinZeta.R")

## Plotting for all 4 Algorithms 
test_params = data.frame(alg8MH3 = test_phi.8MH3[(burnIn+1):numIter],
                         zeta    = test_zeta[(burnIn+1):numIter])  

## Plotting the posterior densit of the Algorithms 
p = ggplot(test_params) + 
  geom_density((aes(x=alg8MH3, color="8MH3"))) + 
  ggtitle("Posterior Density of Theta") + 
  scale_y_continuous("Density values")+
  xlab("Parameter values") +
  labs(color="Algorithm")
print(p)

print(summary(test_params$alg8MH3))

q = ggplot(test_params) + 
  geom_density((aes(x=zeta, color="Zeta"))) + 
  ggtitle("Posterior Density of Zeta") + 
  scale_y_continuous("Density values")+
  xlab("Parameter values") +
  labs(color="Algorithm")
print(q)

print(summary(test_params$zeta))
