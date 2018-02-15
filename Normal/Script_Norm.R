### Script to run the Algorithms ###
####################################

rm(list=ls())

##Loading Library##
library(ggplot2)
library(matrixStats)

##Generating data from the Gaussian Mixtures##
source("./crpNorm.R")
num_cust = 10
alpha = 2

test_obs = sample(1:num_cust,1,replace=F)

dat <- crpNorm(num_cust,alpha)
mixDat <- dat$dirichNorm
# mixDat = c(rnorm(5,1,1), rnorm(5,20,1))

##Parameter initialization##
mu_base = 0 #mean of base distribution
sd_obs = 1 #sd obs
sd_base = 10 #sd of base distribution
var_obs = sd_obs^2
var_base = sd_base^2
m = 2 #number of auxiliary params
xtraSpace = 20

numIter = 5000 #number of observations
burnIn = 3000 #number of burnIn
telliter = 500 #number of iteration to print

##Predictive Distribution##
points = seq(-30, 30, by = .1) #points for plotting
len = length(points)
runs = (numIter-burnIn)

##Predictive probability allocation##
val_alg3 = matrix(rep(0, runs*len), nrow=runs, ncol = len)
val_alg2 = matrix(rep(0, runs*len), nrow=runs, ncol = len)
val_alg8 = matrix(rep(0, runs*len), nrow=runs, ncol = len)
val_alg8MH = matrix(rep(0, runs*len), nrow=runs, ncol = len)

print("Starting Algorithm 2")
source("./Algo2_Norm.R")
print("Starting Algorithm 3")
source("./Algo3_Norm.R")
print("Starting Algorithm 8")
source("./Algo8_Norm.R")
print("Starting Algorithm 8MH")
source("./Algo8MH_Norm.R")

test_phi.2=test_phi.2[(burnIn+1):numIter]
test_phi.3=test_phi.3[(burnIn+1):numIter]
test_phi.8=test_phi.8[(burnIn+1):numIter]
test_phi.8MH=test_phi.8MH[(burnIn+1):numIter]

#test for equal distribution
ksTest <- function(x, y){
  ks.test(x,y,alternative="two.sided")
}

## KS Test for the distribution of the parameters 
ksTest(test_phi.2, test_phi.3)
ksTest(test_phi.2, test_phi.8)
ksTest(test_phi.3, test_phi.8)
ksTest(test_phi.2, test_phi.8MH)

#Plotting the distribution of the posterior
test_params3 = data.frame(alg2 = test_phi.2, alg3 = test_phi.3, alg8 = test_phi.8)

p = ggplot(test_params3) + geom_density(aes(x=alg2, color="2")) + geom_density(aes(x=alg3, color="3")) + geom_density(aes(x=alg8, color="8")) + ggtitle("Posterior Density of Normal Parameter") + scale_y_continuous("Density values") + xlab("Parameter values") + labs(color="Algorithm")
print(p)

#Plotting the distribution of the posterior
test_params4 = data.frame(alg2 = test_phi.2,
                         alg3 = test_phi.3,
                         alg8 = test_phi.8,
                         alg8MH = test_phi.8MH)

p = ggplot(test_params4) + 
  geom_density(aes(x=alg2, color="2")) +
  geom_density(aes(x=alg3, color="3")) +
  geom_density(aes(x=alg8, color="8")) +
  geom_density(aes(x=alg8MH, color="8MH")) +
  ggtitle("Posterior Density of Normal Parameter") +
  scale_y_continuous("Density values") +
  xlab("Parameter values") +
  labs(color="Algorithm")
print(p)







