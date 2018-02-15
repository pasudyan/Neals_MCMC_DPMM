### Script to run the Algorithms ###
####################################

rm(list=ls())

##Loading Library##
library(ggplot2)
library(matrixStats)

##Generating data from the Gaussian Mixtures##
source("./crpNorm.R")
num_cust = 15
alpha = 2

test_obs = sample(1:num_cust,1,replace=F)

# dat <- crpNorm(num_cust,alpha)
# mixDat <- dat$dirichNorm
mixDat = c(rnorm(5,-10,1), rnorm(5,0,1), rnorm(5,10,1))
mixDat = sample(mixDat, num_cust, replace=FALSE)

##Parameter initialization##
mu_base = 0 #mean of base distribution
sd_obs = 1 #sd obs
sd_base = 10 #sd of base distribution
var_obs = sd_obs^2
var_base = sd_base^2
m = 2 #number of auxiliary params
xtraSpace = num_cust

numIter = 5000 #number of observations
burnIn = 2000 #number of burnIn
telliter = 500 #number of iteration to print

##Predictive Distribution##
points = seq(-25, 25, by = .1) #points for plotting
len = length(points)
runs = (numIter-burnIn)

##Predictive probability allocation##
val_alg3 = matrix(rep(0, runs*len), nrow=runs, ncol = len)
val_alg2 = matrix(rep(0, runs*len), nrow=runs, ncol = len)
val_alg8 = matrix(rep(0, runs*len), nrow=runs, ncol = len)
val_alg8MH = matrix(rep(0, runs*len), nrow=runs, ncol = len)

print("Starting Algorithm 2")
source("./Algo2_NormPred.R")
print("Starting Algorithm 3")
source("./Algo3_NormPred.R")
print("Starting Algorithm 8")
source("./Algo8_NormPred.R")
print("Starting Algorithm 8MH")
source("./Algo8MH_NormPred.R")

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
ksTest(test_phi.8, test_phi.8MH)

#Plotting the distribution of the posterior
test_params = data.frame(alg2 = test_phi.2, 
                         alg3 = test_phi.3, 
                         alg8 = test_phi.8,
                         alg8MH = test_phi.8MH)

p = ggplot(test_params) + 
  geom_density(aes(x=alg2, color="2")) + 
  geom_density(aes(x=alg3, color="3")) + 
  geom_density(aes(x=alg8, color="8")) +
  geom_density(aes(x=alg8MH, color="8MH")) +
  ggtitle("Posterior Density of Normal Parameter") + 
  scale_y_continuous("Density values") + 
  xlab("Parameter values") + 
  labs(color="Algorithm")
print(p)

#plotting predictive distribution
datDensity <- data.frame(
  cond = factor(rep("A",each=length(mixDat))), 
  mixDat=mixDat
  )

geom_rug <- function(data, mapping, ...){
  mapping$y = 0
  geom_point(data = data, mapping = mapping, shape=I("|"), ...)
}

g <- ggplot(datPoints2, aes(x=points))
g <- g + geom_line(data=datPoints2, aes(y=mean, color="2"))
g <- g + geom_ribbon(data=datPoints2, 
                     aes(ymin=lower,ymax=upper,fill="2"),
                     alpha=0.1)

g <- g + geom_line(data=datPoints3, aes(y=mean, color="3"))
g <- g + geom_ribbon(data=datPoints3, 
                     aes(ymin=lower,ymax=upper,fill="3"),
                     alpha=0.1) 

g <- g + geom_line(data=datPoints8, aes(y=mean, color="8"))
g <- g + geom_ribbon(data=datPoints8, 
                     aes(ymin=lower,ymax=upper,fill="8"),
                     alpha=0.1)

g <- g + geom_line(data=datPoints8MH, aes(y=mean, color="8MH"))
g <- g + geom_ribbon(data=datPoints8MH, 
                     aes(ymin=lower,ymax=upper,fill="8MH"),
                     alpha=0.1)

g <- g + geom_rug(data = datDensity, 
                  mapping = aes(x=mixDat),
                  alpha=0.8)

g <- g + xlab("Points") + 
  ggtitle("Predictive Density for Normal Distribution") +
  scale_y_continuous("Density values") + 
  labs(color="Algorithm", fill="Algorithm")
print(g)
