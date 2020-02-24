# HW4Q4.R
# Created: 2/24/2020
# Author:  Margaret Swift
# Contact: margaret.swift@duke.edu

# Homework 4 Question 4:
#   Simulate populations of size n = 5, 10, 20 individuals with 
#   positive inbreeding. Use 10 chains to see if you can recover 
#   parameter estimates. Is recovery better, worse, or the same 
#   as with negative inbreeding?

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LIBRARIES
library(rjags)
source('../clarkfunctions2020.R')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SCRATCH

file <- "logisReg.txt" 
n     <- 500
p     <- 2
x     <- matrix(rnorm(n*p),n,p)        # design matrix
x[,1] <- 1                             # intercept
beta  <- matrix(c(-3,3),p,1)           # parameter vector
ltg   <- x %*% beta                    # logit theta
tg    <- invlogit(ltg)                 # initial values for theta
y     <- rbinom(n,1,tg)                # data


cat("model{
  
  for (i in 1:n) {
    y[i] ~ dbern(theta[i])
    logit(theta[i]) <- b0 + b1*x1[i]
  }
  b0 ~ dnorm(0, .001)  
  b1 ~ dnorm(0, .001)
}", file = file)

logisData <- list(y = y, x1 = x[,2], n = length(y))
parNames <- c('b0','b1')

parInit <- function(){
  list(b0 = rnorm(1), b1 = rnorm(1))
}
jagsfit <- jags.model(data=logisData, inits=parInit,file=file)
update(jagsfit)
jagsLogit <- coda.samples(jagsfit, variable.names=c("b0","b1"), 
                          n.iter=5000)
tmp <- summary(jagsLogit)
print(tmp$statistics)


data  <- read.table('../dataFiles/FACEtrees.txt',header=T)
form  <- as.formula(cones ~ nfert*trt + diam)
X     <- model.matrix(form, data=data)
Y     <- model.frame(form, data=data)$cones
p     <- ncol(X)
n     <- nrow(X)
file <- "regModel.txt" 

cat("model{

  for (i in 1:n){
    Y[i] ~ dpois(lambda[i])                  
    lambda[i] <- exp( inprod(beta[],X[i,]) ) 
  }
  beta[1] ~ dnorm(0, .1)
  beta[2] ~ dunif(0, 10)
  beta[3] ~ dunif(0, 10)
  beta[4] ~ dnorm(0, .1)
  beta[5] ~ dnorm(0, .1)
}", file = file)

outjags <- jags.model(file = file, 
                      data = list(Y=Y,n=nrow(X),X = X))
update(outjags, 200)                # burnin 

outjags <- coda.samples(outjags, variable.names="beta", 
                        n.iter=2000)


tmp <- summary(outjags)
print(tmp$statistics)

traceplot(outjags)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# params
n <- c(5, 10, 20)
