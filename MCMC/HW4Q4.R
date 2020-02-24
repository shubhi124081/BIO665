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
rm(list=ls())
library(rjags)
source('../clarkfunctions2020.R')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SCRATCH

data  <- read.table('FACEtrees.txt',header=T)
form  <- as.formula(cones ~ nfert*trt + diam)
X     <- model.matrix(form, data=data)
Y     <- model.frame(form, data=data)$cones
n     <- nrow(X)
p     <- ncol(X)
# beta  <- matrix(c(-3,3),p,1)           # parameter vector
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
# p = freq of allele A; f = inbreeding coef
# pmake out = c(paa, pab, pbb)

pmake(c(.5,0)) #no inbreeding

n <- c










#eof