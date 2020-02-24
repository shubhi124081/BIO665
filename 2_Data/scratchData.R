# scratchData.R
# Created: 1/29/20
# Author:  Margaret Swift
# Contact: margaret.swift@duke.edu

# Simulate data and estimate betas

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOADING
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('../clarkFunctions2020.R')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCTIONS
genData <- function(n, Q, sigma) {
  # Parameters
  epsilon <- rnorm(n, 0, sigma)
  beta <- matrix( rnorm(Q), Q, 1 )
  
  # Simulate x and y
  x <- matrix( rnorm(n*Q), n, Q )
  x[,1] <- 1 #set intercept
  y <- x%*%beta + epsilon
  
  # Return data
  data <- list(x=x, y=y, beta=beta, n=n, Q=Q, sigma=sigma, epsilon=epsilon)
  return(data)
}
gibbs <- function(d, nstep=1000) {
  # Gibbs sampler
  Q=d$Q; x=d$x; y=d$y; n=d$n; sigma=d$sigma
  
  # Set up to find mu
  b <- matrix(0, Q, 1)
  B <- diag(1000, Q) #priors are horrible, so I give it a wide variance
  B.inv <- solve(B)
  
  for (i in 1:nstep) {
    # Find mu and standard errors
    V.inv <- 1/sigma/sigma * crossprod(x) + B.inv
    V <- solve(V.inv) #covariance matrix of beta
    v <- 1/sigma/sigma * crossprod(x,y) + B.inv%*%b
    sderrors <- sqrt(diag(V))
    mu <- V%*%v
    
    # And generating random values from this distn we know now
    bstar <- .rMVN(1000, mu, V)
    
    # new sigma
    u1 <- u2 <- 1
    # sigma <- 1/rgamma(1, u1+(n/2), u2+0.5*crossprod(y-x%*%beta))
    # but we're using beta here... should I be using mu then?
    sigma2 <- 1/rgamma(1, u1+(n/2), u2+0.5*crossprod(y-x%*%mu))
    sigma <- sqrt(sigma2)
  }
  return(colMeans(bstar))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MAIN
n=100; Q=5; sigma=0.1
data <- genData(n=n, Q=Q, sigma=sigma)
newBeta <- gibbs(data)






#EOF