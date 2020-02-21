# HW4.R
# Created: 2/21/2020
# Author:  Margaret Swift
# Contact: margaret.swift@duke.edu

# Homework 4 exercises

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Exercise 2: Metropolis Sampler
n <- 100
lambda <- 4
y <- rpois(n, lambda)

# Metropolis algorithm
nrep <- 10000
Lcur <- runif(1, 0, 10)

# Loop over a number of repetitions
for (i in 1:nrep) {
  # Propose a new lambda
  Lstar <- .tnorm(1, 0, 10, mu=Lcur, sig=2)

  # Find probability of each lambda given y
  num   <- prod(dpois(y, Lstar) * dunif(Lstar, 0, 10))
  denom <- prod(dpois(y, Lcur ) * dunif(Lcur, 0, 10))
  ratio <- num / denom
  
  # Accept or reject:
  u <- runif(1, 0, 1)
  if (u < ratio) Lcur <- Lstar
}
Lcur



#eof