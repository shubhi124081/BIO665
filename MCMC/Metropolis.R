# HW4.R
# Created: 2/21/2020
# Author:  Margaret Swift
# Contact: margaret.swift@duke.edu

# Homework 4 exercises

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOADING
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('clarkFunctions2020.R')
# Rcpp::sourceCpp('../clarkCppFns.cpp')


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCTIONS

getMu <- function(sigma, y, M, m){ # conditional distributions
  
  n <- length(y)
  V <- 1/( n/sigma + 1/M )
  v <- sum(y)/sigma + m/M
  rnorm(1, V*v, sqrt(V))
}
getSigma <- function(mu, y, s1, s2){
  
  n  <- length(y)
  u1 <- s1 + n/2
  u2 <- s2 + 1/2*sum( (y - mu)^2 )
  1/rgamma(1, u1, u2)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Gibbs sampling
par(bty='n')

n  <- 10                          # simulate data
s  <- 1
mu <- -2
y  <- rnorm(n, mu, sqrt(s))

m  <- 0                            # prior parameter values
M  <- 1
s1 <- s2 <- 1

mg <- 0
sg <- 2
G  <- 10

plot(mg, sg, xlim= mg + c(-2.5,2.5) ,ylim = sg + c(-2.5,2.5), 
     xlab='mean', ylab='variance')

# one step at a time

for(g in 1:G){
  mn <- getMu(sg, y, M, m);      lines(c(mg,mn),c(sg,sg),col='blue')
  sn <- getSigma(mn, y, s1, s2); lines(c(mn,mn),c(sg,sn),col='green')
  text(mn,sn,g)
  mg <- mn
  sg <- sn
  readline('return to continue')
}


# Chains
par(mfcol=c(2,2),mar=c(4,4,1,1),bty='n')

G <- 1000
chains <- matrix(NA,G,2)

for(g in 1:G){
  mg <- getMu(sg, y, M, m)
  sg <- getSigma(mg, y, s1, s2)
  chains[g,] <- c(mg, sg)
}
qc <- apply(chains,2,quantile,c(.025,.975))
ml <- 'mean parameter'
vl <- 'variance parameter'
plot(chains[,1],chains[,2],cex=.3, xlab=ml, ylab=vl, log='y')
abline(v=qc[,1],lty=2); abline(h=qc[,2],lty=2)
plot(chains[,1],c(1:G),type='l', xlab=ml, ylab='Index')
abline(v=qc[,1],lty=2)
plot(chains[,2],type='l', ylab=vl, log='y')
abline(h=qc[,2],lty=2)






#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Exercise 1
#   Summarize the estimates as posterior mean colMeans(chains), median apply(chains,2,median), 
#   credible interval apply(chains,2,quantile,c(.025,.975)), and parameter correlation 
#   cor(chains).





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


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Exercise 3: Metropolis Sampler



#eof