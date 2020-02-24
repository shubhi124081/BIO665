# LIBRARIES
rm(list=ls())
library(rjags)
source('../clarkfunctions2020.R')
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

n  <- 10                          # simulate data
s  <- 1
mu <- -2
y  <- rnorm(n, mu, sqrt(s))

m  <- 0                            # prior parameter values
M  <- 1
s1 <- s2 <- 1

mg <- 0
sg <- 2

G <- 1000
chains <- matrix(NA,G,2)

for(g in 1:G){
  mg <- getMu(sg, y, M, m)
  sg <- getSigma(mg, y, s1, s2)
  chains[g,] <- c(mg, sg)
}
qc <- apply(chains,2,quantile,c(.025,.975))

# Z's code
summ <- data.frame(mean = colMeans(chains), 
                   median = apply(chains,2,median), 
                   Quantile1 = apply(chains,2,quantile,c(.025,.975))[1,],
                   Quantile2 = apply(chains,2,quantile,c(.025,.975))[2,],
                   cor_Mean = cor(chains)[1,],
                   cor_Variance = cor(chains)[2,])

row.names(summ)[1] <- c("Mean")
row.names(summ)[2] <- c("Variance")


knitr::kable(
  summ, 
  caption = "Summary for Chains of Mean and Variance Parameter"
)
