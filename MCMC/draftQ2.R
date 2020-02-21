#posterior is poisson 

n <- 10 #number of samples 
#since prior is uniform 

y <- rpois(n, lambda)


lambda <- 2 #mean and variance of posterior + likelihood 
LAMBDA <- NULL

lambda0 <- 1 #random 

mu <- lambda0 #mu is the mean of the proposal normal 
sig <- 1 #sigma is the sd of the proposal 
nrep <- 1000
for(i in 1:nrep){
  
  #simulate from our normal 
  lambdaStar <- rnorm(1, mu, sig)

#get ratio r in Log terms 
logR <- log(sum(dpois(y, lambdaStar)+ dunif(lambdaStar, min = 0, max = 10)) - 
              sum(dpois(y, lmabda0) + dunif(lambda0, min = 0, max = 10)))

if(log(runif(1))< logR) { lambda0 <-lambdaStar }
LAMBDA <- c(LAMBDA, lambdaStar)
  
}





















