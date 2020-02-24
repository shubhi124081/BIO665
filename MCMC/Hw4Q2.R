
# y
n <- 100
lambda <- 3
y <- rpois(n, lambda)

# proposed lambda 
L0 <- 8
sig <- 3  # variance for the proposal distribution
nrep <- 100
LAMBDA <- NULL
par(mfrow=c(2,3))

for(sig in c(0.05,0.5,1,3,6,10)){
  for (i in 1:nrep){
    Lstar <- .tnorm(1,0,10,L0,sig)
    r <- prod(dpois(y,Lstar)*dunif(Lstar, min = 0, max = 10))/prod(dpois(y,L0) *dunif(L0, min = 0, max = 10))
    
    # reject or accept 
    if (r > runif(1,0,1)) L0 <- Lstar
    LAMBDA <- c(LAMBDA, L0)
  }
  plot(density(LAMBDA), main = paste("sigma = ", sig), xlim = c(0,10))
  LAMBDA <- NULL
}

# looks like proposal distribution with higher sigma are more efficient than those with smaller sigma. 
