# regressionExample.R
# 1/13/20
# JS Clark
# Regression example for Bayes class (Intro)
# expand.grid!!!!!!!!!!!!
# image!!!!!!!!!!!!!!!!!!

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LIBRARIES
require(mvtnorm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCTIONS

# Log-likelihood
logLikelihood <- function(y, mu, sigma, log=T){ 
  return(sum( dnorm(y, mu, sigma, log) ))
}

# Surface for betas
betaSurface <- function(X, y, sigma, nx = 50, ny = 50, bmu = NULL, bvar = NULL){
  
  # Set bounds of the surface to encapsulate most of ~N(m,s)
  b0 <- seq(-3, 3, length=nx)
  b1 <- seq(-3, 3, length=ny)
  bsurf <- as.matrix( expand.grid(b0, b1) )
  
  # mean and likelihood surfaces
  msurf <- X%*%t(bsurf) #expected y for all parameter combinations
  lsurf <- dnorm(y, msurf, sigma, log = T)  # logLikelihood
  lsurf <- colSums(lsurf)
  
  # Find the maximum likelihood surface & turn lsurf into matrix
  bmax  <- which.max(lsurf)
  z <- matrix(lsurf, nx, ny)
  
  # Show an image of the beta surface.
  image(b0, b1, z)
  abline(v=bsurf[bmax,1])
  abline(h=bsurf[bmax,2])
  
  # what to return.
  out <- list( betaMu = matrix(bsurf[bmax,], 2), b0 = b0, b1 = b1, lsurf = z )
  
  # from prior, if beta-mu is given
  if(!is.null(bmu)){
    psurf <- dmvnorm(bsurf, bmu, bvar, log=T)
    zp <-  matrix(psurf, nx, ny)
    contour(b0, b1, zp, add=T)
    out <- append( out, list(psurf = zp) )
  }
  return( out )
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Simulate regression data

# scalars
n     <- 200
Q     <- 2
sigma <- 0.1

# Betas and design matrix (X)
beta  <- matrix( rnorm(Q), nrow=2 )
X     <- matrix( rnorm(n*Q), nrow=n, ncol=Q )
X[,1] <- 1 #start at 1 to make it the intercept

# randomly distribute responses; %*% = matrix multiplication! :)
y <- rnorm(n, X%*%beta, sigma)

# Plot y's by x's
plot(X[,2], y)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Beta surface

blist   <- betaSurface(X, y, sigma)
betaHat <- blist$betaMu
mu      <- X%*%betaHat

# Let's see how much better/worse our max likelihood is, compared to 
# the likelihood at (0,0).
logLikelihood(y, mu, sigma, log = T)
logLikelihood(y, mu*0, sigma, log = T)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Bayes part

# Now let's add in the prior distribution surface to the likelihood
bmu <- matrix(0, Q, 1)
bvar <- diag(Q)*.1
plist <- betaSurface(X, y, sigma, bmu = bmu, bvar = bvar)

# likelihood X prior
image(plist$b0, plist$b1, plist$lsurf)
image(plist$b0, plist$b1, plist$psurf + plist$lsurf) # recall log scale
persp(plist$b0, plist$b1, plist$psurf + plist$lsurf) # 3d picture


# analytics: ML
betaHat <- solve(crossprod(X))%*%crossprod(X,y)
betaVar <- solve(crossprod(X))*sigma

# Bayes
V <- solve( crossprod(X)/sigma + solve(bvar) )
v <- crossprod(X, y)/sigma + solve(bvar)%*%bmu
betaMuBayes <- V%*%v
betaVrBayes <- V

# prediction
nsim <- 100
bvalues <- rmvnorm(nsim, betaMuBayes, betaVrBayes)
ypred <- X%*%t(bvalues)
yMu <- rowMeans(ypred)
ySe <- apply(ypred, 1, sd)
yCI <- apply(ypred, 1, quantile, c(.025, .975))
plot(y, yMu, ylim = range(yCI))
segments(y, yCI[1,], y, yCI[2,])


