# Scratch.R
# Created: 2/25/20
# Author:  Margaret Swift
# Contact: margaret.swift@duke.edu

# Description

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOADING
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('../clarkFunctions2020.R')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCTIONS

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MAIN

# MLE for a single data point
mu   <- 0
y    <- rnorm(1,mu)
mseq <- seq(-4,4,length=100)
like <- dnorm(mseq,y)
plot(mseq,like,type='l')
abline(v=y,lty=2)

# 95% CI
plot(mseq,like,type='l')
abline(v=y,lty=2)
alpha <- .05
ci <- qnorm(c(alpha/2,1-alpha/2),y)
abline(v=ci,lty=3)
lines(mseq,dnorm(mseq,ci[1]),col=2)
lines(mseq,dnorm(mseq,ci[2]),col=2)

# Standard Error
n    <- 10
mu   <- 100
sd   <- 20
yseq <- seq(0,2*mu,by=1)
y    <- rnorm(n,mu,sd)
hist(y,breaks=yseq,probability=T)
se   <- sqrt(var(y)/n)
mle  <- mean(y)
ci   <- c(mle - 1.96*se,mle + 1.96*se)
abline(v=ci,lty=2)

# Standard Error with many repetitions
rep  <- 1000
ci   <- matrix(0,rep,2)
for(i in 1:rep){
  y   <- rnorm(n,mu,sd)
  se  <- sqrt(var(y)/n)
  mle <- mean(y)
  ci[i,] <- c(mle - 1.96*se,mle + 1.96*se)
}
clo <- hist(ci[,1],breaks=yseq,plot=F) #lower conf limit
chi <- hist(ci[,2],breaks=yseq,plot=F) #upper conf limit
hist(y,breaks=yseq,probability=T)
lines(clo$mids,clo$density,type='s')
lines(chi$mids,chi$density,type='s')
cover <- length(which(mu > ci[,1] & mu < ci[,2]))/rep

# Likelihood Ratio
n    <- 10
y    <- rexp(n)                                   #simulated counts
ybar <- mean(y)
rseq <- seq(.3,3,length=100)*ybar
R    <- ybar*( log(rseq) - log(ybar) + 1) - rseq    #normed likelihood

par(mfrow=c(3,1), mar=c(4,4,1,4))
plot(rseq, R, type='l', xlab='', ylab="normed likelihood")
abline(v=ybar,lty=2)
D  <- -2*R                #Deviance
plot(rseq,D,type='l', xlab='', ylab="Deviance") 
abline(v=ybar,lty=2)
abline(h=3.84,lty=2)      #value of D at P = 0.05

P <- 1 - pchisq(D,1)      #P value
plot(rseq, P, type='l', xlab='lambda_0', ylab="P-value") 
abline(v=ybar,lty=2)
abline(h=.05, lty=2)

#find the 95% CI:
CI <- range(rseq[P > .05])
abline(v=CI, lty=2)

# Bootstrapping (cone example)
filename <- ('../dataFiles/FACEtrees.txt')
data <- read.table(filename, header=T)

y <- data[,'cones']
x <- data[,'diam']
w <- is.finite(x) & is.finite(y)
ambient  <- which(data[,'trt'] == 0 & w)
elevated <- which(data[,'trt'] == 1 & w)

Y    <- mean(y[ambient])
X2   <- mean(x[ambient]^2)
nlo  <- length(ambient)
bmuA <- Y/X2
bseA  <- sqrt(Y/nlo)/X2

Y    <- mean(y[elevated])
X2   <- mean(x[elevated]^2)
nhi  <- length(elevated)
bmuE <- Y/X2
bseE  <- sqrt(Y/nhi)/X2

estimates <- signif( matrix( cbind(c(bmuA, bseA),c(bmuE, bseE)), 2, 2), 3 )
rownames(estimates) <- c('mean','se')
colnames(estimates) <- c('ambient','elevated')

bseq <- seq(.002,.007,length=1000)
plot(bseq,dnorm(bseq,estimates['mean',1],estimates['se',1]),type='l',
     xlab='beta',ylab='density')
lines(bseq,dnorm(bseq,estimates['mean',2],estimates['se',1]),col=2)

nboot <- 2000         #no. bootstrap estimates
bvals <- matrix(0, nboot, 2) #matrix to hold estimates
colnames(bvals) <- c('ambient','elevated')

for(b in 1:nboot){
  
  bindex <- sample(ambient, nlo, replace=T) #sample with replacement
  Y    <- mean(y[bindex])
  X2   <- mean(x[bindex]^2)
  bvals[b, 1] <- Y/X2
  
  bindex <- sample(elevated, nhi, replace=T) #sample with replacement
  Y    <- mean(y[bindex])
  X2   <- mean(x[bindex]^2)
  bvals[b, 2] <- Y/X2
}

se <- apply(bvals, 2, sd)

hist(bvals[,1],freq=F,nclass=20,xlim=c(0,.01))
lines(density(bvals[,1]))
ci1 <- quantile(bvals[,1],c(0.025,.975))
abline(v=ci1,lty=2)

elev <- hist(bvals[,2],nclass=20, plot=F)
lines(elev$mids, elev$density, type='s')
lines(density(bvals[,2]))
ci2 <- quantile(bvals[,2],c(0.025,.975))
abline(v=ci2,lty=2)

#EOF