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

# params
# p = freq of allele A; f = inbreeding coef
# pmake out = c(paa, pab, pbb)
# minf out = minimum amt inbreeding based on allele A freq

pmake(c(.5,0)) #no inbreeding
pmake(c(.5,0.9)) #lots of inbreeding
pmake(c(.5,-0.9)) #negative inbreeding
minf(.5)

# Simulation

# params
npop <- 20               #no. populations
p  <- 0.6                #frequency of allele a
f  <- -.2                #inbreeding coefficient
pr <- pmake(c(p ,f))     #values for (Paa,Pab,Pbb)

# Create data
popsizes <- c(5, 10, 20)        #population sizes
for (n in popsizes) {
  y <- rmultinom(npop,n,pr)    #n random vectors of size m
}

# Gibbs params
g1       <- 1   #beta parameters for p
g2       <- 1
fg       <- f   #initial value of f
priorF   <- 0   #mean and sd for f
priorFSD <- 1
ng     <- 5000
nchain <- 10
nt     <- 500
thin   <- round(seq(1, ng, length=nt))  
pgibbs <- matrix(0,nt,nchain)
fgibbs <- pgibbs

# Beta chains
for(j in 1:nchain){
  
  pg <- rbeta(1,1,1)   #draw initial values from prior
  lg <- minf(pg)  # lower limit for f
  fg <- .tnorm(1,lg,1,priorF,priorFSD) #random draw for f
  k <- 0 #count variable
  
  for(g in 1:ng){
    
    pf <- update_pf()
    pg <- pf$pg
    fg <- pf$fg
    
    if(g %in% thin){
      k <- k + 1
      pgibbs[k,j] <- pg
      fgibbs[k,j] <- fg
    }
  }
}

# Graph of all chains
par(mfrow=c(2,1), mar=c(4,4,1,1),bty='n')
plot(fgibbs[,1],type='l',ylim=c(-1,1), ylab='f')
for(j in 2:nchain)lines(fgibbs[,j])
abline(h=f, lty=2)
plot(pgibbs[,1],type='l',ylim=c(0,1), ylab='p')
for(j in 2:nchain)lines(pgibbs[,j])
abline(h=p, lty=2)

# Graph of joint 
par(mfrow=c(1,1), bty='n')
pseq <- seq(.01,.99,length=100)
plot(pgibbs,fgibbs,xlim=c(0,1),ylim=c(-1,1), cex=.2,xlab='p',
     ylab='f')
lines(pseq,minf(pseq))

.processPars(as.vector(fgibbs[-c(1:100),]),xtrue=f,DPLOT=T)
# red line is true; does pretty well









#eof