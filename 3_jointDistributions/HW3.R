# HW3.R
# Created: 2/6/2020
# Author:  Margaret Swift
# Contact: margaret.swift@duke.edu

# Homework 3 exercises

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOADING
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('../clarkFunctions2020.R')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCTIONS

# none

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Exercise 1
n <- 100
pi0 <- 0.9
pi1 <- 0.5
t <- .9

infect0 <- ((1 - pi1) * t) / ( (1-pi0)*(1-t) + (1-pi1)*t )
infect1 <- (pi1 * t) / ( pi0*(1-t) + pi1*t )

plot(dbinom(50:150, n, infect0), type='l', col='red')
lines(dbinom(50:150, n, infect1), col='purple')
# not sure that this is right??

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Exercise 2

m  <- 50                      # no. at risk
S  <- 0:m          
pi <- .35                     # survival Pr
b  <- 4:2                     # beta parameter

# Initialize plot
a  <- signif(b[1]/(1/pi - 1),3)  # 2nd beta parameter to give mean value = pi
plot(S/m, dbeta(S/m, a, b[1]), xlab=expression( pi ), 
     ylab=expression( paste("[", pi, "]") ), type='l')

title('beta density for a, b')
cols <- c('black', 'red', 'purple')
locs <- c(1.75, 1.5, 1.25)
for (i in 1:length(b)) {
  a  <- signif(b[i]/(1/pi - 1),3)  # 2nd beta parameter to give mean value = pi
  lines(S/m, dbeta(S/m, a, b[i]), col=cols[i])
  ptext <- paste( "(", a, ", ", b[i], ")", sep="")
  text(1,locs[i], ptext,  pos=2, col=cols[i])
}

# comparing variances?
meanS <- sum( S*dbetaBinom(S, m, mu=pi,b=b) )
varS  <- sum( (S - meanS)^2*dbetaBinom(S, m, mu=pi,b=b) )

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Exercise 3








#eof