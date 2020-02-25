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

mu   <- 0
y    <- rnorm(1,mu)
mseq <- seq(-4,4,length=100)
like <- dnorm(mseq,y)
plot(mseq,like,type='l')
abline(v=y,lty=2)









#EOF