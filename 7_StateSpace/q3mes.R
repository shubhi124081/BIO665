## Exercise 3
# Author: Maggie Swift
# Created: 29 March 2020

# total dissolved nitrogen (TDN) and dissolved organic nitrogen (DON) 

#-------------------------------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('../clarkfunctions2020.R')
Rcpp::sourceCpp('../cppFns.cpp')

#-------------------------------------------------------------------------------
# First I read the data and set up a design matrix for the sites and several 
# covariates:

# Data and formula matrix
data <- read.csv('../dataFiles/combinedLevelChem2007.csv',header=T)
formula <- as.formula( ~ prec + Temp*Site + Cl*Site)
tmp <- model.frame(formula, data, na.action=na.pass)
z   <- model.matrix(formula, tmp)

# Find site levels
slevs  <- attr(data$Site,'levels')
sl     <- match(data$Site,slevs)

# Counts
ngroup <- length(slevs)
ns     <- ngroup - 1
nt     <- nrow(data)
q      <- ncol(z)

# X and Y
xnames <- colnames(z)
ynames <- c('TDN','DON')
vnames <- c('sigma', 'tau')
y <- data[,ynames[1]]        # use TDN as response; change to DON later.

#-------------------------------------------------------------------------------
# multiple time series

last  <- which( diff(sl) == 1 )  # first and last row for each site
first <- c(1,last + 1)
last  <- c(last,nt)

#-------------------------------------------------------------------------------
# missing covariates, missing states

v <- getMissX(z, y, first, last)
prows <- v$prows                 # t rows by site
rrows <- v$rrows                 # t+1 rows by site
missz <- v$missx
missy <- v$missy                 # missing y index in stacked data
missY <- v$missY                 # missing y by site as a list

priorx   <- colMeans(z,na.rm=T)  # prior mean for missing covariates
z[missz] <- priorx[missz[,2]]    # initial missing z

#-------------------------------------------------------------------------------
# Initial values for missing states come from a function initialStatesSS.

sigM <- .1                        #prior mean variances
tauM <- .1
wt <- 100
sp <- varPrior(sigM, wt)          #varPrior is gonna spit out a,b for gamma            
tp <- varPrior(tauM, wt/5)        #weak prior

sg <- 1/rgamma(1,sp[1],sp[2])     #initialize variance values
tg <- 1/rgamma(1,tp[1],tp[2])
bg <- matrix(0,q,1)

tmp <- initialStatesSS(y)         #initial missing y
xg  <- tmp$x

#-------------------------------------------------------------------------------
# Gibbs Sampling
ng     <- 2000
vgibbs <- matrix(0, ng, 2)
xgibbs <- matrix(0, ng, nt)
bgibbs <- matrix(0, ng, q)
zgibbs <- matrix(0, ng, length(missz))
colnames(vgibbs) <- vnames
colnames(bgibbs) <- xnames

for(g in 1:ng){
  
  # parameters
  bg <- .updateBeta(z[prows,], xg[rrows] - xg[prows], sigma=sg, 
                    priorIVB=diag(.001,q), priorB=matrix(0,q)) # 
  mu <- z%*%bg  # mean vector
  
  # latent variables by site
  for(j in 1:ngroup){           
    jj <- first[j]:last[j]
    xg[jj] <- updateSSRW(states=xg[jj], yy=y[jj], zb=mu[jj], missing=missY[[j]], 
                         tg=tg, sg=sg)   
  }
  
  # variances
  sg <- updateVariance(xg[rrows],xg[prows] + mu[prows],sp[1],sp[2])
  tg <- updateVariance(y[-missy],xg[-missy],tp[1],tp[2])
  
  # missing observations
  z[missz] <- imputeX( missx = missz, xx = z[prows,], 
                       yy = xg[prows] - xg[rrows],
                       bg = bg, sg, priorx )[missz]
  vgibbs[g,] <- c(sg,tg)
  xgibbs[g,] <- xg
  bgibbs[g,] <- bg
  zgibbs[g,] <- z[missz]
}

#-------------------------------------------------------------------------------
# Latent states plots

keep  <- 200:ng
xpost <- t( apply(xgibbs[keep,],2,quantile,c(.5,.025,.975)) )
ylim  <- range(xpost)

par(mfrow=c(2,2),bty='n', mar=c(3,2,2,1))
for(j in 1:4){
  ij <- first[j]:last[j]
  shadeInterval(ij,xpost[ij,2:3],add=F,ylim=ylim,xlab='Weeks')  
  points(ij,y[ij],lwd=2,cex=.1, col='blue')                 #obs 
  wm <- missY[[j]]
  points(ij[wm],xpost[ij[wm],1],col='orange',cex=.5,pch=16) #missing
  title(slevs[j])
}

.processPars(bgibbs[keep,],xtrue=bg*0,DPLOT=T, burnin = 200) 
.processPars(vgibbs[keep,],xtrue=c(sigM,tauM),DPLOT=T, burnin = 200) #compare prior mean

vars <- c('Temp', 'Cl')
par(mfrow=c(2,1), bty='n', mar=c(4,4,2,2))
sites <- levels(data$Site)
for(k in 1:2){
  
  wk <- grep(vars[k],colnames(bgibbs))
  bk <- bgibbs[keep,wk]
  
  xlim <- range(bk)
  ymax <- 8/diff(xlim)
  plot(NULL, xlim=xlim,ylim=c(0,ymax),xlab=vars[k],ylab='Density')
  
  for(j in 1:ncol(bk)){
    dk <- density(bk[,j])
    lines(dk[[1]],dk[[2]],col=j, lwd=2)
    wmax <- which.max(dk[[2]])
    text( dk[[1]][wmax], 1.2*dk[[2]][wmax], sites[j] ,col=j)
  }
}

# Before running any more analysis, say what you think should be the variance on 
# the process and the observations. Justify your answer.


# How much 'weight' on the priors for $\sigma^2, \tau^2$ is needed to insure that 
# the posterior distribution is close to your prior distribution?
# if I make the weights really low (wt=10), then all of the estimates come out to a 
# gaussian with mean 0. If I make them ridiculously high (wt=1e6), the same also 
# happens. Odd.

# How does the prior on $\sigma^2, \tau^2$ affect estimates of $\mathbf{\beta}$?


# Repeat the analysis with DON as the response variable and interpret results.
#anything to do with site is flipped for DON.
#sigma and tau are underestimated now, whereas before they were overestimated

# What is the answer to the motivating question: Can total dissolved nitrogen 
# (TDN) and dissolved organic nitrogen (DON) be explained by temperature, 
# precipitation, and chlorine concentation ([Cl]), and does it differs between 
# coastal streams?
