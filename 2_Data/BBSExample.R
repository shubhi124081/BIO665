# BBSExample.R
# Margaret Swift
# 1/15/2020
# margaret.swift@duke.edu

# BBS Example

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOADING
library(rjags)
library(coda)

setwd("/Users/margaretswift/Documents/2020_Spring/BIO665/2_Data")
source('../clarkFunctions2020.R')
load('BBSexampleTime.rdata')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCTIONS


makeSpace <- function(){
  outputFormat = knitr::opts_knit$get("rmarkdown.pandoc.to")
  if(outputFormat == 'latex')
    "\\bigskip"
  else if(outputFormat == 'html')
    "<br>"
  else
    "<br>"
}
colFmt = function(x,color){
  outputFormat = knitr::opts_knit$get("rmarkdown.pandoc.to")
  if(outputFormat == 'latex')
    paste("\\textcolor{",color,"}{",x,"}",sep="")
  else if(outputFormat == 'html')
    paste("<font color='",color,"'>",x,"</font>",sep="")
  else
    x
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MAIN

# Put the points on the map to see the routes
maps::map('county', xlim = c(-85, -75), ylim = c(33.6, 36.8), col='grey')
maps::map('state', xlim = c(-85, -75), ylim = c(33.6, 36.8), add=T)
points(xdata$lon, xdata$lat, pch = 16, cex = 1)

# Select species
spec <- 'WoodThrush'
y <- ydata[,spec]

# Assign bins to each deficit value
def <- xdata$defSite # moisture deficit
def <- (def - mean(def, na.rm=T))/sd(def, na.rm=T)
df <- seq(-min(def), max(def), length=20)
di <- findInterval(def, df, all.inside = T)

# Set colors (moisture deficit) and sizes (#observed)
cex <- 5*ydata[,spec]/max(y, na.rm=T)
colT <- colorRampPalette( c('#8c510a','#d8b365','#01665e','#2166ac') )
cols <- rev( colT(20) )

# Re-draw map with new information
maps::map('county', xlim = c(-85, -75), ylim = c(33.6, 36.8), col='grey')
maps::map('state', xlim = c(-85, -75), ylim = c(33.6, 36.8), add=T)
points(xdata$lon, xdata$lat, pch = 16, cex = cex, col = cols[di] )
title(spec)

#--------------------------------------------------------
# Trends
# route by year

plot(xdata$year, y, col=cols[10])
spec.t <- tapply(y, list( route=xdata$Route, year = xdata$year),
                               mean, na.rm=T)
# average by year
spec.year <- colMeans( spec.t, na.rm=T)
year <- as.numeric(names(spec.year))
spec.quant <- apply( spec.t, 2, quantile, na.rm=T)

# Draw lines for quantiles
lines(year, spec.year, lwd=5, col='white')
lines(year, spec.year, lwd=2)
lines(year, spec.quant[2,], lwd=5, col='white')
lines(year, spec.quant[2,], lty=2, lwd=2)
lines(year, spec.quant[4,], lwd=5, col='white')
lines(year, spec.quant[4,], lty=2, lwd=2)
title(spec)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# JAGS

# Save wood thrush model to text file
file <- "jagsWoodThrush.txt"                # model file
cat("model{
  for(i in 1:n){                             
    Y[i] ~ dpois(lambda[i])                  # stochastic Y
    lambda[i] <- exp( inprod(beta[],X[i,]) ) # lambda deterministic
  }
   for (i in 1:p) {
     beta[i] ~ dnorm(0, 1.0E-5)              
   }
}", file = file)

# Create a factor version for each variable
xdata$routeFactor     <- as.factor(xdata$Route)
xdata$StartWindFactor <- as.factor(xdata$StartWind)
xdata$SkyFactor       <- as.factor(xdata$StartSky)

# Fill in X and Y
xmiss <- which(is.na(xdata),arr.ind=T)
X <- model.matrix(y ~ temp + year + StartWindFactor, data=xdata)
Y <- model.frame(y ~ temp + year + StartWindFactor, data=xdata)$y
dataList <- list(Y = Y, X = X, n = nrow(X), p = ncol(X))
save(dataList,file='jagsData.Rdata')

# computation
outjags <- jags.model(data=dataList, inits=NULL, file=file, 
                      n.chains=3, n.adapt=500 )






par(mfrow=c(1,2), bty='n', mar=c(4,3,1,1))
hist(xdata$year, main='', nclass=60)   # routes by year
title('effort by year')

# observations by route

xdata$Route <- as.numeric( as.character(xdata$Route) )

n_route <- table( xdata$Route )
mm      <- match( as.numeric(names(n_route)), xdata$Route)
routeSummary <- data.frame(n_route, xdata[mm, c('lon','lat')],
                           stringsAsFactors = F)
colnames(routeSummary)[1:2] <- c('route','nobs')

# symbols scaled by observation times per route
map('county',region='north carolina')
points(routeSummary$lon,routeSummary$lat, cex=.05*routeSummary$nobs, pch=16, col='grey')
points(xdata$lon,xdata$lat, cex=.3, pch=16)
title('effort by route')


firstYr <- tapply(xdata$year, xdata$Route, min, na.rm=T)
firstYr <- firstYr[ as.character(routeSummary$route) ]
routeSummary <- cbind( routeSummary, firstYr )

col <- rep('blue', nrow(routeSummary))
col[ routeSummary$firstYr > 1987 ] <- 'green'

map('county',region='north carolina')
points(routeSummary$lon,routeSummary$lat, cex=.05*routeSummary$nobs, pch=16, col=col)



plot(xdata$year, xdata$juneTemp, cex=.3, ylab='degrees C')
title('temperature trend')

routes <- sort(unique(routeSummary$route))
nroute <- length(routes)

for(i in 1:nroute){
  
  wi <- which(xdata$Route == routes[i])
  lines(xdata$year[wi], xdata$juneTemp[wi], col = col[i])
}
tmu <- tapply(xdata$juneTemp, xdata$year, mean, na.rm=T)
years <- as.numeric(names(tmu))
lines( years, tmu, col=2, lwd = 3)


trend <- lm(tmu ~ years)

# local anomaly
tmuByRoute <- tapply(xdata$juneTemp, xdata$Route, mean, na.rm=T)
tanomaly   <- xdata$juneTemp - tmuByRoute[ as.character(xdata$Route) ]

plot(xdata$year, tanomaly)
abline(h = 0)
# EOF