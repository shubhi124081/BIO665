path <- "~/Documents/R/BIO665/"
source(paste0(path,'Functions/clarkFunctions2020.r'))
Rcpp::sourceCpp(paste0(path,'Functions/cppFns.cpp')) #? this file doesnt exist.

data  <- read.table(paste0(path,'Data/FACEtrees.txt'),header=T)
form  <- as.formula(cones ~ nfert*trt + diam)
X     <- model.matrix(form, data=data)
Y     <- model.frame(form, data=data)$cones
p     <- ncol(X)
n     <- nrow(X)

bg <- priorB <- matrix(0,p)
xnames <- colnames(X)
rownames(bg) <- xnames
priorVB <- diag(p)*1000
hi <- bg + 100                                # big enough to not matter
lo <- -hi
lo[xnames %in% c('nfert','trt','diam')] <- 0  # positive main effects
cmat    <- .1*solve(crossprod(X))

ng     <- 10000
bgibbs <- matrix(NA,ng,p)  #for regression parameters
colnames(bgibbs) <- colnames(X)
accept <- 0

cupdate <- c(200, 500, 1000, 2000)
accept  <- 0

for(g in 1:ng){
  
  tmp <- updateBetaGLM(bg, cmat, X, Y, likelihood='dpois',
                       link='log', priorB = priorB, priorVB = priorVB,
                       lo = lo, hi = hi)
  bgibbs[g,] <- bg <- tmp$beta
  accept     <- accept + tmp$accept
  
  if(g %in% cupdate){
    cmat <- .1*var(bgibbs[1:(g-1),])  #adapt proposal
    diag(cmat) <- diag(cmat)*1.001
  }
}

.processPars(bgibbs, CPLOT=T)

.processPars(bgibbs[5000:ng,], CPLOT=T)


################## without truncated normal 


#run this part again -----

bg <- priorB <- matrix(0,p)
xnames <- colnames(X)
rownames(bg) <- xnames
priorVB <- diag(p)*1000
hi <- bg + 100                                # big enough to not matter
lo <- -hi
lo[xnames %in% c('nfert','trt','diam')] <- 0  # positive main effects
cmat    <- .1*solve(crossprod(X))

lo <- NULL #set bounds to null
hi <- NULL 

ng     <- 10000
str2 <- matrix(NA,ng,p)  #for regression parameters
colnames(str2) <- paste0(colnames(X), "- Normal")
accept <- 0

cupdate <- c(200, 500, 1000, 2000)
accept  <- 0

for(g in 1:ng){
  
  tmp <- updateBetaGLM(bg, cmat, X, Y, likelihood='dpois',
                       link='log', priorB = priorB, priorVB = priorVB,
                       lo = NULL, hi = NULL)
  str2[g,] <- bg <- tmp$beta
  accept     <- accept + tmp$accept
  
  if(g %in% cupdate){
    cmat <- .1*var(str2[1:(g-1),])  #adapt proposal
    diag(cmat) <- diag(cmat)*1.001
  }
}

.processPars(str2, CPLOT=T)

.processPars(str2[5000:ng,], CPLOT=T)


plotFn <- function(normal = str2, trunNorm = bgibbs, ng = 1e4, burnin = NULL){
  
  par(mfrow = c(3, 2))
  if(is.null(burnin)){
  for(index in 1:ncol(normal)){
plot(normal[, index], type = "l", col = "red", 
     ylim = c(min(min(normal[,index]), min(trunNorm[,index])), 
              max(max(normal[,index]), max(trunNorm[,index]))), main = colnames(trunNorm)[index])
lines(trunNorm[, index], col = "black")}
  }
  
  else{
    for(index in 1:ncol(normal)){
      plot(normal[burnin:ng, index], type = "l", col = "red", 
           ylim = c(min(min(normal[burnin:ng,index]), min(trunNorm[burnin:ng,index])), 
                    max(max(normal[burnin:ng,index]), max(trunNorm[burnin:ng,index]))),
           main = colnames(trunNorm)[index])
      lines(trunNorm[burnin:ng, index], col = "black")}
  }
  
    
    
  }



plotFn() # this is with truncated normal proposal
plotFn(burnin = 5000) #this is with normal proposal 

# While both chains converge, normal and truncated normal, the truncated 
# normal proposal matrix is more efficient. Without the truncated normal, the 
# MCMC tends to get stuck more often, which may cause issues with inference on
# the posterior. 




