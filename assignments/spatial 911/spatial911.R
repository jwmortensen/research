library(LatticeKrig)
library(MASS)
library(FNN)

source("AMCMCUpdate.R")
load("Spatial911PtPtrn.RData")

# Create matrix to plot values
y.grid <- matrix(hrldas.grid[,1],nrow=125)
x.grid <- matrix(hrldas.grid[,2],nrow=125)

# Calculate distance between elements to use in Gaussian process
D <- rdist(pp.grid)

# The values for nu and phi are given in the assignment
nu <- 3.5  
phi <- 343

# Get number of prediction locations and observed locations to use throughout
K <- nrow(pp.grid)

# Break 1428 grid points into B blocks in order to speed up updating
num.blocks <- 51 # choose 51 because it factors evenly into 1428
b.pts <- pp.grid[seq.int(from=1, to=K, length=num.blocks),]

close.pts <- vector("list", length=num.blocks)
close.pts.index <- vector("list", length=num.blocks)
pp.copy <- pp.grid

# This loop allows us to track both the close points and their indices, but requires
# a bit of finagling
for (i in 1:num.blocks) {
  nn <- get.knnx(pp.copy, b.pts[i,], K/num.blocks)
  close.pts[[i]] <- pp.copy[nn$nn.index,]
  close.pts.index[[i]] <- which(paste(pp.grid$Latitude, pp.grid$Longitude) 
    %in% paste(close.pts[[i]]$Latitude, close.pts[[i]]$Longitude))
  pp.copy <- pp.copy[-nn$nn.index,]
}

# Generate Matern correlation matrix
M <- Matern(D, alpha=phi, nu=nu)
M.inv <- solve(M)

# Find nearest prediction location for each observed location
call.locs <- calls[,1:2]
nn <- get.knnx(pp.grid, call.locs, 1)
nn.index <- nn$nn.index
Nk <- rep(0, K)
for (i in nn.index) Nk[i] <- Nk[i] + 1

# Adaptive MCMC Stuff
amcmc <- vector("list", length=num.blocks)
for (i in 1:num.blocks) {
  amcmc[[i]]$mn <- matrix(0, ncol=1, nrow=length(close.pts.index[[i]]))
  amcmc[[i]]$var <- matrix(0, ncol=length(close.pts.index[[i]]), nrow=length(close.pts.index[[i]]))
}
amcmc.it <- 100

# Gibbs sampler to estimate sigma2 and lambda
mh.gibbs <- function(ndraws, sigma.start, sigma2.a, sigma2.b) {
  # initialize containers to hold draws
  lambda.star <- matrix(NA, nrow=ndraws, ncol=K)
  sigma2 <- numeric(ndraws)
  delta <- rgamma(ndraws, shape=nrow(calls)+5, rate=1.2)

  # create initial proposal and initialize variables 
  sigma2[1] <- sigma.start
  lambda.star[1,] <- log((Nk + 1)/(sum(Nk)+1))

  # Initialize functions for use in M-H
  log.like <- function(lambda) {
    sum(Nk*log(lambda))
  }
  
  lambda.prior <- function(lambda.star,sig2){
    -0.5*(t(lambda.star)%*%M.inv%*%lambda.star)/sig2
  }

  for (i in 2:ndraws) {
    # Get draws for sigma2 using the complete conditional
    a <- sigma2.a + K/2
    b <- sigma2.b + (1/2)*t(lambda.star[i-1,])%*%M.inv%*%(lambda.star[i-1,]) 
    sigma2[i] <- 1/rgamma(1, shape=a, rate=b)

    # Here we implement metropolis hastings to get draws for lambda
    lambda.star[i,] <- lambda.star[i-1,]
    for (j in 1:num.blocks) {
      # If enough iterations have been completed, begin using adaptive MCMC techniques
      if (i < amcmc.it) {
        prop.var <- 0.01*diag(K/num.blocks)
      }
      else {
        prop.var <- (2.4^2/(K/num.blocks))*(0.01*diag(K/num.blocks)+amcmc[[j]]$var)
      }

      prop.lstar <- mvrnorm(1, lambda.star[i,close.pts.index[[j]]], prop.var)
      prop.lstar.vec <- lambda.star[i,]
      prop.lstar.vec[close.pts.index[[j]]] <- prop.lstar
      prop.lvec <- exp(prop.lstar.vec)/sum(exp(prop.lstar.vec))
      cur.lvec <- exp(lambda.star[i,])/sum(exp(lambda.star[i,]))
      
      log.MH <- log.like(prop.lvec) - log.like(cur.lvec) 
      log.MH <- log.MH + lambda.prior(prop.lstar.vec,sigma2[i]) - lambda.prior(lambda.star[i,],sigma2[i])
      
      if (log(runif(1)) < log.MH) {
        lambda.star[i,close.pts.index[[j]]] <- prop.lstar
      }
      new.amcmc <- AMCMC.update(lambda.star[i, close.pts.index[[j]]], amcmc[[j]]$mn, amcmc[[j]]$var, i-1)
      amcmc[[j]] <- new.amcmc
    }
  }
  lambda <- exp(lambda.star) / apply(exp(lambda.star), 1, sum)  
  return(list(delta=delta, sigma2=sigma2, lambda=lambda))
}

ndraws <- 10
op.time <- system.time(draws <- mh.gibbs(ndraws, 15, 0.01, 0.01))

lambda.post.mn <- apply(draws$lambda, 2, mean)
plot.grid <- matrix(NA, nrow=125, ncol=125)
plot.grid[kp.gp] <- lambda.post.mn

pdf("callLocationsObs.pdf")
plot(houston)
image.plot(x.grid, y.grid, plot.grid, add=TRUE)
plot(houston, add=TRUE)
dev.off()

pdf("trace.plots.pdf")
plot(1:500, draws$lambda[,1], type="l", xlab="Iteration", main="Small Lambda")
plot(1:500, draws$lambda[,707], type="l", xlab="Iteration", main="Large Lambda")
plot(1:500, draws$lambda[,1231], type="l", xlab="Iteration", main="Medium Lambda")
dev.off()

