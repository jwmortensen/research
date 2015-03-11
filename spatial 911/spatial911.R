library(LatticeKrig)
library(FNN)
load("Spatial911PtPtrn.RData")

# Calculate distance between elements to use in Gaussian process
D <- rdist(pp.grid)

# The values for nu and phi are given in the assignment
nu <- 3.5  
phi <- 343

# Get number of prediction locations and observed locations to use throughout
K <- nrow(pp.grid)

# Generate Matern correlation matrix
M <- Matern(D, alpha=phi, nu=nu)
M.inv <- solve(M)

# Find nearest prediction location for each observed location
call.locs <- calls[,1:2]
nn <- get.knnx(pp.grid, call.locs, 1)
nn.index <- nn$nn.index
Nk <- rep(0, K)
for (i in nn.index) Nk[i] <- Nk[i] + 1

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
    for (j in 1:K) {
      prop.lstar <- rnorm(1, lambda.star[i,j], 0.01)
      prop.lstar.vec <- lambda.star[i,]
      prop.lstar.vec[j] <- prop.lstar
      prop.lvec <- exp(prop.lstar.vec)/sum(exp(prop.lstar.vec))
      cur.lvec <- exp(lambda.star[i,])/sum(exp(lambda.star[i,]))
      
      log.MH <- log.like(prop.lvec) - log.like(cur.lvec) + lambda.prior(prop.lstar.vec,sigma2[i]) - lambda.prior(lambda.star[i,],sigma2[i])
      
      if (log(runif(1)) < log.MH) {
        lambda.star[i,j] <- prop.lstar
      }
    }
  }
  lambda <- exp(lambda.star) / apply(exp(lambda.star), 1, sum)  
  return(list(sigma2=sigma2, lambda=lambda))
}

op.time <- system.time(draws <- mh.gibbs(10, 15, 0.01, 0.01))

quilt.plot(pp.grid$Longitude, pp.grid$Latitude, apply(draws$lambda, 2, mean))
plot(houston, add=TRUE)
points(calls$Longitude, calls$Latitude, col="white", pch=19, cex=0.1)
