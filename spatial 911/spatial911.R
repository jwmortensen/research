library(LatticeKrig)
library(FNN)
load("Spatial911PtPtrn.RData")

# Calculate distance between elements to use in Gaussian process
call.locs <- calls[,1:2]
D <- rdist(pp.grid)

# The values for nu and phi are given in the assignment
nu <- 3.5  
phi <- 343

# Get number of prediction locations and observed locations to use throughout
N <- nrow(calls)
K <- nrow(pp.grid)

# Generate Matern correlation matrix
M <- Matern(D, alpha=phi, nu=nu)
M.inv <- solve(M)
M.chol <- t(chol(M))

# Find nearest prediction location for each observed location
nn <- get.knnx(pp.grid, call.locs, 1)
nn.index <- nn$nn.index

mh.gibbs <- function(ndraws, start.point, mu.mean, mu.cov, sigma2.a, sigma2.b) {
  # initialize containers to hold draws
  lambda <- matrix(NA, nrow=ndraws, ncol=K)
  mu <- numeric(ndraws)
  sigma2 <- numeric(ndraws)
  
  X <- rep(1, length=K)
  mu[1] <- start.point[1]
  sigma2[1] <- start.point[1]
  # How do we initialize lambda? Do we use a random draw? Or is there a deterministic way to calculate from the Gaussian process?
  lambda[1,] <- mu[1] + sqrt(sigma2[1])*M.chol%*%rnorm(K)

  for (i in 2:ndraws) {
    s.star <- solve(t(X)%*%(sigma2[i-1]^(-1)*M.inv)%*%X + 100^(-2))
    m.star <- s.star%*%(t(X)%*%(sigma2[i-1]^(-1)*M.inv)%*%lambda[i-1] + 100^(-2))
    mu[i] <- rnorm(1, mean=m.star, sd=sqrt(s.star))
    
    a <- sigma2.a + K/2
    b <- sigma2.b + (1/2)*t(lambda[i-1]-mu[i])%*%M.inv%*%(lambda[i-1]-mu[i]) 
    sigma2[i] <- 1/rgamma(1, shape=a, rate=b)

    # Here we need to implement metropolis hastings to get draws for lambda
  }
  return(list(mu=mu, sigma2=sigma2, lambda=lambda)
}
