library("LatticeKrig")
load("Spatial911PtPtrn.RData")

# Calculate distance between elements to use in Gaussian process
call.locs <- calls[,1:2]
D <- rdist(rbind(pp.grid, call.locs))

# The values for nu and phi are given in the assignment
nu <- 3.5  
phi <- 343

# Get number of prediction locations and observed locations to use throughout
N <- nrow(calls)
K <- nrow(pp.grid)


