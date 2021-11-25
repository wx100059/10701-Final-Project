## try to compare the recover effectiveness of direct convex biclustering and svd based combined biclustering

library(Matrix)
library(igraph)
library(seriation)
library(RMThreshold)
library(fossil)
library(lpSolve)
library(mcclust)
library(cvxclustr)
library(ggplot2)
library(sp)
library(maps)
library(shapefiles)
library(foreign)
library(CVXR)

source("C:/Users/wx100/OneDrive/科研/Master thesis/simulation/svd_comb_biclustr.r")
source("C:/Users/wx100/OneDrive/科研/Master thesis/simulation/svd_sep_biclustr.r")

## Set the iterator time
start_time <- Sys.time()
iterate <- 1;
accuracy_sum1 <- vector(mode = "numeric", length = 30)
accuracy_sum2 <- vector(mode = "numeric", length = 30)
accuracy_sum3 <- vector(mode = "numeric", length = 30)
## Create a matrix of 50 rows and 50 cols with checkboard pattern
X=c(rep(c(rep(1, 10), rep(2, 10), rep(3, 10), rep(4, 10), rep(5, 10)), 10),
    rep(c(rep(2, 10), rep(3, 10), rep(4, 10), rep(5, 10), rep(1, 10)), 10),
    rep(c(rep(3, 10), rep(4, 10), rep(5, 10), rep(1, 10), rep(2, 10)), 10),
    rep(c(rep(4, 10), rep(5, 10), rep(1, 10), rep(2, 10), rep(3, 10)), 10),
    rep(c(rep(5, 10), rep(1, 10), rep(2, 10), rep(3, 10), rep(4, 10)), 10))
X <- matrix(X, nrow = 50, ncol = 50)
type <- c(rep(x = "RI_cvxbiclustr", times = 20),
          rep(x = "ARI_cvxbiclustr", times = 20),
          rep(x = "VI_cvxbiclustr", times = 20),
          rep(x = "RI_svdbiclustr", times = 20),
          rep(x = "ARI_svdbiclustr", times = 20),
          rep(x = "VI_svdbiclustr", times = 20))
rank <- rep(x = 1:20, times = 6)

r <- 2
k <- 5
phi <- 0.5

X_noise <- add.Gaussian.noise(X, stddev = 2.0, symm = FALSE)
s <- svd(X_noise)
Sigma_matrix <- diag(x = s$d[1:r], nrow = r, ncol = r)
X_layer1 <- matrix(nrow = 50, ncol = 50)
X_layer2 <- matrix(nrow = 50, ncol = 50)
X_layer <- list(X_layer1, X_layer2)
sol <- list()

for(r_index in 1:r)
{
  X_layer[[r_index]] <- Sigma_matrix[r_index, r_index] * (s$u[,r_index] %*% t(s$v[,r_index]))
  wts <- cvxbiclustr::gkn_weights(X_layer[[r_index]],phi=phi,k_row=k,k_col=k)
  w_row <- wts$w_row
  w_col <- wts$w_col
  E_row <- wts$E_row
  E_col <- wts$E_col
  ## Connected Components of Row and Column Graphs
  wts$nRowComp
  wts$nColComp
  #### Initialize path parameters and structures
  nGamma <- 5
  gammaSeq <- 10**seq(0,3,length.out=nGamma)
  ## Generate solution path
  sol[[r_index]] <- cvxbiclustr::cobra_validate(X,E_row,E_col,w_row,w_col,gammaSeq)
}

temp_a <- matrix(0, nrow = 50, ncol = 50)
temp_b <- matrix(0, nrow = 50, ncol =50)

for (r_index in 1:r)
{
  temp_b <- temp_b + X_layer[[1]] + X_layer[[2]]
}
wts <- cvxbiclustr::gkn_weights(temp_b,phi=phi,k_row=k,k_col=k)
w_row <- wts$w_row
w_col <- wts$w_col
E_row <- wts$E_row
E_col <- wts$E_col
## Connected Components of Row and Column Graphs
wts$nRowComp
wts$nColComp
#### Initialize path parameters and structures
nGamma <- 5
gammaSeq <- 10**seq(0,3,length.out=nGamma)
## Generate solution path
sol_b <- cvxbiclustr::cobra_validate(temp_b,E_row,E_col,w_row,w_col,gammaSeq)
for (r_index in 1:r) 
{
  temp_a <- temp_a + sol[[r_index]]$U[[1]] 
}