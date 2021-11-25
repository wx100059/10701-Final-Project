#' Use svd to separate the clustering information and separately detect the clustering information
#'
#' Input: 
#' @param X m x n matrix wait to be clustered
#' @param X_baseline the ground truth of matrix X
#' @param phi weight that controls the penalty term of checkerboard pattern
#' @param k hyperparameter for k-nearest neighborhood
#' @param nGamma how many regularization parameters are tried
#' @param r control the rank of singular matrix
#' @param method Algorithm to use: "admm" or "ama"
#' @param flag if auto-tuning parameter. 0 -> off, 1 -> on
#'
# Output: 
#' @param cluster cluster of row vectors
#' @param RI rand index to measure the clustering accuracy
#' @param ARI adjusted rand index to measure the clustering accuracy
#' @param VI variance of information to measure the clustering accuracy
svd_separated_bicluster<-function(X, X_baseline, phi = 0.5,k = 5, nGamma = 5, r = 3, method = 'admm', flag = 0)
{
  m <- dim(X)[1]
  n <- dim(X)[2]
  
  ## Use svd to seperate the row and column information 
  s <- svd(X)
  Sigma <- diag(sqrt(s$d[1:r]), r, r)
  UT <- Sigma %*% (t(s$u)[1:r, ])
  VT <- Sigma %*% (t(s$v)[1:r, ])
  
  # for r equals the specific value, gamma parameter need to be generated in different way
  if((r >= 9) & (r <= 14) & (flag == 1))
  {
    ix <- 1
    gammaSeq = 1000  
    weight_of_U <- cvxclustr::kernel_weights(UT, phi)
    weight_of_U <- cvxclustr::knn_weights(weight_of_U, k, m)
    # weight_of_V <- cvxclustr::kernel_weights(VT, phi)
    # weight_of_V <- cvxclustr::knn_weights(weight_of_V, k, n)
    sol_U <- cvxclust(UT, weight_of_U, gammaSeq, method = 'admm') ## the gamma parameter need to be explored later
    # sol_V <- cvxclust(VT, weight_of_V, gammaSeq, method = 'admm')
    
    U <- t(sol_U$U[[1]])
    # V <- t(sol_V$U[[1]])
    recovered <- U %*% t(V)
    
    ## Cluster assignment
    A_row <- cvxclustr::create_adjacency(sol_U$V[[ix]],weight_of_U[[ix]],m,method = 'admm')
    clu_row <- cvxclustr::find_clusters(A_row)
    
    # A_col <- cvxclustr::create_adjacency(sol_V$V[[ix]],weight_of_V[[ix]], n, method = 'admm')
    # clu_col <- cvxclustr::find_clusters(A_col)
    return(list(groups_row=clu_row, groups_col=clu_col, A_row=A_row, A_col=A_col)) 
  }
  else
  {
    gammaSeq = 10**seq(0, 4,length.out=nGamma)  
    weight_of_U <- cvxclustr::kernel_weights(UT, phi)
    weight_of_U <- cvxclustr::knn_weights(weight_of_U, k, m)
    # weight_of_V <- cvxclustr::kernel_weights(VT, phi)
    # weight_of_V <- cvxclustr::knn_weights(weight_of_V, k, n)
    ix <- 1
    sol_U <- cvxclust(UT, weight_of_U, gammaSeq[ix], method = 'admm') ## the gamma parameter need to be explored later
    # sol_V <- cvxclust(VT, weight_of_V, gammaSeq[ix], method = 'admm')
    
    U <- t(sol_U$U[[1]])
    # V <- t(sol_V$U[[1]])
    # recovered <- U %*% t(V)
    
    ## Cluster assignment
    A_row <- cvxclustr::create_adjacency(sol_U$V[[ix]], weight_of_U[[ix]], m, method = 'admm')
    clu_row <- cvxclustr::find_clusters(A_row)
    # A_col <- cvxclustr::create_adjacency(sol_V$V[[ix]],weight_of_V[[ix]], n, method = 'admm')
    # clu_col <- cvxclustr::find_clusters(A_col)
    
    ## Compute the clustering accuracy
    row_cluster <- as.integer(clu_row$cluster)
    # col_cluster <- as.integer(clu_col)
    # row_length <- length(row_cluster)
    # col_length <-length(col_cluster)
    # groups_block <- matrix(nrow = 50, ncol = 50)
    # for(i in 1:50)
    #   for(j in 1:50)
    #   {
    #     groups_block[i,j] <- (row_cluster[i]-1) * length(groups_col$size) + col_cluster[j]
    #   }
    # RI[r] <- rand.index(X_baseline,groups_block)
    # ARI[r] <- adj.rand.index(X_baseline,groups_block)
    # VI[r] <- vi.dist(X_baseline,groups_block)
    RI <- rand.index(X_baseline,row_cluster)
    ARI <- adj.rand.index(X_baseline,row_cluster)
    VI <- vi.dist(X_baseline,row_cluster)
    return(list("cluster" = clu_row$cluster, "RI" = RI, "ARI" = ARI, "VI" = VI)) 
  }
}