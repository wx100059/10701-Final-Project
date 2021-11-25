#' Detect the cluster pattern in row and column direction
#' @param X the matrix waited to be biclustered
#' @param X_baseline the groundtruth to compute the clustering accuracy
#' @return clu_row: cluster of row vectors
#' @return clu_col: cluster of column vectors
#' @return RI: rand index to measure the clustering accuracy
#' @return ARI: adjusted rand index to measure the clustering accuracy
#' @return VI: variance of information to measure the clustering accuracy

cvx_cluster <- function(X, k = 5, phi = 0.5, gamma)
{
  X <- t(scale(X,center=TRUE,scale=FALSE))
  n <- ncol(X)
  
  ## Pick some weights and a sequence of regularization parameters.
  w <- cvxclustr::kernel_weights(X, phi)
  w <- cvxclustr::knn_weights(w,k,n)
  
  ## Perform convex clustering
  nu <- AMA_step_size(w,n)
  sol <- cvxclust_path_ama(X,w,gamma,nu=nu)
  
  ## Construct adjacency matrix and make hard assignment
  A <- cvxclustr::create_adjacency(sol$V[[1]],w,n)
  sol <- cvxclustr::find_clusters(A)
  
  return(sol$cluster) 
}

#' Use convex bicluster algorithm to detect the row pattern and column pattern
#' @param X the matrix waited to be biclustered
#' @param X_baseline the groundtruth to compute the clustering accuracy
#' @return clu_row: cluster of row vectors
#' @return clu_col: cluster of column vectors
#' @return RI: rand index to measure the clustering accuracy
#' @return ARI: adjusted rand index to measure the clustering accuracy
#' @return VI: variance of information to measure the clustering accuracy
cvx_bicluster <- function(X, X_baseline, k_row = 5, k_col = 5, phi = 0.5, nGamma = 10)
{
  # The number of nearest neighbors to build the Gaussian weight matrix
  
  ## Data normalization
  # X <- X - mean(X)
  # X <- X/norm(X,'f')
    
  ## Construct weights and edge-incidence matrices
  wts <- cvxbiclustr::gkn_weights(X = X, phi = phi, k_row = k_row, k_col = k_col)
  w_row <- wts$w_row
  w_col <- wts$w_col
  E_row <- wts$E_row
  E_col <- wts$E_col
    
  ## Connected Components of Row and Column Graphs
  wts$nRowComp
  wts$nColComp
    
  ## Initialize path parameters and structures
  # gammaSeq <- 10**seq(0, 4,length.out=nGamma)
  gammaSeq <- 4*(10**4)
  ## Apply the convex biclustering algorithm
  sol <- cvxbiclustr::cobra_validate(X, E_row, E_col, w_row, w_col, gammaSeq)
  X_clustered <- matrix(data = 0, nrow = nrow(X), ncol = ncol(X))
  ix <- which.min(sol$validation_error)
  for(i in 1:nrow(X))
    for(j in 1:ncol(X))
      {
        X_clustered[i,j] <- (sol$groups_row[[ix]]$cluster[i] - 1) * length(sol$groups_col[[ix]]$size) + sol$groups_col[[ix]]$cluster[j]
      }
  RI <- rand.index(X_baseline, X_clustered)
  ARI <- adj.rand.index(X_baseline, X_clustered)
  VI <- vi.dist(X_baseline, X_clustered)
  return(list("clu_row" = sol$groups_row[[ix]]$cluster, "clu_col" = sol$groups_col[[ix]]$cluster, "RI" = RI,"ARI" = ARI, "VI" = VI))
}

#' Use convex clustering algorithm to separately detect the row pattern and column pattern
#' @param X the matrix waited to be biclustered
#' @param X_baseline the groundtruth to compute the clustering accuracy
#' @return clu_row: cluster of row vectors
#' @return clu_col: cluster of column vectors
#' @return RI: rand index to measure the clustering accuracy
#' @return ARI: adjusted rand index to measure the clustering accuracy
#' @return VI: variance of information to measure the clustering accuracy

cvx_separate_cluster <- function(X, X_baseline, row_k = 5, col_k = 5, phi = 0.5, row_gamma = 5, col_gamma = 5)
{
  ## Data normalization
  # X <- t(scale(X,center=TRUE,scale=FALSE))
  
  ## Separate cluster in column and row directions
  clu_row <- cvx_cluster(X = t(X), k = row_k, phi = phi, gamma = row_gamma)
  clu_col <- cvx_cluster(X = X, k = col_k, phi = phi, gamma = col_gamma)
  
  m <- nrow(X)
  n <- ncol(X)
  X_clustered <- matrix(data = 0, nrow = m, ncol = n)
  for(i in 1:m)
    for(j in 1:n)
    {
      X_clustered[i,j] <- (clu_row[i] - 1) * col_k + clu_col[j]
    }
  
  ## Compute clustering accuracy index
  RI <- rand.index(X_baseline,X_clustered)
  ARI <- adj.rand.index(X_baseline,X_clustered)
  VI <- vi.dist(X_baseline,X_clustered)
  return(list("clu_row" = clu_row, "clu_col" = clu_col, "RI" = RI, "ARI" = ARI, "VI" = VI)) 
}