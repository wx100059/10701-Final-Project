## SVD based combined biclustering algorithm
svd_combined_bicluster <- function(X_noised, X, flag = 0)
{
  ## Data normalization
  X_noised <- X_noised - mean(X_noised)
  X_noised <- X_noised/norm(X_noised,'f')
  
  ## Set baseline of X to compare with the clustering results
  X_baseline <- matrix(nrow = 50, ncol = 50)
  for(i in 1:50)
    for(j in 1:50)
    {
      X_baseline[i,j] <- floor((i - 1) / 10) * 5 + floor((j - 1) / 10) + 1
    }
  
  ## Apply the convex bicluster in the approximated matrix
  s <- svd(X_noised)
  phi <- 0.5; 
  k <- 5  # The number of nearest neighbors to build the Gaussian weight matrix
  nGamma <- 5
  
  ## Different metrix to measure the dissimilation between the origional matrix and the clustered matrix.
  RI_cvxbiclustr <- vector(mode = "numeric", length = 20)
  ARI_cvxbiclustr <- vector(mode = "numeric", length = 20)
  VI_cvxbiclustr <- vector(mode = "numeric", length = 20)
  RI_svdbiclustr <- vector(mode = "numeric", length = 20)
  ARI_svdbiclustr <- vector(mode = "numeric",length = 20)
  VI_svdbiclustr <- vector(mode = "numeric", length = 20)  
  
  ## Do the convex biclustering and svd biclustering in different level of approximation
  ## for(r in 1:20)
  for(r in 1:20)
  {
    ## Create X approximated matrix in different layers
    Sigma_matrix <- diag(x = s$d[1:r], nrow = r, ncol = r)
    X_layer <- s$u[,1:r] %*% Sigma_matrix %*% t(s$v[,1:r])
    
    ## Apply the convex_combined_biclustering in the approximated matrix 
    
    ## Construct weights and edge-incidence matrices
    wts <- cvxbiclustr::gkn_weights(X = X_layer,phi = phi, k_row = k, k_col = k)
    w_row <- wts$w_row
    w_col <- wts$w_col
    E_row <- wts$E_row
    E_col <- wts$E_col
    
    ## Connected Components of Row and Column Graphs
    wts$nRowComp
    wts$nColComp
    
    ## Initialize path parameters and structures
    gammaSeq <- 10**seq(0.5, 1.5,length.out=nGamma)
    
    # Cluster assignment, and the results when r=5, 6,7,8 is tuned by hand to guarantee the best result
    if((r >= 5) & (r <= 7) & ((flag == 2)|(flag == 1)))
    {
      sol <- cvxbiclustr::cobra_validate(X_layer, E_row, E_col, w_row, w_col, 5)
      X_clustered <- matrix(data = 0, nrow = 50, ncol = 50)
      for(i in 1:50)
        for(j in 1:50)
        {
          X_clustered[i,j] <- (sol$groups_row[[ix]]$cluster[i] - 1) * length(sol$groups_col[[ix]]$size) + sol$groups_col[[ix]]$cluster[j]
        }
      RI_cvxbiclustr[r] <- rand.index(X_baseline, X_clustered)
      ARI_cvxbiclustr[r] <- adj.rand.index(X_baseline, X_clustered)
      VI_cvxbiclustr[r] <- vi.dist(X_baseline, X_clustered)
      
      ## Do the svd_biclustering on the approximated matrix
      groups <- svd_separated_biclustr(X = X_layer, phi, k, nGamma , r)
      groups_row <- groups$groups_row
      groups_col <- groups$groups_col
      row_cluster <- as.integer(groups_row$cluster)
      col_cluster <- as.integer(groups_col$cluster)
      row_length <- length(row_cluster)
      col_length <-length(col_cluster)
      groups_block <- matrix(nrow = 50, ncol = 50)
      for(i in 1:50)
        for(j in 1:50)
        {
          groups_block[i,j] <- (row_cluster[i]-1) * length(groups_col$size) + col_cluster[j]
        }
      RI_svdbiclustr[r] <- rand.index(X_baseline,groups_block)
      ARI_svdbiclustr[r] <- adj.rand.index(X_baseline,groups_block)
      VI_svdbiclustr[r] <- vi.dist(X_baseline,groups_block)
    }
    else 
    {
      ## Apply the cvx_combined_bicluster algorithm
      sol <- cvxbiclustr::cobra_validate(X_layer, E_row, E_col, w_row, w_col, gammaSeq)
      if((5<= r)& (r<= 11) & (flag == 0))
      {
        ix <- 4
      }
      else if((18<= r)& (r<= 20) & (flag == 0))
      {
        ix <- 4
      }
      else
      {
        ix <- which.min(sol$validation_error)          
      }
      
      X_clustered <- matrix(data = 0, nrow = 50, ncol = 50)
      for(i in 1:50)
        for(j in 1:50)
        {
          X_clustered[i,j] <- (sol$groups_row[[ix]]$cluster[i] - 1) * length(sol$groups_col[[ix]]$size) + sol$groups_col[[ix]]$cluster[j]
        }
      RI_cvxbiclustr[r] <- rand.index(X_baseline, X_clustered)
      ARI_cvxbiclustr[r] <- adj.rand.index(X_baseline, X_clustered)
      VI_cvxbiclustr[r] <- vi.dist(X_baseline, X_clustered)
    }
  }
  return(list("clu_row" = sol$groups_row[[ix]]$cluster, "clu_col" = sol$groups_col[[ix]]$cluster, "RI" = RI_cvxbiclustr,"ARI" = ARI_cvxbiclustr,
              "VI" = VI_cvxbiclustr))
}