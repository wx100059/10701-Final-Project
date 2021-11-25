#' objective: compute wasserstain distance between two Gaussian distributions
#' @param mu1 the mean vectors of distribution1
#' @param mu2 the mean vectors of distribution2
#' @param cov1 the covariance matrix of distribution1
#' @param cov2 the covariance matrix of distribution2
#' @return Wasserstain distance between two distributions
mean_vector <- function(mu, k = 5, weight = rep(x = 1, times = ncol(mu)))
{
  ## m: the number of row vectors
  ## distance_matrix: stores the Wasserstain distance between each pair of distributions
  m <- nrow(mu)
  n <- ncol(mu)
  distance_matrix <- matrix(data = 0, nrow = m, ncol = m)
  for (i in 2:m)
    for (j in 1:(i-1))
    {
      temp <- 0
      for(p in 1:length(weight))
      {
        temp <- temp + ((mu[i,p] - mu[j,p]) ** 2) * weight[p]
      }
      distance_matrix[i,j] <- sqrt(temp)
      sigma <- 0.25 * mean(distance_matrix)
    }
  
  ## compute the similarity between each pair of distributions
  adjacency_matrix <- exp(- distance_matrix / sigma)
  for (i in 2:m)
    for (j in 1:(i - 1))
    {
      adjacency_matrix[j,i] <- adjacency_matrix[i,j]
    }
  
  cluster <- anocva::spectralClustering(adjacency_matrix, k)
  return(cluster) 
}


#' objective: compute wasserstain distance between two Gaussian distributions
#' @param mu1 the mean vectors of distribution1
#' @param mu2 the mean vectors of distribution2
#' @param cov1 the covariance matrix of distribution1
#' @param cov2 the covariance matrix of distribution2
#' @return Wasserstain distance between two distributions
wasserstain_compute <- function(mu1, mu2, cov1, cov2)
{
  mu_dist = stats::dist(rbind(mu1, mu2))^2
  cov1_half_power <- Matpow(cov2, 0.5)
  cov2_half_power <- Matpow(cov2, 0.5)
  cov_dist <- sum(diag(cov1 + cov2 - 2 * Matpow(cov1_half_power %*% cov2 %*% cov1_half_power, 0.5)))
  return(mu_dist + cov_dist)
}


#' objective: Use Wasserstain distance based spectral clustering to cluster distributions.
#' @param mu a m x n matrix where each row is a n-dimensional mean vector of a distribution
#' @param cov a list of m n x n matrices where each matrix is a covariance matrix 
#' of a distribution
#' @param k number of clusters
#' @return clusters: a group of clusters c_1, c_2,...c_k where c_i is the ith cluster of distributions
#' @return RI: rand index to measure the clustering accuracy
#' @return ARI: adjusted rand index to measure the clustering accuracy
#' @return VI: variance of information to measure the clustering accuracy
wasserstain_cluster <- function(mu, cov, k = 5)
{
  ## m: the number of row vectors
  ## distance_matrix: stores the Wasserstain distance between each pair of distributions
  m <- nrow(mu)
  n <- ncol(mu)
  distance_matrix <- matrix(data = 0, nrow = m, ncol = m)
  for (i in 2:m)
    for (j in 1:(i-1))
    {
      distance_matrix[i,j] <- wasserstain_compute(mu[i, ], mu[j, ], cov[[i]], cov[[j]])
      sigma <- 0.25 * mean(distance_matrix)
    }
  
  ## compute the similarity between each pair of distributions
  adjacency_matrix <- exp(- distance_matrix / sigma)
  for (i in 2:m)
    for (j in 1:(i - 1))
    {
      adjacency_matrix[j,i] <- adjacency_matrix[i,j]
    }
  
  clusters <- anocva::spectralClustering(adjacency_matrix, k)
  return(cluster) 
}

#' Use Bhattacharyya distance based spectral clustering to cluster distributions.
#' @param mu: a m x n matrix where each row is a n-dimensional mean vector of a distribution
#' @param cov a list of m n x n matrices where each matrix is a covariance matrix 
#' of a distribution
#' @param k number of clusters
#' @return clusters: a group of clusters c_1, c_2,...c_k where c_i is the ith cluster of distributions
#' @return RI: rand index to measure the clustering accuracy
#' @return ARI: adjusted rand index to measure the clustering accuracy
#' @return VI: variance of information to measure the clustering accuracy
bhattacharyya_cluster <- function(mu, cov, k = 5)
{
  ## k: number of clusters
  ## m: the number of row vectors
  ## distance_matrix: stores the Bhattacharyya distance between each pair of distributions
  k <- 5 # number of clusters
  m <- nrow(mu)
  n <- ncol(mu)
  
  ## Compute the Bhattacharyya distance matrix
  distance_matrix <- matrix(data = 0, nrow = m, ncol = m)
  for(i in 2:m)
    for(j in 1:(i - 1))
    {
      distance_matrix[i,j] <- fpc::bhattacharyya.dist(mu[i, ], mu[j, ], cov[[i]], cov[[j]])
      sigma <- 0.5 * mean(distance_matrix)
    }
  
  ## compute the similarity between each pair of distributions
  adjacency_matrix <- exp(- distance_matrix / sigma)
  for (i in 2:m)
    for (j in 1:(i - 1))
    {
      adjacency_matrix[j,i] <- adjacency_matrix[i,j]
    }
  
  clusters <- anocva::spectralClustering(adjacency_matrix, k)
  return(clusters) 
}

#' Achieve clustering in one dimension based on the clustering information of another dimension. 
## Each row can be seen as a k dimensional distribution.
#' @param X a m x n noised matrix with additive Gaussian noise. Elements of a row vectors belongs to k clusters.
#' The elements in the same cluster can be view as a i.i.d sample of a hidden distribution.
#' These k distributions are independent
#' @param partition the clustering information of columns
#' @param  k number of clusters in row/column direction
#' @param type choose the corresponded multiple sample clustering algorithm,'wasserstain distance',
#'  'bhattacharyya distance', and 'mean vector'
#' @param flag can be "row" or "col", achieve row or column clustering based on the value of flag
#'
#' @return cluster: the cluster for column/row vectors
spectral_multiple_cluster <- function(X, partition, k = 5, flag = 'row', type = "BH")
{
  ## store the indexes of vectors belonging to the same cluster
  cluster_index <- list()  
  for(i in 1:max(partition))
  {
    cluster_index[[i]] <- which(partition == i)
  }
  weight <- rep(x = 0, times = length(cluster_index))
  for(i in 1:length(cluster_index))
  {
    weight[i] <- length(cluster_index[[i]])
  }
  # ## Generate a low rand rebuilt matirx X_layer of X
  # s <- svd(X)
  # Sigma_matrix <- diag(x = s$d[1:r], nrow = r, ncol = r)
  # X_layer <- s$u[ ,1:r] %*% Sigma_matrix %*% t(s$v[ ,1:r])
  
  
  ## n is the number of clusters in row direction
  ## mu is a m x n matrix where each row is a n-dimensional mean vector of a distribution
  ## cov is a list of m n x n matrices where each matrix is a covariance matrix 
  if(flag == "row")
  {
    m <- nrow(X)
    n <- length(cluster_index)
    mu <- matrix(data = NA, nrow = m, ncol = n)
    cov <- list()
    cov_temp <- matrix(data = 0, nrow = n, ncol = n)
    for(i in 1:m)
    {
      for(j in 1:n)
      {
        ## Compute the length weighted mean vector
        mu[i,j] <- mean(X[i, cluster_index[[j]]])
        if(length(cluster_index[[j]]) == 1)
        {
          cov_temp[j,j] <- 0.1
        }
        else
        {
          cov_temp[j,j] <- var(X[i, cluster_index[[j]]])          
        }

      }
      cov[[i]] <- cov_temp
    }
    if(type == 'MEAN')
    {
      return(mean_vector(mu, k, weight))
    }
    else if(type == "WA")
    {
      return(wasserstain_cluster(mu, cov, k))
    }
    
    else if(type == "BH")
    {
      return(bhattacharyya_cluster(mu, cov, k))
    }
    
    else
    {
      print("type undefined!")
      exit()
    }
  }

  else if (flag == "col")
  {
    m <- length(cluster_index)
    n <- ncol(X)
    mu <- matrix(data = 0, nrow = n, ncol = m)
    cov <- list()
    cov_temp <- matrix(data = 0, nrow = m, ncol = m)
    for(i in 1:n)
    {
      for(j in 1:m)
      {
        mu[i,j] <- mean(X[cluster_index[[j]], i])
        if(length(cluster_index[[j]]) == 1)
        {
          cov_temp[j,j] <- 0.1
        }
        else
        {
          cov_temp[j,j] <- var(X[cluster_index[[j]], i])
        }
      }
      cov[[i]] <- cov_temp
    }
    
    if(type == 'MEAN')
    {
      return(mean_vector(mu, k, weight))
    }
    else if(type == "WA")
    {
      return(wasserstain_cluster(mu, cov, k))
    }
    
    else if(type == "BH")
    {
      return(bhattacharyya_cluster(mu, cov, k))
    }
    
    else
    {
      print("type undefined!")
      exit()
    }
  }
  else
  {
    print("flag undefined!")
    exit()
  }
  
}



#' spectral clustering for the column vectors
#' @param X noised matrix
#' @param k the number of clusters
#' @return cluster: the cluster assignment for column vectors
col_spectral <- function(X, k = 5)
{
  ## Create the distance matrix
  n <- ncol(X)
  distance_matrix <- matrix(data = 0, nrow = n, ncol = n)
  for (i in 2:n)
    for (j in 1:(i-1))
    {
      distance_matrix[i,j] = pracma::Norm(X[ ,i] - X[ ,j])
    }
  sigma <- 2 * mean(distance_matrix)
  
  ## Create the adjacency matrix
  adjacency_matrix <- exp(- distance_matrix / sigma)
  for (i in 2:n)
    for (j in 1:(i - 1))
    {
      adjacency_matrix[j,i] <- adjacency_matrix[i,j]
    }
  
  ## Apply the convex clustering algorithm on this
  cluster <- anocva::spectralClustering(adjacency_matrix, k)
  return(cluster) 
}

#' spectral clustering for the row vectors
#' @param X noised matrix
#' @param k the number of clusters
#' @return cluster: the cluster assignment for row vectors
row_spectral <- function(X, k = 5)
{
  ## Create the distance matrix
  m <- nrow(X)
  distance_matrix <- matrix(data = 0, nrow = m, ncol = m)
  for (i in 2:m)
    for (j in 1:(i-1))
    {
      distance_matrix[i,j] = pracma::Norm(X[i, ] - X[j, ])
    }
  sigma <- 2 * mean(distance_matrix)
  
  ## Create the adjacency matrix
  adjacency_matrix <- exp(- distance_matrix / sigma)
  for (i in 2:m)
    for (j in 1:(i - 1))
    {
      adjacency_matrix[j,i] <- adjacency_matrix[i,j]
    }
  
  ## Apply the convex clustering algorithm on this
  cluster <- anocva::spectralClustering(adjacency_matrix, k)
  return(cluster) 
}

#' Separately apply spectral clustering on the column and row vectors to detect the row and column partitions
#' @param X noised matrix
#' @param X_baseline the ground truth to compute the clustering accuracy
#' @param row_k the number of clusters in row direction
#' @param col_k the number of clusters in column direction
#' @return clu_row: cluster of row vectors
#' @return clu_col: cluster of column vectors
#' @return RI: rand index to measure the clustering accuracy
#' @return ARI: adjusted rand index to measure the clustering accuracy
#' @return VI: variance of information to measure the clustering accuracy
spectral_separate_bicluster <- function(X, X_baseline, row_k = 5, col_k = 5)
{
  ## Apply the spectral clustering in column and row direction separately
  col_cluster <- col_spectral(X = X, k = col_k)
  row_cluster <- row_spectral(X = X, k = row_k)
  
  ## Compute the clustering accuracy
  m <- nrow(X)
  n <- ncol(X)
  X_clustered <- matrix(data = 0, nrow = m, ncol = n)
  for(i in 1:m)
    for(j in 1:n)
    {
      X_clustered[i,j] <- (row_cluster[i] - 1) * col_k + col_cluster[j]
    }
  RI <- fossil::rand.index(X_baseline, X_clustered)
  ARI <- fossil::adj.rand.index(X_baseline, X_clustered)
  VI <- mcclust::vi.dist(X_baseline, X_clustered)
  return(list("row_cluster" = row_cluster, "col_cluster" = col_cluster, "RI" = RI, "ARI" = ARI, "VI" = VI))
}

#' Apply multiple spectral clustering to detect the row and column partitions
#' @param X noised matrix
#' @param X_baseline the ground truth to compute the clustering accuracy
#' @param row_k the number of clusters in row direction
#' @param col_k the number of clusters in column direction
#' @return clu_row: cluster of row vectors
#' @return clu_col: cluster of column vectors
#' @return RI: rand index to measure the clustering accuracy
#' @return ARI: adjusted rand index to measure the clustering accuracy
#' @return VI: variance of information to measure the clustering accuracy
spectral_multiple_bicluster <- function(X, X_baseline, row_k = 5, col_k = 5, type = "BH")
{
  ## Compute the initial clusters
  col_cluster <- col_spectral(X = X, k = col_k)
  col_cluster_record <- seq(from = 0, to = 0, length.out = ncol(X))
  row_cluster <- spectral_multiple_cluster(X = X, partition = col_cluster, k = row_k, type = type)
  row_cluster_record <- seq(from = 0, to = 0, length.out = nrow(X))
  
  ## Iterate to compute the row cluster and column clusters until converge
  if ((typeof(all.equal(col_cluster_record, col_cluster)) == "character") || (typeof(all.equal(row_cluster_record, row_cluster)) == "character"))
  # if (isTURE(all.equal(col_cluster_record, col_cluster)) || isTURE(all.equal(row_cluster_record, row_cluster)))
  {
    col_cluster_record <- col_cluster
    row_cluster_record <- row_cluster
    col_cluster <- spectral_multiple_cluster(X = X, partition = row_cluster, flag = "col", type = type)
    row_cluster <- spectral_multiple_cluster(X = X, partition = col_cluster, flag = "row", type = type)
  }
  
  ## Compute the clustering accuracy
  m <- nrow(X)
  n <- ncol(X)
  X_clustered <- matrix(data = 0, nrow = m, ncol = n)
  for(i in 1:m)
    for(j in 1:n)
    {
      X_clustered[i,j] <- (row_cluster[i] - 1) * col_k + col_cluster[j]
    }
  RI <- fossil::rand.index(X_baseline, X_clustered)
  ARI <- fossil::adj.rand.index(X_baseline, X_clustered)
  VI <- mcclust::vi.dist(X_baseline, X_clustered)
  return(list("row_cluster" = row_cluster, "col_cluster" = col_cluster, "RI" = RI, "ARI" = ARI, "VI" = VI))
}