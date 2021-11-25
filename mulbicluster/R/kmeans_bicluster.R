#' Apply k means clustering based separated biclustering algorithm to detect the row and column partitions
#' @param X noised matrix
#' @param X_baseline the ground truth to compute the clustering accuracy
#' @param row_k the number of clusters in row direction
#' @param col_k the number of clusters in column direction
#' @return clu_row: cluster of row vectors
#' @return clu_col: cluster of column vectors
#' @return RI: rand index to measure the clustering accuracy
#' @return ARI: adjusted rand index to measure the clustering accuracy
#' @return VI: variance of information to measure the clustering accuracy
kmeans_separate_bicluster <- function(X, X_baseline, row_k = 5, col_k = 5)
{
  ## Apply the spectral clustering in column and row direction separately
  col_cluster <- stats::kmeans(x = t(X), centers = col_k)
  col_cluster <- col_cluster$cluster
  row_cluster <- stats::kmeans(x = X, centers = row_k)
  row_cluster <- row_cluster$cluster
  
  ## Compute the clustering accuracy
  m <- nrow(X)
  n <- ncol(X)
  X_clustered <- matrix(data = 0, nrow = m, ncol = n)
  for(i in 1:m)
    for(j in 1:n)
    {
      X_clustered[i,j] <- (row_cluster[i] - 1) * col_k + col_cluster[j]
    }
  RI <- rand.index(X_baseline, X_clustered)
  ARI <- adj.rand.index(X_baseline, X_clustered)
  VI <- vi.dist(X_baseline, X_clustered)
  return(list("row_cluster" = row_cluster, "col_cluster" = col_cluster, "RI" = RI, "ARI" = ARI, "VI" = VI))
}

#' objective: Use KL divergence based k means clustering to cluster distributions.
#' @param mu a m x n matrix where each row is a n-dimensional mean vector of a distribution
#' @param cov a list of m n x n matrices where each matrix is a covariance matrix 
#' of a distribution
#' @param k number of clusters
#' @return clusters: a group of clusters c_1, c_2,...c_k where c_i is the ith cluster of distributions
#' @return RI: rand index to measure the clustering accuracy
#' @return ARI: adjusted rand index to measure the clustering accuracy
#' @return VI: variance of information to measure the clustering accuracy
kl_cluster <- function(mu, cov, k = 5, max_it = 1000)
{
  m <- nrow(mu)
  n <- ncol(mu)
  loop_time <- 0
  ## Choose the initial clustering center randomly
  center <- sample(x = 1:m, size = k, replace = FALSE)
  mean_center <- mu[center, ]
  cov_center <- cov[center]
  
  ## assignment stores the cluster that each distribution belongs to
  assignment <- rep(x = 0, times = m)
  for (i in 1:m)
  {
    kl_divergence <- rep(x = 0, times = k)
    for (j in 1:k)
    {
      kl_divergence[j] <- gaussDiff::normdiff(mu1 = mu[i, ], sigma1 = cov[[i]], mu2 = mean_center[j, ], sigma2 = cov_center[[j]], method = "KL")
    }
    assignment[i] <- which.min(kl_divergence)    
  }
  
  ## assignment_record is used to verify if k means converge
  assignment_record <- rep(x = 0, times = m)
  while(typeof(all.equal(assignment_record, assignment)) == "character")
  {
    assignment_record <- assignment
    loop_time <- loop_time + 1
    ## Update cluster means and covarance matrices
    for(j in 1:k)
    {
      ## cluster_j is the index of vectors belonging to the same cluster
      cluster_j <- which(j == assignment)
      
      ## If no elements in this cluster, then randomly choose an element from another cluster
      if(length(cluster_j) == 0)
      {
        cluster_j <- sample(x = 1:m, size = 1)       
      }

      ## Update the mean of the clustering center of these distribution
      mean_center[j, ] <- apply(X = matrix(mu[cluster_j, ], ncol = n), MARGIN = 2, FUN = mean)
      cov_center[[j]] <- matrix(data = 0, nrow = n, ncol = n)
      for (i in cluster_j)
      {
        cov_center[[j]] <- cov_center[[j]] + cov[[i]] + (mu[i, ] - mean_center[j]) %*% t((mu[i, ] - mean_center[j]))
      }
      cov_center[[j]] <- cov_center[[j]] / length(cluster_j)
    }
    
    ## Update the cluster assignment
    for (i in 1:m)
    {
      kl_divergence <- rep(x = 0, times = k)
      for (j in 1:k)
      {
        kl_divergence[j] <- gaussDiff::normdiff(mu1 = mu[i, ], sigma1 = cov[[i]], mu2 = mean_center[j, ], sigma2 = cov_center[[j]], method = "KL")
      }
      assignment[i] <- which.min(kl_divergence)    
    }
    if(loop_time >= max_it)
    {
      print("too many iterations, restart kl_cluster function!")
      kl_cluster(mu = mu, cov = cov, k = k)
    }
  }
  return(assignment)
}

#' Achieve k means multiple sample clustering in one dimension based on the clustering information of another dimension. 
## Each row can be seen as a k dimensional distribution.
#' @param X a m x n noised matrix with additive Gaussian noise. Elements of a row vectors belongs to k clusters.
#' The elements in the same cluster can be view as a i.i.d sample of a hidden distribution.
#' These k distributions are independent
#' @param partition the clustering information of columns
#' @param  k number of clusters in row/column direction
#' @param type choose the corresponded multiple sample clustering algorithm
#' @param flag can be "row" or "col", achieve row or column clustering based on the value of flag
#'
#' @return cluster: the cluster for column/row vectors
kmeans_multiple_cluster <- function(X, partition, k = 5, flag = 'row', type = 'KL')
{
  ## store the indexes of vectors belonging to the same cluster
  cluster_index <- list()  
  for(i in 1:max(partition))
  {
    cluster_index[[i]] <- which(partition == i)
  }
  
  
  if(type == 'KL')
  {
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
          mu[i,j] <- mean(X[i, cluster_index[[j]]])
          if(length(cluster_index[[j]]) <= 1)
          {
            cov_temp[j,j] <- 0.01
          }
          else
          {
            cov_temp[j,j] <- var(X[i, cluster_index[[j]]])          
          }
          
        }
        cov[[i]] <- cov_temp
      }
      cluster <- kl_cluster(mu, cov, k)
      return(cluster)
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
          if(length(cluster_index[[j]]) <= 1)
          {
            cov_temp[j,j] <- 0.01
          }
          else
          {
            cov_temp[j,j] <- var(X[cluster_index[[j]], i])
          }
        }
        cov[[i]] <- cov_temp
      }
      cluster <- kl_cluster(mu, cov, k)
      return(cluster) 
    }
    else
    {
      print("flag undefined!")
      exit()
    }
  }
  
  else if (type == 'mean')
  {
    ## n is the number of clusters in row direction
    ## mu is a m x n matrix where each row is a n-dimensional mean vector of a distribution
    ## cov is a list of m n x n matrices where each matrix is a covariance matrix 
    if(flag == "row")
    {
      m <- nrow(X)
      n <- length(cluster_index)
      mu <- matrix(data = NA, nrow = m, ncol = n)
      for(i in 1:m)
      {
        for(j in 1:n)
        {
          mu[i,j] <- mean(X[i, cluster_index[[j]]])
        }
      }
      row_cluster <- stats::kmeans(mu, centers = k)
      return(row_cluster$cluster)
    }
    
    else if (flag == "col")
    {
      m <- length(cluster_index)
      n <- ncol(X)
      mu <- matrix(data = 0, nrow = n, ncol = m)
      for(i in 1:n)
      {
        for(j in 1:m)
        {
          mu[i,j] <- mean(X[cluster_index[[j]], i])
        }
      }
      col_cluster <- stats::kmeans(mu, centers = k)
      return(col_cluster$cluster)
    }
    else
    {
      print("flag undefined!")
      exit()
    }
  }
}

#' Apply k means clustering based biclustering algorithm to detect the row and column partitions
#' @param X noised matrix
#' @param X_baseline the ground truth to compute the clustering accuracy
#' @param type the clusterin algorithms used, 'KL' represents the KL divergence, 'mean_vecotr' represents the traditional k means
#' @param row_k the number of clusters in row direction
#' @param col_k the number of clusters in column direction
#' @return clu_row: cluster of row vectors
#' @return clu_col: cluster of column vectors
#' @return RI: rand index to measure the clustering accuracy
#' @return ARI: adjusted rand index to measure the clustering accuracy
#' @return VI: variance of information to measure the clustering accuracy
kmeans_multiple_bicluster <- function(X, X_baseline, row_k = 5, col_k = 5, type = 'KL')
{
  ## Compute the initial clusters
  col_cluster <- stats::kmeans(x = t(X), centers = col_k)
  col_cluster <- col_cluster$cluster
  col_cluster_record <- seq(from = 0, to = 0, length.out = ncol(X))
  row_cluster <- kmeans_multiple_cluster(X = X, partition = col_cluster, k = row_k, flag = "row", type = type)
  row_cluster_record <- seq(from = 0, to = 0, length.out = nrow(X))
  
  ## Iterate to compute the row cluster and column clusters until converge
  if ((typeof(all.equal(col_cluster_record, col_cluster)) == "character") || (typeof(all.equal(row_cluster_record, row_cluster)) == "character"))
  {
    col_cluster_record <- col_cluster
    row_cluster_record <- row_cluster
    col_cluster <- kmeans_multiple_cluster(X = X, partition = row_cluster, flag = "col", type = type)
    row_cluster <- kmeans_multiple_cluster(X = X, partition = col_cluster, flag = "row", type = type)
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
  ARI <- adj.rand.index(X_baseline, X_clustered)
  VI <- mcclust::vi.dist(X_baseline, X_clustered)
  return(list("row_cluster" = row_cluster, "col_cluster" = col_cluster, "RI" = RI, "ARI" = ARI, "VI" = VI))
}